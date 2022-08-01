#' Identifies background from a concentration time-series.
#'
#' This process takes 3 steps:
#' 1. Obtaining a smoothed rolling standard deviation of the background
#' 2. Identifying background measurements as those that lie within a set
#' threshold
#' 3. Interpolating values outside of this background threshold in a
#' linear manner and then using a rolling mean so that 'background'
#' measurements are available for the entire time-series.
#'
#' @param concentration Concentration time-series as a vector.
#' @param method The background identification method to be employed.
#' \itemize{
#'     \item{\code{gam}: Fits a Generalised Additive Model (gam), which is a non-linear
#'     regression technique. Requires parameter \code{k} to be tuned.}
#'     \item{\code{rolling}: The original 2-step rolling average method, which
#'     firstly fits a rolling standard deviation and then a rolling mean to extract
#'     the background. Requires parameters \code{bg_sd_window}, \code{bg_sd_threshold},
#'     and \code{bg_mean_window} to be tuned.}
#' }
#' @param k The number of spline points in the GAM, a higher number
#'       should provide a better fit but at the cost of overfitting to
#'       noise. Used only
#'       when \code{method='gam'}.
#' @param bg_sd_window Window size for the rolling standard deviation
#'       smooth to identify the background, as an integer. Used only
#'       when \code{method='rolling'}.
#' @param bg_sd_threshold Background measurements are considered as
#'       those whose rolling sd is within this threshold
#'       when \code{method='rolling'}. Used only
#' @param bg_mean_window The rolling mean to smooth the interpolated
#'       background, as an integer. Used only
#'       when \code{method='rolling'}.
#'
#' @return A vector with the same length as `concentration` containing the
#' smoothed background.
#' @export
identify_background <- function(concentration,
                                method = c("gam", "rolling"),
                                k = 10,
                                bg_sd_window = 180,
                                bg_sd_threshold = 0.5,
                                bg_mean_window = 660) {
    method <- match.arg(method)
    if (method == "rolling") {
        # Smooth concentration to identify background values (low SD threshold)
        bg <- RcppRoll::roll_sd(data.table::nafill(concentration, type = "locf"),
            n = bg_sd_window,
            align = "center",
            fill = NA
        )
        # With even windows there is a non-even number of NAs around the output
        # Pandas the extra NA to the start, R places it at the end
        # This line make R consistent with pandas
        if (bg_sd_window %% 2 == 0) {
            bg <- c(NA, bg[1:length(bg) - 1])
        }
        # Average
        bg <- RcppRoll::roll_mean(bg, n = bg_sd_window, align = "center", fill = NA)
        # Again move surplus NA to start
        if (bg_sd_window %% 2 == 0) {
            bg <- c(NA, bg[1:length(bg) - 1])
        }

        is_bg <- !is.na(concentration) & !is.na(bg) & bg <= bg_sd_threshold

        # Now take the original raw concentration and linearly interpolate any non-background
        # values, and then take a final smooth
        output <- rep(NA, length(concentration))
        output[is_bg] <- concentration[is_bg]
        # If want to interpolate everything, so including the NAs introduced by the
        # rolling functions then use rule=2
        output[!is_bg] <- approx(output, xout = which(!is_bg), rule = 1)$y
        bg_out <- RcppRoll::roll_mean(output, n = bg_mean_window, align = "right", fill = NA)
    } else if (method == "gam") {
        # Interpolate missing values, unlike the first method that uses LastOneCarryForward
        conc <- forecast::na.interp(conc)
        x <- seq_along(conc)

        # Fit baseline
        mod <- mgcv::gam(conc ~ s(x, bs = "cs", k = k), method = "REML", select = TRUE)
        bg_out <- mod$fitted.values
    }
    bg_out
}

#' Detects plumes in a concentration time series.
#'
#' @inheritParams identify_background
#' @param background The smoothed background time-series, as can be
#' obtained from `identify_background`. Must have the same
#' length as `concentration`.
#' @param time The time-series as a vector of POSIX or nanotime.
#' Must have the same length as `concentration`.
#' @param plume_sd_threshold The number of standard deviations that a sample
#' must exceed to be defined as a plume.
#' @param plume_sd_starting Once a plume has been identified due to it
#' crossing \code{plume_sd_threshold}, its duration is considered as the
#' times when it is more than this many standard deviatiations above the
#' background.
#' @param plume_buffer A buffer in seconds applied to plumes, so
#' that if they are overlapping they are merged into the same plume.
#'
#' @return A Data Frame where each row corresponds to a unique plume, whose
#' time boundaries are contained in the 2 columns: `start` and `end`.
#' @import data.table
#' @export
detect_plumes <- function(concentration,
                          background,
                          time,
                          plume_sd_threshold = 3,
                          plume_sd_starting = 2,
                          plume_buffer = 10) {
    # Convert into nanotime (64-bit int) if in POSIX (32-bit double)
    time <- posix_to_nanotime(time)
    residual_sd <- sd(concentration - background, na.rm = T)
    dt <- data.table::data.table(time = time, concentration = concentration, background = background)
    dt[, is_plume_starting := !is.na(concentration) & !is.na(background) & concentration > (background + plume_sd_starting * residual_sd)]
    dt[, is_plume := !is.na(concentration) & !is.na(background) & concentration > (background + plume_sd_threshold * residual_sd)]
    dt[, plume_group_starting := cumsum((is_plume_starting != data.table::shift(is_plume_starting, fill = F, type = "lag")))]

    plume_groups_dt <- dt[is_plume_starting == TRUE, list("has_plume" = sum(is_plume), start = min(time), end = max(time)), by = plume_group_starting][has_plume > 0]

    # Find overlapping plumes within the buffer
    setkey(plume_groups_dt, start, end)
    overlaps <- foverlaps(plume_groups_dt,
        plume_groups_dt[, .(start, end = end + nanotime::nanoduration(hours = 0, minutes = 0, seconds = plume_buffer, nanoseconds = 0))],
        type = "any",
        which = TRUE
    )

    overlaps[, c("x_prev_seen", "y_prev_seen") := list(duplicated(xid), duplicated(yid))]
    overlaps[, combined_plume := cumsum(!(x_prev_seen | y_prev_seen))]
    # Find group for each row number
    new_plumes <- unique(overlaps, by = c("xid", "combined_plume"))
    # Join back into overlapping plumes then elongate plumes with new members
    plumes_final <- plume_groups_dt[, xid := 1:nrow(plume_groups_dt)][new_plumes, on = "xid"][, list(
        start = min(start),
        end = max(end)
    ),
    by = combined_plume
    ]
    # Then add on baseline and convert back to datetime
    plumes_final[, combined_plume := NULL]
    setorder(plumes_final, start)
    as.data.frame(plumes_final)
}

#' Integrate the Area Under a Plume (aup) using a trapezoidal method.
#'
#' @inheritParams detect_plumes
#' @param plumes A Data Frame with 'start' and 'end' columns
#' containing plume boundaries, as returned by `detect_plumes`.
#' @param dx Sampling period, passed onto the dz argument of
#' np.trapz. I.e. the time between consecutive measurements.
#'
#' @return A Data Frame with one row per plume and 3 columns `start`, `end`, and
#' `area`. The first 2 are the same as in the input `plumes`, while `area`
#' contains the integrated area.
#' @import data.table
#' @export
integrate_aup_trapz <- function(concentration, time, plumes, dx = 1) {
    # Ensure both time columns are in the same format. Could stick with both in POSIX
    # but might as we use nanotime
    time <- posix_to_nanotime(time)
    plumes$start <- posix_to_nanotime(plumes$start)
    plumes$end <- posix_to_nanotime(plumes$end)
    # Join plumes into the main concentration so can calculate area by plume
    dt <- data.table(concentration = concentration, time = time)
    plumes_dt <- data.table(plumes)
    plumes_dt[, plume_id := 1:nrow(plumes_dt)]
    areas <- dt[plumes_dt, on = c("time >= start", "time <= end")][, .(area = pracma::trapz(seq(length.out = .N, by = dx), concentration)), by = list(time, time.1, plume_id)]

    areas[, plume_id := NULL]
    setnames(areas, c("time", "time.1"), c("start", "end"))
    setcolorder(areas, c("start", "end", "area"))
    areas
}

#' Diagnostic plot to aid background extraction
#'
#' Plots the concentration time-series highlighting the extracted background
#' (red), the threshold for what is considered a plume (orange), and when
#' plumes will be determined to have started at (blue).
#'
#' @inheritParams detect_plumes
#' @param ylabel y-axis label
#' @param xlabel x-axis label
#' @param date_fmt Format of the x-axis labels, see \code{date_labels}
#' argument of \code{ggplot2::scale_x_datetime}.
#' @param bg_alpha The alpha level of the background concentration.
#' @return A \code{ggplot2} object, so it can be modified further.
#' @export
plot_background <- function(concentration,
                            time,
                            background,
                            plume_sd_threshold = 3,
                            plume_sd_starting = 2,
                            ylabel = "Concentration",
                            xlabel = "Time (UTC)",
                            date_fmt = "%H:%M",
                            bg_alpha = 0.5) {
    time <- nanotime_to_posix(time) # Can't plot nanotime
    dt <- data.table(concentration = concentration, time = time, bg = background)
    sd_residual <- sd(concentration - background, na.rm = T)
    dt[, bg_starting := bg + plume_sd_starting * sd_residual]
    dt[, bg_threshold := bg + plume_sd_threshold * sd_residual]
    dt <- melt(dt, id.vars = "time")
    dt[, variable := factor(variable,
        levels = c("concentration", "bg", "bg_starting", "bg_threshold"),
        labels = c("Concentration", "Mean background", "Plume starting point", "Plume threshold")
    )]

    ggplot2::ggplot(dt, ggplot2::aes(x = time, y = value, colour = variable, alpha = variable)) +
        ggplot2::geom_line() +
        ggplot2::labs(x = xlabel, y = ylabel) +
        ggplot2::scale_x_datetime(date_labels = date_fmt) +
        ggplot2::scale_colour_manual("", values = c("gray", "red", "steelblue", "orange")) +
        ggplot2::scale_alpha_manual("", values = c(bg_alpha, 1, 1, 1)) +
        ggplot2::theme_minimal()
}

#' Plots plumes against the background concentration.
#'
#' @inheritParams integrate_aup_trapz
#' @inheritParams plot_background
#'
#' @return A \code{ggplot2} object, so it can be modified further.
#' @export
plot_plumes <- function(concentration,
                        time,
                        plumes,
                        ylabel = "Concentration",
                        xlabel = "Time (UTC)",
                        date_fmt = "%H:%M",
                        bg_alpha = 0.5) {
    dt <- data.table(concentration = concentration, time = time)
    plumes_dt <- as.data.table(plumes)
    plumes_dt[, plume_id := 1:nrow(plumes_dt)]
    dt <- plumes_dt[dt, on = c("start <= time", "end >= time"), .(time = as.POSIXct(time), concentration, plume_id)]
    dt[, time := nanotime_to_posix(time)] # Can't plot nanotime

    # Manually repeat high contrast colour palette to match number of plumes,
    # as palette only has 9 values. This mimics matplotlib's behaviour
    n_plumes <- nrow(plumes_dt)
    max_set_1 <- 7 # The last colour in the palette is grey, which clashes with background. Remove it.
    palette <- "Dark2"
    if (n_plumes <= max_set_1) {
        colours <- RColorBrewer::brewer.pal(n_plumes, palette)
    } else {
        n_repeats <- floor(n_plumes / max_set_1)
        mod <- max(n_plumes %% max_set_1, 3) # Can't request < 3 colours
        colours <- c(
            rep(RColorBrewer::brewer.pal(max_set_1, palette), times = n_repeats),
            RColorBrewer::brewer.pal(mod, palette)
        )
    }

    ggplot2::ggplot(dt, ggplot2::aes(x = time, y = concentration, colour = as.factor(plume_id), alpha = as.factor(is.na(plume_id)))) +
        ggplot2::geom_line() +
        ggplot2::labs(x = xlabel, y = ylabel) +
        ggplot2::scale_x_datetime(date_labels = date_fmt) +
        ggplot2::scale_colour_manual("", values = colours) +
        ggplot2::scale_alpha_manual("", values = c(1, bg_alpha)) +
        ggplot2::guides(colour = "none", alpha = "none") +
        ggplot2::theme_minimal()
}

posix_to_nanotime <- function(x) {
    if (any(grepl("POSIX", class(x)))) {
        tz <- attr(x, "tz")
        if (tz == "") {
            tz <- "UTC"
        }
        x <- nanotime::as.nanotime(strftime(x, "%Y-%m-%d %H:%M:%OS3"), format = "%Y-%m-%d %H:%M:%ES", tz = tz)
    }
    x
}

nanotime_to_posix <- function(x) {
    if (any(grepl("nanotime", class(x)))) {
        tz <- "UTC"
        x <- as.POSIXct(x, tz = tz)
    }
    x
}
