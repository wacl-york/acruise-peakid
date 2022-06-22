identify_background <- function(concentration,
                                bg_sd_window = 180,
                                bg_sd_threshold = 0.5,
                                bg_mean_window = 660) {
    # Smooth concentration to identify background values (low SD threshold)
    # TODO could do this in data.table too?!
    df <- dplyr::tibble(raw = concentration) |>
        tidyr::fill(raw, .direction = "down")

    bg <- RcppRoll::roll_sd(df$raw, n = bg_sd_window, align = "center", fill = NA)
    # Pandas rolling functions adds extra NA to the start, R adds to the end
    # Make consistent
    if (bg_sd_window %% 2 == 0) {
        bg <- c(NA, bg[1:length(bg) - 1])
    }
    # Average
    bg <- RcppRoll::roll_mean(bg, n = bg_sd_window, align = "center", fill = NA)
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
    RcppRoll::roll_mean(output, n = bg_mean_window, align = "right", fill = NA)
}

#' @import data.table
detect_plumes <- function(conc,
                          bg,
                          time,
                          plume_sd_threshold = 3,
                          plume_sd_starting = 2,
                          plume_buffer = 10) {
    dt <- data.table::data.table(time = time, conc = conc, bg = bg)
    dt[, is_plume_starting := !is.na(conc) & !is.na(bg) & conc > (bg + plume_sd_starting * sd(bg, na.rm = T))]
    dt[, is_plume := !is.na(conc) & !is.na(bg) & conc > (bg + plume_sd_threshold * sd(bg, na.rm = T))]
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

#' @import data.table
integrate_aup_trapz <- function(conc, times, plumes, dx = 1) {
    # Join plumes into the main concentration so can calculate area by plume
    dt <- data.table(conc = conc, time = times)
    plumes_dt <- data.table(plumes)
    plumes_dt[, plume_id := 1:nrow(plumes_dt)]
    areas <- dt[plumes_dt, on = c("time >= start", "time <= end")][, .(area = pracma::trapz(seq(length.out = .N, by = dx), conc)), by = list(time, time.1, plume_id)]

    areas[, plume_id := NULL]
    setnames(areas, c("time", "time.1"), c("start", "end"))
    setcolorder(areas, c("start", "end", "area"))
    areas
}
