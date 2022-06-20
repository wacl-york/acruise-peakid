identify_background <- function(concentration,
                                bg_sd_window=180,
                                bg_sd_threshold=0.5,
                                bg_mean_window=660) {
    # Smooth concentration to identify background values (low SD threshold)
    df <- dplyr::tibble(raw = concentration) |>
        tidyr::fill(raw, .direction="down")

    bg <- RcppRoll::roll_sd(df$raw, n=bg_sd_window, align="center", fill=NA)
    # Pandas rolling functions adds extra NA to the start, R adds to the end
    # Make consistent
    if (bg_sd_window %% 2 == 0) {
        bg <- c(NA, bg[1:length(bg)-1])
    }
    # Average
    bg <- RcppRoll::roll_mean(bg, n=bg_sd_window, align="center", fill=NA)
    if (bg_sd_window %% 2 == 0) {
        bg <- c(NA, bg[1:length(bg)-1])
    }

    is_bg <- !is.na(concentration) & !is.na(bg) & bg <= bg_sd_threshold

    # Now take the original raw concentration and linearly interpolate any non-background
    # values, and then take a final smooth

    output <- rep(NA, length(concentration))
    output[is_bg] <- concentration[is_bg]
    # If want to interpolate everything, so including the NAs introduced by the
    # rolling functions then use rule=2
    output[!is_bg] <- approx(output, xout=which(!is_bg), rule=1)$y
    RcppRoll::roll_mean(output, n=bg_mean_window, align="right", fill=NA)


    # TODO look at edges!

    #roll_std = (
    #    conc.fillna(method="ffill")
    #    .rolling(bg_sd_window, center=True)
    #    .std()
    #    .rolling(bg_sd_window, center=True)
    #    .mean()
    #)
    ## Define background as when rolling std deviation is below a threshold
    #is_bg = (~conc.isna()) & (roll_std <= bg_sd_threshold)
    ## Interpolate background values when don't have background
    #bg = (
    #    conc.loc[is_bg].reindex(conc.index).interpolate().rolling(bg_mean_window).mean()
    #)
}
