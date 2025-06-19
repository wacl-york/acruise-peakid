test_that("Works with arguments 100, 1, 600", {
    output_r <- identify_background(df_sim$conc,
        method="rolling",
        bg_sd_window = 100,
        bg_sd_threshold = 1,
        bg_mean_window = 600
    )

    reticulate::py_run_string("from acruisepy import peakid")
    reticulate::py_run_string("df = r.df_sim.set_index('time')")
    reticulate::py_run_string("output_py = peakid.identify_background(df['conc'],
                                                                      bg_sd_window=100,
                                                                      bg_sd_threshold=1,
                                                                      bg_mean_window=600)")
    # Parse Python array into same format as R vector
    bg_py <- as.numeric(reticulate::py$output_py)
    bg_py[is.nan(bg_py)] <- NA

    expect_length(output_r$bg, nrow(df_sim))
    expect_equal(length(output_r$bg), length(bg_py))
    expect_equal(output_r$bg, bg_py)
})
