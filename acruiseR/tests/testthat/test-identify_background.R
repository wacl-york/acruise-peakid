test_that("Works with arguments 100, 0.5, 600", {
    output_r <- identify_background(df_co2$CO2_ppm,
                                    bg_sd_window=100,
                                    bg_sd_threshold=0.5,
                                    bg_mean_window=600)
    reticulate::py_run_string("output_py = peakid.identify_background(df['conc'],
                                                          bg_sd_window=100,
                                                          bg_sd_threshold=0.5,
                                                          bg_mean_window=600)")
    expect_length(output_r, nrow(df_co2))
    expect_equal(length(output_r), length(reticulate::py$output_py))
    expect_equal(output_r, reticulate::py$output_py, ignore_attr=TRUE)
})
