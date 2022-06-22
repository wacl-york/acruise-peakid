test_that("Works with arguments 5, 3, 5", {
    bg <- identify_background(df_co2$CO2_ppm,
                              bg_sd_window=100,
                              bg_sd_threshold=0.5,
                              bg_mean_window=600)
    output_r <- detect_plumes(df_co2$CO2_ppm,
                              bg,
                              df_co2$UTC_time,
                              plume_sd_threshold=5,
                              plume_sd_starting=3,
                              plume_buffer=5)
    reticulate::py_run_string("bg = peakid.identify_background(df['conc'],
                                                               bg_sd_window=100,
                                                               bg_sd_threshold=0.5,
                                                               bg_mean_window=600)")
    reticulate::py_run_string("output_py = peakid.detect_plumes(df['conc'],
                                                                bg,
                                                                plume_sd_threshold=5,
                                                                plume_sd_starting=3,
                                                                plume_buffer=5)")
    # To compare datetimes need to explicitly parse Python as UTC.
    # Otherwise it will get implicitly parsed without a TZ
    output_py_parsed <- reticulate::py$output_py
    output_py_parsed$start <- lubridate::as_datetime(output_py_parsed$start, tz="UTC")
    output_py_parsed$end <- lubridate::as_datetime(output_py_parsed$end, tz="UTC")
    expect_equal(nrow(output_r), nrow(output_py_parsed))
    expect_equal(output_r, output_py_parsed, ignore_attr=TRUE)
})
