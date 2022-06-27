test_that("Works with arguments 5, 3, 5", {
    bg <- identify_background(df_co2$CO2_ppm,
        bg_sd_window = 100,
        bg_sd_threshold = 0.5,
        bg_mean_window = 600
    )
    output_r <- detect_plumes(df_co2$CO2_ppm,
        bg,
        df_co2$time_nano,
        plume_sd_threshold = 5,
        plume_sd_starting = 3,
        plume_buffer = 5
    )
    reticulate::py_run_string("bg = peakid.identify_background(df['conc'],
                                                               bg_sd_window=100,
                                                               bg_sd_threshold=0.5,
                                                               bg_mean_window=600)")
    reticulate::py_run_string("output_py = peakid.detect_plumes(df['conc'],
                                                                bg,
                                                                plume_sd_threshold=5,
                                                                plume_sd_starting=3,
                                                                plume_buffer=5)")
    output_py_parsed <- reticulate::py$output_py
    # To compare datetimes need to explicitly parse Python as int64 from nanotime package
    output_py_parsed$start <- nanotime::nanotime(output_py_parsed$start)
    output_py_parsed$end <- nanotime::nanotime(output_py_parsed$end)
    expect_equal(nrow(output_r), nrow(output_py_parsed))
    expect_equal(output_r, output_py_parsed, ignore_attr = TRUE)
})

test_that("Works with POSIX time", {
    df_co2[, time_posix := as.POSIXct(UTC_time,
        format = "%d/%m/%Y %H:%M:%OS"
    )]
    bg <- identify_background(df_co2$CO2_ppm,
        bg_sd_window = 50,
        bg_sd_threshold = 0.5,
        bg_mean_window = 500
    )
    output_r <- detect_plumes(df_co2$CO2_ppm,
        bg,
        df_co2$time_posix,
        plume_sd_threshold = 6,
        plume_sd_starting = 2,
        plume_buffer = 10
    )
    reticulate::py_run_string("bg = peakid.identify_background(df['conc'],
                                                               bg_sd_window=50,
                                                               bg_sd_threshold=0.5,
                                                               bg_mean_window=500)")
    reticulate::py_run_string("output_py = peakid.detect_plumes(df['conc'],
                                                                bg,
                                                                plume_sd_threshold=6,
                                                                plume_sd_starting=2,
                                                                plume_buffer=10)")
    output_py_parsed <- reticulate::py$output_py
    # To compare datetimes need to explicitly parse Python as int64 from nanotime package
    output_py_parsed$start <- nanotime::nanotime(output_py_parsed$start)
    output_py_parsed$end <- nanotime::nanotime(output_py_parsed$end)
    expect_equal(nrow(output_r), nrow(output_py_parsed))
    expect_equal(output_r, output_py_parsed, ignore_attr = TRUE)
})

test_that("Works with POSIX time with explicit timezone", {
    df_co2[, time_posix_utc := as.POSIXct(UTC_time,
        format = "%d/%m/%Y %H:%M:%OS",
        tz = "UTC"
    )]
    bg <- identify_background(df_co2$CO2_ppm,
        bg_sd_window = 50,
        bg_sd_threshold = 0.5,
        bg_mean_window = 500
    )
    output_r <- detect_plumes(df_co2$CO2_ppm,
        bg,
        df_co2$time_posix_utc,
        plume_sd_threshold = 6,
        plume_sd_starting = 2,
        plume_buffer = 10
    )
    reticulate::py_run_string("bg = peakid.identify_background(df['conc'],
                                                               bg_sd_window=50,
                                                               bg_sd_threshold=0.5,
                                                               bg_mean_window=500)")
    reticulate::py_run_string("output_py = peakid.detect_plumes(df['conc'],
                                                                bg,
                                                                plume_sd_threshold=6,
                                                                plume_sd_starting=2,
                                                                plume_buffer=10)")
    output_py_parsed <- reticulate::py$output_py
    # To compare datetimes need to explicitly parse Python as int64 from nanotime package
    output_py_parsed$start <- nanotime::nanotime(output_py_parsed$start)
    output_py_parsed$end <- nanotime::nanotime(output_py_parsed$end)
    expect_equal(nrow(output_r), nrow(output_py_parsed))
    expect_equal(output_r, output_py_parsed, ignore_attr = TRUE)
})
