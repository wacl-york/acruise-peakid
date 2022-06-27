test_that("Works with arguments 5, 3, 5", {
    bg <- identify_background(df_co2$CO2_ppm,
        bg_sd_window = 100,
        bg_sd_threshold = 0.5,
        bg_mean_window = 600
    )
    plumes <- detect_plumes(df_co2$CO2_ppm,
        bg,
        df_co2$time_nano,
        plume_sd_threshold = 5,
        plume_sd_starting = 3,
        plume_buffer = 5
    )
    areas_r <- integrate_aup_trapz(df_co2$CO2_ppm,
        df_co2$time_nano,
        plumes,
        dx = 0.1
    )
    reticulate::py_run_string("bg = peakid.identify_background(df['conc'],
                                                               bg_sd_window=100,
                                                               bg_sd_threshold=0.5,
                                                               bg_mean_window=600)")
    reticulate::py_run_string("plumes = peakid.detect_plumes(df['conc'],
                                                                bg,
                                                                plume_sd_threshold=5,
                                                                plume_sd_starting=3,
                                                                plume_buffer=5)")
    reticulate::py_run_string("areas_py = peakid.integrate_aup_trapz(df['conc'],
                                                                      plumes,
                                                                      dx=0.1)")
    reticulate::py_run_string("areas_py.reset_index(drop=True, inplace=True)")
    # To compare datetimes need to explicitly parse Python as int64 from nanotime package
    output_py_parsed <- reticulate::py$areas_py
    output_py_parsed$start <- nanotime::as.nanotime(output_py_parsed$start)
    output_py_parsed$end <- nanotime::as.nanotime(output_py_parsed$end)
    expect_equal(nrow(areas_r), nrow(output_py_parsed))
    expect_equal(areas_r, output_py_parsed, ignore_attr = TRUE)
})

test_that("Works with POSIX timestamps with implied timezone", {
    df_co2[, time_posix := as.POSIXct(UTC_time,
        format = "%d/%m/%Y %H:%M:%OS"
    )]
    bg <- identify_background(df_co2$CO2_ppm,
        bg_sd_window = 100,
        bg_sd_threshold = 0.5,
        bg_mean_window = 600
    )
    plumes <- detect_plumes(df_co2$CO2_ppm,
        bg,
        df_co2$time_posix,
        plume_sd_threshold = 5,
        plume_sd_starting = 3,
        plume_buffer = 5
    )
    areas_r <- integrate_aup_trapz(df_co2$CO2_ppm,
        df_co2$time_posix,
        plumes,
        dx = 0.1
    )
    reticulate::py_run_string("bg = peakid.identify_background(df['conc'],
                                                               bg_sd_window=100,
                                                               bg_sd_threshold=0.5,
                                                               bg_mean_window=600)")
    reticulate::py_run_string("plumes = peakid.detect_plumes(df['conc'],
                                                                bg,
                                                                plume_sd_threshold=5,
                                                                plume_sd_starting=3,
                                                                plume_buffer=5)")
    reticulate::py_run_string("areas_py = peakid.integrate_aup_trapz(df['conc'],
                                                                      plumes,
                                                                      dx=0.1)")
    reticulate::py_run_string("areas_py.reset_index(drop=True, inplace=True)")
    # To compare datetimes need to explicitly parse Python as int64 from nanotime package
    output_py_parsed <- reticulate::py$areas_py
    output_py_parsed$start <- nanotime::as.nanotime(output_py_parsed$start)
    output_py_parsed$end <- nanotime::as.nanotime(output_py_parsed$end)
    expect_equal(nrow(areas_r), nrow(output_py_parsed))
    expect_equal(areas_r, output_py_parsed, ignore_attr = TRUE)
})

test_that("Works with POSIX timestamps with explicit timezone", {
    df_co2[, time_posix_utc := as.POSIXct(UTC_time,
        format = "%d/%m/%Y %H:%M:%OS",
        tz = "UTC"
    )]
    bg <- identify_background(df_co2$CO2_ppm,
        bg_sd_window = 100,
        bg_sd_threshold = 0.5,
        bg_mean_window = 600
    )
    plumes <- detect_plumes(df_co2$CO2_ppm,
        bg,
        df_co2$time_posix_utc,
        plume_sd_threshold = 5,
        plume_sd_starting = 3,
        plume_buffer = 5
    )
    areas_r <- integrate_aup_trapz(df_co2$CO2_ppm,
        df_co2$time_posix_utc,
        plumes,
        dx = 0.1
    )
    reticulate::py_run_string("bg = peakid.identify_background(df['conc'],
                                                               bg_sd_window=100,
                                                               bg_sd_threshold=0.5,
                                                               bg_mean_window=600)")
    reticulate::py_run_string("plumes = peakid.detect_plumes(df['conc'],
                                                                bg,
                                                                plume_sd_threshold=5,
                                                                plume_sd_starting=3,
                                                                plume_buffer=5)")
    reticulate::py_run_string("areas_py = peakid.integrate_aup_trapz(df['conc'],
                                                                      plumes,
                                                                      dx=0.1)")
    reticulate::py_run_string("areas_py.reset_index(drop=True, inplace=True)")
    # To compare datetimes need to explicitly parse Python as int64 from nanotime package
    output_py_parsed <- reticulate::py$areas_py
    output_py_parsed$start <- nanotime::as.nanotime(output_py_parsed$start)
    output_py_parsed$end <- nanotime::as.nanotime(output_py_parsed$end)
    expect_equal(nrow(areas_r), nrow(output_py_parsed))
    expect_equal(areas_r, output_py_parsed, ignore_attr = TRUE)
})
