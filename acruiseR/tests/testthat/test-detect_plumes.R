test_that("Works with arguments 10, 5, 5", {
    bg <- identify_background(df_sim$conc,
        method='rolling',
        bg_sd_window = 100,
        bg_sd_threshold = 1,
        bg_mean_window = 600
    )
    output_r <- detect_plumes(df_sim$conc,
        bg,
        df_sim$time,
        plume_sd_threshold = 10,
        plume_sd_starting = 5,
        plume_buffer = 5
    )
    reticulate::py_run_string("from acruisepy import peakid")
    reticulate::py_run_string("df = r.df_sim.set_index('time')")
    reticulate::py_run_string("bg = peakid.identify_background(df['conc'],
                                                               bg_sd_window=100,
                                                               bg_sd_threshold=1,
                                                               bg_mean_window=600)")
    reticulate::py_run_string("output_py = peakid.detect_plumes(df['conc'],
                                                                bg,
                                                                plume_sd_threshold=10,
                                                                plume_sd_starting=5,
                                                                plume_buffer=5)")
    output_py_parsed <- as.data.table(reticulate::py$output_py)
    # To compare datetimes need to explicitly parse Python as int64 from nanotime package
    output_py_parsed[, time := nanotime::nanotime(rownames(reticulate::py$output_py),
                                                  format="%Y-%m-%d %H:%M:%S") ]
    output_py_wide <- output_py_parsed[, .(start=min(time), end=max(time)), by=plume_id]
    output_py_wide[, plume_id := NULL ]
    expect_equal(nrow(output_r), nrow(output_py_wide))
    expect_equal(output_r, output_py_wide, ignore_attr = TRUE)
})

test_that("Works with nanotime", {
    # Change sample rate to be 10Hz instead so use nanotime for higher resolution
    start <- nanotime::nanotime("2020-05-05 15:03:49.200", tz = "UTC")
    additions <- seq(nrow(df_sim))
    time_nano <- start + nanotime::nanoduration(hours = 0, minutes = 0, seconds = 0, nanoseconds = additions * 1e8)
    bg <- identify_background(df_sim$conc,
        method="rolling",
        bg_sd_window = 100,
        bg_sd_threshold = 1,
        bg_mean_window = 600
    )
    output_r <- detect_plumes(df_sim$conc,
        bg,
        time_nano,
        plume_sd_threshold = 10,
        plume_sd_starting = 5,
        plume_buffer = 5
    )
    reticulate::py_run_string("from acruisepy import peakid")
    reticulate::py_run_string("import pandas as pd")
    reticulate::py_run_string("df = r.df_sim")
    reticulate::py_run_string("df['time_nano'] = pd.date_range(start='2020-05-05 15:03:49.200', tz='UTC', freq='100ms', periods=df.shape[0])")
    reticulate::py_run_string("df.set_index('time_nano', inplace=True)")
    reticulate::py_run_string("bg = peakid.identify_background(df['conc'],
                                                               bg_sd_window=100,
                                                               bg_sd_threshold=1,
                                                               bg_mean_window=600)")
    reticulate::py_run_string("output_py = peakid.detect_plumes(df['conc'],
                                                                bg,
                                                                plume_sd_threshold=10,
                                                                plume_sd_starting=5,
                                                                plume_buffer=5)")
    reticulate::py_run_string("time_str = pd.Series(output_py.index).dt.strftime('%Y-%m-%d %H:%M:%S.%f').values")

    output_py_parsed <- as.data.table(reticulate::py$output_py)
    # To compare datetimes need to explicitly parse Python as int64 from nanotime package
    output_py_parsed[, time := nanotime::nanotime(as.character(reticulate::py$time_str), format="%Y-%m-%d %H:%M:%E6S")]
    output_py_wide <- output_py_parsed[, .(start=min(time), end=max(time)), by=plume_id]
    output_py_wide[, plume_id := NULL ]
    expect_equal(nrow(output_r), nrow(output_py_wide))
    expect_equal(output_r, output_py_wide, ignore_attr = TRUE)
})
