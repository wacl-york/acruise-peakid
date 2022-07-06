test_that("Works with arguments 10, 5, 5", {
    bg <- identify_background(df_sim$conc,
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
    output_py_parsed <- reticulate::py$output_py
    # To compare datetimes need to explicitly parse Python as int64 from nanotime package
    output_py_parsed$start <- nanotime::nanotime(output_py_parsed$start)
    output_py_parsed$end <- nanotime::nanotime(output_py_parsed$end)
    expect_equal(nrow(output_r), nrow(output_py_parsed))
    expect_equal(output_r, output_py_parsed, ignore_attr = TRUE)
})

test_that("Works with nanotime", {
    # Change sample rate to be 10Hz instead so use nanotime for higher resolutionjk
    start <- nanotime::nanotime("2020-05-05 15:03:49.200", tz = "UTC")
    additions <- seq(nrow(df_sim))
    time_nano <- start + nanotime::nanoduration(hours = 0, minutes = 0, seconds = 0, nanoseconds = additions * 1e8)
    bg <- identify_background(df_sim$conc,
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
    reticulate::py_run_string("df['time_nano'] = pd.date_range(start='2020-05-05 15:03:49.200', tz='UTC', freq='100L', periods=df.shape[0])")
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
    output_py_parsed <- reticulate::py$output_py
    # To compare datetimes need to explicitly parse Python as int64 from nanotime package
    output_py_parsed$start <- nanotime::nanotime(output_py_parsed$start)
    output_py_parsed$end <- nanotime::nanotime(output_py_parsed$end)
    expect_equal(nrow(output_r), nrow(output_py_parsed))
    expect_equal(output_r, output_py_parsed, ignore_attr = TRUE)
})
