df_co2 <- readr::read_csv("../../../data/C258_FGGA_FAAM.txt", show_col_types = FALSE) |>
    dplyr::mutate(UTC_time = lubridate::as_datetime(UTC_time, format="%d/%m/%Y %H:%M:%OS")) |>
    dplyr::filter(UTC_time >= lubridate::as_datetime("2021-10-04 09:40:00"),
                  UTC_time <= lubridate::as_datetime("2021-10-04 13:30:00"))

# Load Python dataset + acruisepy code but keep in Python environment rather
# than attempt to convert into R. This is because the code makes use of pandas
# Indexes, which get stripped out when converting into R dataframes
reticulate::source_python("../load_test_data_pandas.py", convert=FALSE)

test_that("Works with arguments 100, 0.5, 600", {
    output_r <- identify_background(df_co2$CO2_ppm,
                                    bg_sd_window=100,
                                    bg_sd_threshold=0.5,
                                    bg_mean_window=600)
    # Run Python equivalent, again in Python environment so as to not
    # try and force pd.Dataframe -> R dataframe conversion
    reticulate::py_run_string("output_py = peakid.identify_background(df['conc'],
                                                          bg_sd_window=100,
                                                          bg_sd_threshold=0.5,
                                                          bg_mean_window=600)")
    expect_length(output_r, nrow(df_co2))
    expect_equal(length(output_r), length(reticulate::py$output_py))
    expect_equal(output_r, reticulate::py$output_py, ignore_attr=TRUE,
                 tolerance=1e1)  # Reduce tolerance from default as have extra wrapper to deal with
})
