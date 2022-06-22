# Load Python dataset + acruisepy code but keep in Python environment rather
# than attempt to convert into R. This is because the code makes use of pandas
# Indexes, which get stripped out when converting into R dataframes
reticulate::source_python("load_test_data_pandas.py", convert=FALSE)

df_co2 <- readr::read_csv("../../../data/C258_FGGA_FAAM.txt", show_col_types = FALSE) |>
    dplyr::mutate(UTC_time = lubridate::as_datetime(UTC_time, format="%d/%m/%Y %H:%M:%OS")) |>
    dplyr::filter(UTC_time >= lubridate::as_datetime("2021-10-04 09:40:00"),
                  UTC_time <= lubridate::as_datetime("2021-10-04 13:30:00"))
