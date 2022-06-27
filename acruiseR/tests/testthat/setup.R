# Load Python dataset + acruisepy code but keep in Python environment rather
# than attempt to convert into R. This is because the code makes use of pandas
# Indexes, which get stripped out when converting into R dataframes
reticulate::source_python("load_test_data_pandas.py", convert = FALSE)

df_co2 <- data.table::fread("../../../data/C258_FGGA_FAAM.txt")
df_co2[, time_nano := nanotime::nanotime(UTC_time,
    format = "%d/%m/%Y %H:%M:%E1S"
)]
study_period <- nanotime::as.nanoival(paste("+2021-10-04 09:40:00UTC",
    "2021-10-04 13:30:00UTC+",
    sep = "->"
))
df_co2 <- df_co2[nanotime::`%in%`(time_nano, study_period)]
