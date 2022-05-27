import datetime
import pandas as pd
from detect_peaks import detect_peaks

start_time = datetime.datetime(2021,10,4,9,40,0)
end_time = datetime.datetime(2021,10,4,13,30,0)

df_co2 = pd.read_csv('C258_FGGA_FAAM.txt', names=["conc", "a", "b"], parse_dates=True, header=1, date_parser=lambda x:datetime.datetime.strptime(x, '%d/%m/%Y %H:%M:%S.%f'))
time_filter = (df_co2.index >= start_time) & (df_co2.index <= end_time)
df_co2 = df_co2.loc[time_filter]

co2_conc = df_co2.conc.values
co2_times = df_co2.index.values
detect_peaks(co2_conc, co2_times)


# Does this work on SO2?
df_so2 = pd.read_csv('C258_SO2_FAAM.txt', names=["conc"], parse_dates=True, header=1, date_parser=lambda x:datetime.datetime.strptime(x, '%d/%m/%Y %H:%M:%S'))
time_filter = (df_so2.index >= start_time) & (df_so2.index <= end_time)
df_so2 = df_so2.loc[time_filter]

so2_conc = df_so2.conc.values
so2_times = df_so2.index.values
detect_peaks(so2_conc, so2_times)

