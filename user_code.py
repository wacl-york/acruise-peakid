import datetime
import time
import pandas as pd
from detect_peaks import detect_peaks

start_time = datetime.datetime(2021,10,4,9,40,0)
end_time = datetime.datetime(2021,10,4,13,30,0)

df_co2 = pd.read_csv('C258_FGGA_FAAM.txt', names=["conc", "a", "b"], parse_dates=True, header=1, date_parser=lambda x:datetime.datetime.strptime(x, '%d/%m/%Y %H:%M:%S.%f'))
time_filter = (df_co2.index >= start_time) & (df_co2.index <= end_time)
df_co2 = df_co2.loc[time_filter]

co2_conc = df_co2.conc.values
co2_times = df_co2.index.values

t_0 = time.time()
co2_plumes = detect_peaks(co2_conc, co2_times, trapz_dx=0.1)
t_e = time.time()
# Original script took 1.808s, now at 0.908s
print(f"Time taken: {t_e - t_0}s")

# Now to regression test
co2_plumes_original = pd.read_csv("plumes_co2.csv")
# Same number of plumes, ignore the difference in # columns
co2_plumes.shape
co2_plumes_original.shape
# Have the same start and end times
all(co2_plumes['start'].values == pd.to_datetime(co2_plumes_original['start']).values)
all(co2_plumes['end'].values == pd.to_datetime(co2_plumes_original['end']).values)
# And all areas are the same within float precision
all((co2_plumes['area'].values - co2_plumes_original['area'].values) < 1e9)

# Does this work on SO2?
df_so2 = pd.read_csv('C258_SO2_FAAM.txt', names=["conc"], parse_dates=True, header=1, date_parser=lambda x:datetime.datetime.strptime(x, '%d/%m/%Y %H:%M:%S'))
time_filter = (df_so2.index >= start_time) & (df_so2.index <= end_time)
df_so2 = df_so2.loc[time_filter]

so2_conc = df_so2.conc.values
so2_times = df_so2.index.values

t_0 = time.time()
so2_plumes = detect_peaks(so2_conc, so2_times)
t_e = time.time()
# Original script took 0.615s, now at 0.313s
print(f"Time taken: {t_e - t_0}s")
so2_plumes_original = pd.read_csv("plumes_so2.csv")

# Same number of plumes (caveat with when using time_filter)
so2_plumes.shape
so2_plumes_original.shape
# Have the same start and end times
all(so2_plumes['start'].values == pd.to_datetime(so2_plumes_original['start']).values)
all(so2_plumes['end'].values == pd.to_datetime(so2_plumes_original['end']).values)
## And all areas are the same within float precision
all((so2_plumes['area'].values - so2_plumes_original['area'].values) < 1e9)
