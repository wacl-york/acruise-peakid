import datetime
import time
import pandas as pd
import detect_peaks

start_time = datetime.datetime(2021, 10, 4, 9, 40, 0)
end_time = datetime.datetime(2021, 10, 4, 13, 30, 0)

df_co2 = pd.read_csv(
    "C258_FGGA_FAAM.txt",
    names=["conc", "a", "b"],
    parse_dates=True,
    header=1,
    date_parser=lambda x: datetime.datetime.strptime(x, "%d/%m/%Y %H:%M:%S.%f"),
)
time_filter = (df_co2.index >= start_time) & (df_co2.index <= end_time)
df_co2 = df_co2.loc[time_filter]

co2_conc = df_co2.conc

t_0 = time.time()
co2_bg = detect_peaks.identify_background(co2_conc)
co2_plumes = detect_peaks.detect_plumes(co2_conc, co2_bg)
detect_peaks.plot_plumes(co2_conc, co2_plumes)
co2_aups = detect_peaks.integrate_aup_trapz(co2_conc - co2_bg, co2_plumes, dx=0.1)
t_e = time.time()
# Original script took 1.808s, immediate refactor at 0.908s, now at 0.48s
print(f"Time taken: {t_e - t_0}s")

# Now to regression test
co2_plumes_original = pd.read_csv("plumes_co2.csv")
co2_plumes.shape
co2_plumes_original.shape
# Same number of plumes
co2_plumes.shape[0] == co2_plumes_original.shape[0]
# Don't have the same start and end times, unsurprising given that I've removed the buffer
all(co2_plumes["start"].values == pd.to_datetime(co2_plumes_original["start"]).values)
all(co2_plumes["end"].values == pd.to_datetime(co2_plumes_original["end"]).values)
# Areas are different too
all((co2_aups["area"].values - co2_plumes_original["area"].values) < 1e-9)
# Mean area difference is -9
(co2_aups["area"].values - co2_plumes_original["area"].values).mean()

# Does this work on SO2?
df_so2 = pd.read_csv(
    "C258_SO2_FAAM.txt",
    names=["conc"],
    parse_dates=True,
    header=1,
    date_parser=lambda x: datetime.datetime.strptime(x, "%d/%m/%Y %H:%M:%S"),
)
time_filter = (df_so2.index >= start_time) & (df_so2.index <= end_time)
df_so2 = df_so2.loc[time_filter]

so2_conc = df_so2.conc

t_0 = time.time()
so2_bg = detect_peaks.identify_background(so2_conc)
so2_plumes = detect_peaks.detect_plumes(so2_conc, so2_bg)
detect_peaks.plot_plumes(so2_conc, so2_plumes)
so2_aups = detect_peaks.integrate_aup_trapz(so2_conc - so2_bg, so2_plumes, dx=0.1)
t_e = time.time()
# Original script took 0.615s, then 0.313s, now 0.26s
print(f"Time taken: {t_e - t_0}s")

so2_plumes_original = pd.read_csv("plumes_so2.csv")
so2_plumes.shape
so2_plumes_original.shape
# Found more plumes, unsurprising given that it has a more strict buffer so harder for plumes to be overlapping
so2_plumes.shape[0] == so2_plumes_original.shape[0]
# Have the same start and end times
all(so2_plumes["start"].values == pd.to_datetime(so2_plumes_original["start"]).values)
all(so2_plumes["end"].values == pd.to_datetime(so2_plumes_original["end"]).values)
# And don't have the same areas
all((so2_aups["area"].values - so2_plumes_original["area"].values) < 1e-9)
