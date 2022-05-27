# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 09:34:21 2022

@author: Jake
"""

import time
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.dates as mdates
mpl.rcParams['figure.figsize'] = (14.0, 8.0)
myFmt = mdates.DateFormatter('%H:%M')

#Housekeeping
columns = ["conc"]
data = pd.read_csv('C258_SO2_FAAM.txt', names=columns, parse_dates=True, header=1, date_parser=lambda x:datetime.datetime.strptime(x, '%d/%m/%Y %H:%M:%S'))
#identify background
ROLLING_PERIOD = 180

t_0 = time.time()

start_time = datetime.datetime(2021,10,4,9,40,0)
end_time = datetime.datetime(2021,10,4,13,30,0)
time_filter = (data.index >= start_time) & (data.index <= end_time)

# TODO REMOVE THIS WHEN FINISHED WITH TIME FILTER
data = data.loc[time_filter]

nans = data.conc.isna()
#std = data.conc.fillna(method='ffill').rolling(ROLLING_PERIOD, center=True).std().rolling(ROLLING_PERIOD, center=True).mean().loc[time_filter]
# TODO REVERT WHEN FINISHED WITH TIME FILTER
std = data.conc.fillna(method='ffill').rolling(ROLLING_PERIOD, center=True).std().rolling(ROLLING_PERIOD, center=True).mean()
std.loc[nans] = np.nan
std_low = std.copy()
std_low.loc[std_low > .5] = np.nan

#derive background stats
s = pd.DataFrame(std_low.dropna())
s['i'] = s.index.values.astype(int) / 1e9
groups = s.groupby((s.i.diff().abs() > 2).cumsum())
print(f"# groups {len(groups)}")
longest_group = sorted(groups, key=lambda x: len(x[1]))[-1]
print(len(groups))
group_times = [i[1].index for i in groups]
idx = group_times[0]
for gt in groups.apply(lambda x: x.index):
    idx = idx.union(gt)

bg = data.conc.loc[idx].reindex(data.index).interpolate().rolling(660).mean()
bg_mean = data.conc.loc[longest_group[1].index].mean()
bg_std = data.conc.loc[longest_group[1].index].std()

print(f'Background mean: {bg_mean}')
print(f'Background std: {bg_std}')
#first iteration of plume
N_STDS = 3

# TODO Found difference between this and CO2!
# In CO2 we subtract 'bg' here
# NB In the final peak integration we also subtract BG in both CO2 and SO2,
# so why isn't it subtracted here?
#conc = data.conc.loc[time_filter].rolling(10).mean()                   # SO2 version
#conc = (data.conc.loc[time_filter]-bg[time_filter]).rolling(10).mean()  # CO2 version
# TODO REVERT TO CO2 VERSION WHEN FINISHED WITH TIME FILTER
conc = (data.conc - bg).rolling(10).mean()  # CO2 version
conc_plume = conc.copy()
# conc_plume[conc_plume < bg_mean + N_STDS*bg_std] = np.nan  # SO2 version
conc_plume[conc_plume < N_STDS*bg_std] = np.nan          # CO2 version


#isolate each bit of plume data   
conc_plume = pd.DataFrame(conc_plume)
conc_plume['i'] = [int(i)/1e9 for i in conc_plume.index.values]
groups = [
    i[1] for i in conc_plume.groupby((conc_plume.dropna().i.diff().abs() > 1).cumsum())
]
# NB: Orig gets 21, new gets 20 here
# See if applying time-filter up front makes the difference
print(f"# groups after conc_plume {len(groups)}")

#extend each plume
CUTOFF_N_STD = 2

plumes = []
for grp in groups:
    start = data.conc.loc[data.index < grp.index[0]]
    start[start > (bg_mean + CUTOFF_N_STD*bg_std)] = np.nan
    start = start.dropna().tail(1)
    plumes.append(data.conc.loc[(data.index>=start.index[0]) & (data.index<=grp.index[-1])])
print(f"# plumes {len(plumes)}")

#buffer time
BUFFER_TIME = datetime.timedelta(seconds=5)

expanded_plumes = []
for plume in plumes:
    expanded_plumes.append(
        data.conc.loc[(data.index >= plume.index[0] - BUFFER_TIME) & (data.index <= plume.index[-1] + BUFFER_TIME)]
    )
print(f"# expanded_plumes {len(expanded_plumes)}")

#merge plumes with overlapping start and end points
index = expanded_plumes[0].index
for plume in expanded_plumes[1:]:
    index = index.union(plume.index)
plumes = data.conc.copy() * np.nan
plumes.loc[index] = data.conc[index]
plumes = pd.DataFrame(plumes)
plumes['i'] = plumes.index.astype(int) / 1e9
groups = plumes.groupby((plumes.dropna().i.diff().abs() > 1).cumsum())

print(f"# groups {len(groups)}")

t_e = time.time()
print(f"Time taken: {t_e - t_0}s")

fig, gx = plt.subplots()
#gx.plot(data.conc.loc[time_filter], color='gray', alpha=.5)
# TODO When finished with time filter revert this
gx.plot(data.conc, color='gray', alpha=.5)

for group in groups:
    plt.plot(group[1].conc)
gx.xaxis.set_major_formatter(myFmt)
plt.ylabel('SO$_{2}$ conc / ppb')
plt.show()

areas = []
for i, df in groups:
    bg_removed = df.conc - bg[df.index]
    print(f'plume {int(i)+1}: {df.index[0]} - {df.index[-1]} -> {np.trapz(bg_removed)}')   
    this_df = pd.DataFrame([{
        'start': df.index[0],
        'end': df.index[-1],
        'area': np.trapz(bg_removed)
    }])
    areas.append(this_df)
    
pd.concat(areas).to_csv("plumes_so2.csv")
    

### ignore: finding coordinates ###
#from netCDF4 import Dataset
#from cftime import num2date
#raise Exception
#first_name='core_faam_20211004_v005_r1_c258_1hz.nc'
#with Dataset(first_name,'r') as start_fh:
#    lat = start_fh['LAT_GIN'][:]
#    lon = start_fh['LON_GIN'][:]
#    time = num2date(start_fh['Time'][:], units=start_fh['Time'].units, only_use_cftime_datetimes=False)
#
#data = pd.DataFrame({'lon': lon, 'lat': lat}, index=time)
#plume_time = datetime.datetime(2021,10,4,12,52,30)
#print(data.lon[plume_time], data.lat[plume_time])
    
