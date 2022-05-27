# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 17:34:49 2022

@author: Jake
"""

import time
import datetime
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.dates as mdates
mpl.rcParams['figure.figsize'] = (14.0, 8.0)
columns = ["conc", "a", "b"]

# Load data
data = pd.read_csv('C258_FGGA_FAAM.txt', names=columns, parse_dates=True, header=1, date_parser=lambda x:datetime.datetime.strptime(x, '%d/%m/%Y %H:%M:%S.%f'))

t_0 = time.time()

#background
ROLLING_PERIOD = 180

start_time = datetime.datetime(2021,10,4,9,40,0)
end_time = datetime.datetime(2021,10,4,13,30,0)
time_filter = (data.index >= start_time) & (data.index <= end_time)


std = data.conc.fillna(method='ffill').rolling(ROLLING_PERIOD, center=True).std().rolling(ROLLING_PERIOD, center=True).mean().loc[time_filter]
std.loc[data.conc.isna()] = np.nan
std_low = std.copy()
# TODO Extract 0.5 as magic number
std_low.loc[std_low > .5] = np.nan

#derive background stats
s = pd.DataFrame(std_low.dropna())
# TODO magic number
s['i'] = s.index.values.astype(int) / 1e9
# TODO magic number
groups = s.groupby((s.i.diff().abs() > 2).cumsum())
longest_group = sorted(groups, key=lambda x: len(x[1]))[-1]
print(len(groups))
group_times = [i[1].index for i in groups]
idx = group_times[0]
for gt in groups.apply(lambda x: x.index):
    idx = idx.union(gt)

# TODO magic number
bg = data.conc.loc[idx].reindex(data.index).interpolate().rolling(660).mean()
bg_mean = data.conc.loc[longest_group[1].index].mean()
bg_std = data.conc.loc[longest_group[1].index].std()

print(f'Background mean: {bg_mean}')
print(f'Background std: {bg_std}')
#first iteration of plume
N_STDS = 3

# TODO 10 Magic number ALSO THIS IS DIFFERENT IN SO2! DON'T SUBTRACT BG THERE!
conc = (data.conc.loc[time_filter]-bg[time_filter]).rolling(10).mean()
conc_plume = conc.copy()
conc_plume[conc_plume < N_STDS*bg_std] = np.nan

#isolate each plume
conc_plume = pd.DataFrame(conc_plume)
# TODO warning
conc_plume['i'] = [int(i)/1e9 for i in conc_plume.index.values]#.astype(int) / 1e9
groups = [
    i[1] for i in conc_plume.groupby((conc_plume.dropna().i.diff().abs() > 1).cumsum())
]

#extend plumes
CUTOFF_N_STD = 2

plumes = []
for grp in groups:
    start = data.conc.loc[data.index < grp.index[0]]
    start[start > (bg_mean + CUTOFF_N_STD*bg_std)] = np.nan
    start = start.dropna().tail(1)
    plumes.append(data.conc.loc[(data.index>=start.index[0]) & (data.index<=grp.index[-1])])

#add buffer
# TODO Magic number to parameter
BUFFER_TIME = datetime.timedelta(seconds=5)
expanded_plumes = []
for plume in plumes:
    expanded_plumes.append(
        data.conc.loc[(data.index >= plume.index[0] - BUFFER_TIME) & (data.index <= plume.index[-1] + BUFFER_TIME)]
    )

#merge plumes with overlapping start and end points
index = expanded_plumes[0].index
for plume in expanded_plumes[1:]:
    index = index.union(plume.index)
# TODO What does this do?!
plumes = data.conc.copy() * np.nan
plumes.loc[index] = data.conc[index]
plumes = pd.DataFrame(plumes)
# TODO warning
plumes['i'] = plumes.index.astype(int) / 1e9
groups = plumes.groupby((plumes.dropna().i.diff().abs() > 1).cumsum())

t_e = time.time()
print(f"Time taken: {t_e - t_0}s")

# TODO This is a big plot!
# This is probably the final plot
fig, ax = plt.subplots()
ax.plot(data.conc.loc[time_filter], color='gray', alpha=.5)
myFmt = mdates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(myFmt)
plt.ylabel('SO$_{2}$ conc / ppb')
plt.xlabel('Time / UTC')
plt.plot(data.conc.loc[time_filter], color='gray', alpha=.5)
for group in groups:
    plt.plot(group[1].conc)
areas = []
for i, df in groups:
    bg_removed = df.conc - bg[df.index]
    print(f'plume {int(i)+1}: {df.index[0]} - {df.index[-1]} -> {np.trapz(bg_removed, dx=0.1)}')    
    this_df = pd.DataFrame([{
        'start': df.index[0],
        'end': df.index[-1],
        'area': np.trapz(bg_removed, dx=0.1)
    }])
    areas.append(this_df)
print(bg_removed)
plt.show()
    
pd.concat(areas).to_csv("plumes_co2.csv")
