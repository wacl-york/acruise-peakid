# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 17:34:49 2022

@author: Jake
"""

import pandas as pd
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime
import numpy as np
import matplotlib.dates as mdates
mpl.rcParams['figure.figsize'] = (14.0, 8.0)
columns = ["conc", "a", "b"]

# Load data
data = pd.read_csv('C258_FGGA_FAAM.txt', names=columns, parse_dates=True, header=1, date_parser=lambda x:datetime.datetime.strptime(x, '%d/%m/%Y %H:%M:%S.%f'))
#plot data
plt.figure()
plt.plot(data.conc)
plt.ylabel('CO2 conc.')
#background
ROLLING_PERIOD = 180

start_time = datetime.datetime(2021,10,4,9,40,0)
end_time = datetime.datetime(2021,10,4,13,30,0)
time_filter = (data.index >= start_time) & (data.index <= end_time)
plt.figure()
plt.plot(data.conc.loc[time_filter])

nans = data.conc.isna()
std = data.conc.fillna(method='ffill').rolling(ROLLING_PERIOD, center=True).std().rolling(ROLLING_PERIOD, center=True).mean().loc[time_filter]
std.loc[nans] = np.nan
plt.plot(std)
std_low = std.copy()
std_low.loc[std_low > .5] = np.nan
plt.gca().twinx().plot(std, color='tab:green', linewidth=2)
plt.plot(std_low, 'tab:red', linewidth=4)
#derive background stats
s = pd.DataFrame(std_low.dropna())
s['i'] = s.index.values.astype(int) / 1e9
groups = s.groupby((s.i.diff().abs() > 2).cumsum())
longest_group = sorted(groups, key=lambda x: len(x[1]))[-1]
print(len(groups))
group_times = [i[1].index for i in groups]
idx = group_times[0]
for gt in groups.apply(lambda x: x.index):
    idx = idx.union(gt)
    
bg = data.conc.loc[idx].reindex(data.index).interpolate().rolling(660).mean()
# plt.plot(bg.rolling(180).mean(), 'k', linewidth=3)
plt.figure()
plt.plot(data.conc[time_filter] - bg.loc[time_filter])
for group in groups:
    plt.plot(data.conc[group[1].index] - bg.loc[group[1].index], color='tab:red')
plt.plot(data.conc.loc[longest_group[1].index] - bg.loc[longest_group[1].index], color='magenta')

group_times = [i[1].index.mean() for i in groups]
group_means = [i[1].conc.mean() for i in groups]

background = pd.Series(group_means, index=group_times).reindex(data.index).interpolate()
plt.plot(background, 'k', linewidth=3)

    
bg_mean = data.conc.loc[longest_group[1].index].mean()
bg_std = data.conc.loc[longest_group[1].index].std()

print(f'Background mean: {bg_mean}')
print(f'Background std: {bg_std}')
#first iteration of plume
N_STDS = 3

conc = (data.conc.loc[time_filter]-bg[time_filter]).rolling(10).mean()
plt.figure()
plt.plot(conc)
conc_plume = conc.copy()
conc_plume[conc_plume < N_STDS*bg_std] = np.nan
plt.plot(conc_plume)
# plt.plot([conc.index[0], conc.index[-1]], [bg_mean + N_STDS*bg_std]*2, '--')
# plt.gca().set_ylim([-1 ,10])

#isolate each plume
conc_plume = pd.DataFrame(conc_plume)
conc_plume['i'] = [int(i)/1e9 for i in conc_plume.index.values]#.astype(int) / 1e9
groups = [
    i[1] for i in conc_plume.groupby((conc_plume.dropna().i.diff().abs() > 1).cumsum())
]

plt.figure()
for grp in groups:
    plt.plot(grp.conc)
    plt.ylabel('SO2 conc / ppb')
#extend plumes
CUTOFF_N_STD = 2

plumes = []
for grp in groups:
    start = data.conc.loc[data.index < grp.index[0]]
    start[start > (bg_mean + CUTOFF_N_STD*bg_std)] = np.nan
    start = start.dropna().tail(1)
    plumes.append(data.conc.loc[(data.index>=start.index[0]) & (data.index<=grp.index[-1])])
plt.figure()
for plume in plumes:
    plt.plot(plume)
plt.plot([conc.index[0], conc.index[-1]], [bg_mean + 2*bg_std]*2, '--')
#add buffer
BUFFER_TIME = datetime.timedelta(seconds=5)

expanded_plumes = []
for plume in plumes:
    expanded_plumes.append(
        data.conc.loc[(data.index >= plume.index[0] - BUFFER_TIME) & (data.index <= plume.index[-1] + BUFFER_TIME)]
    )

plt.figure()
plt.plot(data.conc.loc[time_filter], color='gray', alpha=.5)
for plume in expanded_plumes:
    plt.plot(plume)
#merge plumes with overlapping start and end points
index = expanded_plumes[0].index
for plume in expanded_plumes[1:]:
    index = index.union(plume.index)
plumes = data.conc.copy() * np.nan
plumes.loc[index] = data.conc[index]
plumes = pd.DataFrame(plumes)
plumes['i'] = plumes.index.astype(int) / 1e9
groups = plumes.groupby((plumes.dropna().i.diff().abs() > 1).cumsum())
fig, ax = plt.subplots()
ax.plot(data.conc.loc[time_filter], color='gray', alpha=.5)
myFmt = mdates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(myFmt)
plt.ylabel('SO$_{2}$ conc / ppb')
plt.xlabel('Time / UTC')
plt.plot(data.conc.loc[time_filter], color='gray', alpha=.5)
for group in groups:
    plt.plot(group[1].conc)
for i, df in groups:
    bg_removed = df.conc - bg[df.index]
    print(f'plume {int(i)+1}: {df.index[0]} - {df.index[-1]} -> {np.trapz(bg_removed, dx=0.1)}')    
print(bg_removed)


