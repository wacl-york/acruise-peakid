import datetime
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.dates as mdates
mpl.rcParams['figure.figsize'] = (14.0, 8.0)

NANO_SECS_TO_SECS = 1e9

def datetime_to_seconds(dt):
    return dt.view(int) / NANO_SECS_TO_SECS

def detect_peaks(concentration, times, std_low_threshold=0.5, plume_cutoff=3, extended_plume_cutoff=2, buffer=5, rolling_period=180, group_difference=2, plume_difference=1, 
                 background_window=660, bg_subtraction_window=10, group_difference_2=1):
    df = pd.DataFrame({'conc': concentration}, index=times)
    
    buffer_time = datetime.timedelta(seconds=buffer)
    
    std = df.conc.fillna(method='ffill').rolling(rolling_period, center=True).std().rolling(rolling_period, center=True).mean()
    std.loc[df.conc.isna()] = np.nan
    std_low = std.copy()
    std_low.loc[std_low > std_low_threshold] = np.nan
    
    #derive background stats
    # TODO Does this need to be a DataFrame? Doesn't seem to actually use std_low,
    # just the times
    # TODO What is this line doing?
    # Finds groups that are separated by 2s difference when dropping NAs
    # And groups by them
    # Uses this for 2 reasons:
    #    - Find the times of the longest group
    #    - Something else...
    
    # Need to keep:
    #   - timestamps
    #   - seconds
    #   - group
    blah = pd.DataFrame({'dt': std_low.dropna().index.values})
    blah['seconds'] = datetime_to_seconds(blah['dt'])
    blah['group'] = (blah['seconds'].diff().abs() > group_difference).cumsum()
    group_counts = blah.groupby('group').group.transform('count')
    longest_group_times = blah.loc[group_counts.eq(group_counts.max()), "dt"]
    bg = df.conc.loc[blah['dt']].reindex(df.index).interpolate().rolling(background_window).mean()
    
    
    # TODO original
    #s = pd.DataFrame(std_low.dropna())
    #s['i'] = datetime_to_seconds(s.index.values)
    #groups = s.groupby((s.i.diff().abs() > group_difference).cumsum())
    
    #print(f"Found {len(groups)} groups")
    
    # Isn't this just extracting the unique datetimes?
    #group_times = [i[1].index for i in groups]
    #idx = group_times[0]
    #for gt in groups.apply(lambda x: x.index):
    #    idx = idx.union(gt)
    #bg = df.conc.loc[idx].reindex(df.index).interpolate().rolling(background_window).mean()
    
    bg_mean = df.conc.loc[longest_group_times].mean()
    bg_std = df.conc.loc[longest_group_times].std()
    
    print(f'Background mean: {bg_mean}')
    print(f'Background std: {bg_std}')
    #first iteration of plume
    N_STDS = 3
    
    foo = df.copy()
    foo['conc'] = (foo['conc']-bg).rolling(bg_subtraction_window).mean()
    foo['conc'].loc[foo['conc'] < plume_cutoff*bg_std] = np.nan
    foo['seconds'] = datetime_to_seconds(foo.index.values)
    foo['group'] = (foo.dropna().seconds.diff().abs() > group_difference_2).cumsum()
    # Get start and end time by group
    groups = (
        foo
           .reset_index()
           .groupby('group')
           .agg(
               start_time = pd.NamedAgg(column="index", aggfunc="min"),
               end_time = pd.NamedAgg(column="index", aggfunc="max")
           )
    )
    # Just want start and end times now
    
    #extend plumes
    plumes = []
    for i in range(groups.shape[0]):
        # Get the last time before the group start time where the value is less than the BG
        start = (df
           .loc[(df.index < groups.iloc[i].loc['start_time']) & (df['conc'] <= (bg_mean + extended_plume_cutoff * bg_std))]
           .tail(1)
           .index
        )[0]
        plumes.append(df['conc'].loc[(df.index>=start) & (df.index <= groups.iloc[i].loc['end_time'])])
    
    #add buffer
    expanded_plumes = []
    for plume in plumes:
        expanded_plumes.append(
            df.conc.loc[(df.index >= plume.index[0] - buffer_time) & (df.index <= plume.index[-1] + buffer_time)]
        )
    
    #merge plumes with overlapping start and end points
        
        
    # TODO test this and if works remove these comments
    plume_times = pd.concat(expanded_plumes).index.unique()
    plumes = df.copy()
    plumes.loc[~plumes.index.isin(plume_times)] = np.nan
    
    #index = expanded_plumes[0].index
    #for plume in expanded_plumes[1:]:
    #    index = index.union(plume.index)
    #plumes = df.conc.copy() * np.nan
    #plumes.loc[unique_times] = df.conc[unique_times]
    #plumes = pd.DataFrame(plumes)
    
    # TODO warning
    #plumes['i'] = plumes.index.astype(int) / NANO_SECS_TO_SECS
    plumes['seconds'] = datetime_to_seconds(plumes.index.values)
    groups = plumes.groupby((plumes.dropna().seconds.diff().abs() > plume_difference).cumsum())
    
    # TODO This is a big plot!
    # This is probably the final plot
    fig, ax = plt.subplots()
    myFmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    # TODO get label from args
    ax.set_ylabel('Concentration')
    ax.set_xlabel('Time / UTC')
    ax.plot(df.conc, color='gray', alpha=.5)
    for group in groups:
        ax.plot(group[1].conc)
    for i, dff in groups:
        bg_removed = df.conc - bg[df.index]
        print(f'plume {int(i)+1}: {dff.index[0]} - {dff.index[-1]} -> {np.trapz(bg_removed, dx=0.1)}')    
    print(bg_removed)
    plt.show()
    
    # TODO what to return?
    return plumes
    
    # TODO:
    #  - How many of these plots are useful? I.e. how many functions can I refactor
    #  it into?
    #  - Detect background
    #  - Detect plumes
    #  - Integrate area under plumes?
    
