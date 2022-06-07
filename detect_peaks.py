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

def detect_peaks(concentration,
                 times,
                 std_low_threshold=0.5,    # Background std deviation
                 plume_cutoff=3,           # Plume cutoff
                 buffer=5,                 # Overlapping plumes
                 rolling_period=180,       # Time to smooth background
                 group_difference=1,       # Could probably work with 1s, shouldn't needed to be changed unless using different frequency
                 background_window=660,    # Could try with different, could be lower
                 bg_subtraction_window=10, # Trying to remove bit of noise
                 trapz_dx=1.0):
                     
    buffer_time = datetime.timedelta(seconds=buffer)
    df = pd.DataFrame({'conc': concentration}, index=times)
    
    # Define background as when rolling std deviation is below a threshold
    roll_std = df.conc.fillna(method='ffill').rolling(rolling_period, center=True).std().rolling(rolling_period, center=True).mean()
    is_bg = (~df['conc'].isna()) & (roll_std <= std_low_threshold)
    
    bg = df.conc.iloc[is_bg.values].reindex(df.index).interpolate().rolling(background_window).mean()
    bg_mean = df.conc.iloc[is_bg.values].mean()
    bg_std = df.conc.iloc[is_bg.values].std()
    
    # first iteration of plume detection
    conc_bg_removed = (df['conc'] - bg).rolling(bg_subtraction_window).mean()
    is_plume = conc_bg_removed >= plume_cutoff * bg_std
    plume_groups = (is_plume != is_plume.shift()).cumsum()
    plume_groups = (
        plume_groups
            .iloc[is_plume.values]
            .reset_index()
            .groupby('conc')
            .agg(
                   start = pd.NamedAgg(column="index", aggfunc="min"),
                   end = pd.NamedAgg(column="index", aggfunc="max")
            )
            .assign(
                start = lambda x: x['start'] - buffer_time,
                end = lambda x: x['end'] + buffer_time
            )
    )
    
    # Find all unique time points that are considered within a plume
    plume_times = [df.loc[start:end] for start, end in zip(plume_groups['start'], plume_groups['end'])]
    plume_times = pd.concat(plume_times).index.unique()
    
    # Link back to main data set and find groups
    plume_times = df.loc[plume_times].copy()
    plume_times['seconds'] = datetime_to_seconds(plume_times.index)
    groups = plume_times.groupby((plume_times.seconds.diff().abs() > group_difference).cumsum())
    
    fig, ax = plt.subplots()
    myFmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    # TODO get label from args
    ax.set_ylabel('Concentration')
    ax.set_xlabel('Time / UTC')
    ax.plot(df.conc, color='gray', alpha=.5)
    for group in groups:
        ax.plot(group[1].conc)
    areas = []
    for i, dff in groups:
        bg_removed = dff.conc - bg[dff.index]
        this_df = pd.DataFrame([{
            'start': dff.index[0],
            'end': dff.index[-1],
            'area': np.trapz(bg_removed, dx=trapz_dx)
        }])
        areas.append(this_df)
    plt.show()
    
    return pd.concat(areas)

# TODO:
#  - How many of these plots are useful? I.e. how many functions can I refactor
#  it into?
#  - Detect background
#  - Detect plumes
#  - Integrate area under plumes?

