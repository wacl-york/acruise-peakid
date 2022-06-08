from functools import reduce
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
                 bg_sd_window=180,         # Time to smooth background
                 bg_sd_threshold=0.5,      # Background std deviation
                 bg_mean_window=660,       # Could try with different, could be lower
                 bg_subtraction_window=10, # Trying to remove bit of noise
                 plume_sd_threshold=3,     # Plume cutoff
                 plume_buffer=10,          # Overlapping plumes in seconds
                 trapz_dx=1.0):
                     
    df = pd.DataFrame({'conc': concentration}, index=times)
    
    # Define background as when rolling std deviation is below a threshold
    roll_std = df.conc.fillna(method='ffill').rolling(bg_sd_window, center=True).std().rolling(bg_sd_window, center=True).mean()
    is_bg = (~df['conc'].isna()) & (roll_std <= bg_sd_threshold)
    
    bg = df.conc.iloc[is_bg.values].reindex(df.index).interpolate().rolling(bg_mean_window).mean()
    bg_mean = df.conc.iloc[is_bg.values].mean()
    bg_std = df.conc.iloc[is_bg.values].std()
    
    # first iteration of plume detection
    conc_bg_removed = (df['conc'] - bg).rolling(bg_subtraction_window).mean()
    is_plume = conc_bg_removed >= plume_sd_threshold * bg_std
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
            .sort_values(['start', 'end'])
    )
    plume_intervals = [pd.Interval(s, e, closed="both") for s, e in zip(plume_groups['start'], plume_groups['end'])]
    
    # Combine overlapping plumes
    buffer_time = datetime.timedelta(seconds=plume_buffer)
    def reduce_intervals(acc, el):
        plume_with_buffer = pd.Interval(acc[0].left, acc[0].right + buffer_time)
        if plume_with_buffer.overlaps(el):
            return [pd.Interval(acc[0].left, el.right)] + acc[1:]
        else:
            return [el] + acc
    
    plume_overlap = reduce(reduce_intervals, plume_intervals[1::], [plume_intervals[0]])
    
    # Just turn into DataFrame for user friendliness
    plumes_condensed = pd.DataFrame([{'start': x.left, 'end': x.right} for x in plume_overlap]).sort_values(['start'])
    
    fig, ax = plt.subplots()
    myFmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    # TODO get label from args
    ax.set_ylabel('Concentration')
    ax.set_xlabel('Time / UTC')
    ax.plot(df.conc, color='gray', alpha=.5)
    
    areas = []
    for row in plumes_condensed.itertuples():
        raw = df.conc.loc[row.start:row.end]
        bg_removed = raw - bg.loc[row.start:row.end]
        ax.plot(raw)
        this_df = pd.DataFrame([{
            'start': row.start,
            'end': row.end,
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

