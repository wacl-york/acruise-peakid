from functools import reduce
import datetime
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.dates as mdates


def identify_background(
    conc: pd.Series,
    bg_sd_window: int = 180,
    bg_sd_threshold: float = 0.5,
    bg_mean_window: int = 660,
) -> pd.Series:
    """
    Identifies background from a concentration time-series.

    This process takes 3 steps:

        - Obtaining a smoothed rolling standard deviation of the background
        - Identifying background measurements as those that lie within a set
          sd threshold
        - Interpolating values outside of this background threshold in a
          linear manner and then using a rolling mean so that 'background'
          measurements are available for the entire time-series.

    Args:
        - conc (pd.Series): Concentration time-series.
        - bg_sd_window (int): Window size for the rolling standard deviation
          smooth to identify the background.
        - bg_sd_threshold (float): Background measurements are considered as
          those whose rolling sd is within this threshold.
        - bg_mean_window (int): The rolling mean to smooth the interpolated
          background

    Returns:
        The smoothed background covering the full times series as a pd.Series
        object.
    """
    # Smooth concentration to get smoothed rolling SD
    roll_std = (
        conc.fillna(method="ffill")
        .rolling(bg_sd_window, center=True)
        .std()
        .rolling(bg_sd_window, center=True)
        .mean()
    )
    # Define background as when rolling std deviation is below a threshold
    is_bg = (~conc.isna()) & (roll_std <= bg_sd_threshold)
    # Interpolate background values when don't have background
    # If want to interpolate everything, so including the NAs introduced by the
    # rolling functions then use limit_direction="both"
    bg = (
        conc.loc[is_bg]
        .reindex(conc.index)
        .interpolate(limit_area="inside")
        .rolling(bg_mean_window)
        .mean()
    )

    return bg


def detect_plumes(
    conc: pd.Series,
    bg: pd.Series,
    plume_sd_threshold: float = 3,
    plume_sd_starting: float = 2,
    plume_buffer: float = 10,
) -> pd.DataFrame:
    """
    Detects plumes in a concentration time series.

    Args:
        - conc (pd.Series): The concentration time-series. Should have a
            Datetime index.
        - bg (pd.Series): The smoothed background time-series, as can be
            obtained from identify_background(). Should have the same
            Datetime index as conc.
        - plume_sd_threshold (float): Plumes are identified as samples
            greater than certain number of standard deviations from a
            smoothed background.
        - plume_sd_starting (float): For any identified plumes, it is determined
            to have started once it reached this many standard deviations above
            the background.
        - plume_buffer (float): A buffer in seconds applied to plumes, so
            that if they are overlapping they are merged into the same plume.

    Returns:
        A pd.DataFrame where each row corresponds to a unique plume, whose
        time boundaries are contained in the 2 columns: `start` and `end`.
    """

    df = pd.DataFrame({"conc": conc, "bg": bg})

    # Derive useful values for identifying plumes
    df["is_plume"] = df["conc"] > (df["bg"] + plume_sd_threshold * df["bg"].std())
    df["is_plume_starting"] = df["conc"] > (
        df["bg"] + plume_sd_starting * df["bg"].std()
    )
    df["plume_group_starting"] = (
        df["is_plume_starting"] != df["is_plume_starting"].shift()
    ).cumsum()

    # Find all groups where concentration are > than the lower starting threshold, and that
    # also have at least one value greater than the higher threshold required to say it is a plume
    plume_groups = (
        df.loc[df["is_plume_starting"]]
        .reset_index()
        .groupby("plume_group_starting")
        .agg(
            has_plume=pd.NamedAgg(column="is_plume", aggfunc="sum"),
            start=pd.NamedAgg(column="index", aggfunc="min"),
            end=pd.NamedAgg(column="index", aggfunc="max"),
        )
        .query("has_plume > 0")
        .drop("has_plume", axis=1)
        .sort_values(["start", "end"])
    )

    # Combine overlapping plumes
    buffer_time = datetime.timedelta(seconds=plume_buffer)
    # Use Interval data structure here which is effectively syntatical sugar
    # around a tuple
    plume_intervals = [
        pd.Interval(s, e, closed="both")
        for s, e in zip(plume_groups["start"], plume_groups["end"])
    ]

    def reduce_intervals(acc, el):
        plume_with_buffer = pd.Interval(acc[0].left, acc[0].right + buffer_time)
        if plume_with_buffer.overlaps(el):
            return [pd.Interval(acc[0].left, el.right)] + acc[1:]
        else:
            return [el] + acc

    plume_overlap = reduce(reduce_intervals, plume_intervals[1::], [plume_intervals[0]])

    # Convert into DataFrame for user friendliness
    plumes_condensed = pd.DataFrame(
        [{"start": x.left, "end": x.right} for x in plume_overlap]
    ).sort_values(["start"])
    return plumes_condensed


def integrate_aup_trapz(
    conc: pd.Series, plumes: pd.DataFrame, dx: float = 1.0
) -> pd.DataFrame:
    """
    Integrate the Area Under a Plume (aup) using a trapezoidal method.

    Args:
        - conc (pd.Series): The concentration time-series with the background
          removed. Must have a Datetime index.
        - plumes (pd.DataFrame): A DataFrame with 'start' and 'end' columns
          containing plume boundaries, as returned by detect_plumes()
        - dx (float): Sampling time, passed onto the dz argument of
          np.trapz.

    Returns:
        A pd.DataFrame with one row per plume and 3 columns `start`, `end`, and
        `area`. The first 2 are the same as in the input `plumes`, while `area`
        contains the integrated area.
    """
    areas = []
    for row in plumes.itertuples():
        this_df = pd.DataFrame(
            [
                {
                    "start": row.start,
                    "end": row.end,
                    "area": np.trapz(conc.loc[row.start : row.end], dx=dx),
                }
            ]
        )
        areas.append(this_df)
    return pd.concat(areas)


def plot_background(
    conc: pd.Series,
    background: pd.DataFrame,
    plume_sd_threshold: float = 3,
    plume_sd_starting: float = 2,
    ylabel: str = "Concentration",
    xlabel: str = "Time (UTC)",
    date_fmt: str = "%H:%M",
    bg_alpha: float = 0.5,
) -> None:
    """
    Plots the concentration time-series highlighting the extracted background
    (red), the threshold for what is considered a plume (orange), and when
    plumes will be determined to have started at (blue).

    Args:
        - conc (pd.Series): The concentration time-series. Must have a Datetime
            index.
        - background (pd.Series): The background time-series, as obtained
          from identify_background(). Must have a Datetime index.
        - plume_sd_threshold (float): Plumes are identified as samples
            greater than certain number of standard deviations from a
            smoothed background.
        - plume_sd_starting (float): For any identified plumes, it is determined
            to have started once it reached this many standard deviations above
            the background.
        - ylabel (str): y-axis label
        - xlabel (str): x-axis label
        - date_fmt (str): How to display the x-axis datetime breaks
        - bg_alpha (float): Alpha level of the background concentration.

    Returns:
        None, plots a figure as a side-effect.
    """
    fig, ax = plt.subplots()
    myFmt = mdates.DateFormatter(date_fmt)
    ax.xaxis.set_major_formatter(myFmt)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.plot(conc, color="gray", alpha=bg_alpha)
    ax.plot(background, color="red", label="Mean background")
    ax.plot(
        background + plume_sd_threshold * background.std(),
        color="orange",
        label="Plume threshold",
    )
    ax.plot(
        background + plume_sd_starting * background.std(),
        color="steelblue",
        label="Plume starting point",
    )
    ax.legend()
    plt.show()


def plot_plumes(
    conc: pd.Series,
    plumes: pd.DataFrame,
    ylabel: str = "Concentration",
    xlabel: str = "Time (UTC)",
    date_fmt: str = "%H:%M",
    bg_alpha: float = 0.5,
) -> None:
    """
    Plots plumes against the background concentration.

    Args:
        - conc (pd.Series): The concentration time-series with the background
          removed. Must have a Datetime index.
        - plumes (pd.DataFrame): A DataFrame with 'start' and 'end' columns
          containing plume boundaries, as returned by detect_plumes()
        - ylabel (str): y-axis label
        - xlabel (str): x-axis label
        - date_fmt (str): How to display the x-axis datetime breaks
        - bg_alpha (float): Alpha level of the background concentration.

    Returns:
        None, plots a figure as a side-effect.
    """
    fig, ax = plt.subplots()
    myFmt = mdates.DateFormatter(date_fmt)
    ax.xaxis.set_major_formatter(myFmt)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.plot(conc, color="gray", alpha=bg_alpha)
    for row in plumes.itertuples():
        ax.plot(conc.loc[row.start : row.end])
    plt.show()
