from functools import reduce
import datetime
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.dates as mdates
import pywt
from typing import Optional
from warnings import warn


def identify_background(
    concentration: pd.Series,
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
        - concentration (pd.Series): Concentration time-series.
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
    warn('The rolling window method is deprecated in favour of the wavelet detection. Please refer to the README.', DeprecationWarning, stacklevel=2)
    # Smooth concentration to get smoothed rolling SD
    roll_std = (
        concentration.fillna(method="ffill")
        .rolling(bg_sd_window, center=True)
        .std()
        .rolling(bg_sd_window, center=True)
        .mean()
    )
    # Define background as when rolling std deviation is below a threshold
    is_bg = (~concentration.isna()) & (roll_std <= bg_sd_threshold)
    # Interpolate background values when don't have background
    # If want to interpolate everything, so including the NAs introduced by the
    # rolling functions then use limit_direction="both"
    bg = (
        concentration.loc[is_bg]
        .reindex(concentration.index)
        .interpolate(limit_area="inside")
        .rolling(bg_mean_window)
        .mean()
    )

    return bg


def detect_plumes(
    concentration: pd.Series,
    background: pd.Series,
    plume_sd_threshold: float = 3,
    plume_sd_starting: float = 2,
    plume_buffer: float = 10,
) -> pd.DataFrame:
    """
    Detects plumes in a concentration time series.

    Args:
        - concentration (pd.Series): The concentration time-series. Must have a
            Datetime index.
        - background (pd.Series): The smoothed background time-series, as can be
            obtained from identify_background(). Must have the same
            Datetime index as concentration.
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
    warn('The rolling window method is deprecated in favour of the wavelet detection. Please refer to the README.', DeprecationWarning, stacklevel=2)
    df = pd.DataFrame({"concentration": concentration, "background": background})
    # Rename index so can reliably refer to it later
    df.index.rename("index", inplace=True)

    # Derive useful values for identifying plumes
    df["is_plume"] = df["concentration"] > (
        df["background"] + plume_sd_threshold * df["background"].std()
    )
    df["is_plume_starting"] = df["concentration"] > (
        df["background"] + plume_sd_starting * df["background"].std()
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

    output = []
    plumes_condensed['plume_id'] = range(1, plumes_condensed.shape[0]+1)
    for row in plumes_condensed.itertuples():
        sub_df = df.loc[row.start : row.end].copy(deep=True)
        sub_df['plume_id'] = row.plume_id
        output.append(sub_df)
    output = pd.concat(output)
    return output[['plume_id', 'concentration']]


def integrate_aup_trapz(
    concentration: pd.Series,
    plumes: pd.DataFrame,
    background: Optional[pd.Series] = None,
    dx: Optional[float] = 1.0
) -> pd.DataFrame:
    """
    Integrate the Area Under a Plume (aup) using a trapezoidal method.

    Args:
        - concentration (pd.Series): The concentration time-series. 
          Must have a Datetime index.
        - plumes (pd.DataFrame): A DataFrame as returned by `detect_plumes()`.
        - background (pd.Series): A time-series of background measurements, the
          same length as `concentration`. Used to subtract the background before
          peak integration. If not provided interpolates linearly over the plume
          duration. Must share an index with `concentration`.
        - dx (float): Sampling time, passed onto the dx argument of
          np.trapezoid.

    Returns:
        A pd.DataFrame with one row per plume and 4 columns:
            - `plume_id`: Integer plume ID, as obtained in `detect_plumes_wavelets`
            - `start`: Start time of this plume
            - `end`: End time of this plume
            - `area`: Total area under the plume.
    """
    if background is None:
        background = concentration.copy(deep=True)
        for plume_id, sub_df in plumes.groupby('plume_id'):
            background.loc[sub_df.index.min():sub_df.index.max()] = np.nan
        background = background.interpolate(method='linear')
    concentration = concentration - background

    areas = []
    for plume_id, sub_df in plumes.groupby('plume_id'):
        start = sub_df.index.min()
        end = sub_df.index.max()
        this_df = pd.DataFrame(
            [
                {
                    "plume_id": plume_id,
                    "start": start,
                    "end": end,
                    "area": np.trapezoid(concentration.loc[start:end], dx=dx),
                }
            ]
        )
        areas.append(this_df)
    return pd.concat(areas)


def plot_background(
    concentration: pd.Series,
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
        - concentration (pd.Series): The concentration time-series. Must have a Datetime
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
    warn('The rolling window method is deprecated in favour of the wavelet detection. Please refer to the README.', DeprecationWarning, stacklevel=2)
    fig, ax = plt.subplots()
    myFmt = mdates.DateFormatter(date_fmt)
    ax.xaxis.set_major_formatter(myFmt)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.plot(concentration, color="gray", alpha=bg_alpha)
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
    concentration: pd.Series,
    plumes: pd.DataFrame,
    ylabel: str = "Concentration",
    xlabel: str = "Time (UTC)",
    date_fmt: str = "%H:%M",
    bg_alpha: float = 0.5,
) -> None:
    """
    Plots plumes against the background concentration.

    Args:
        - concentration (pd.Series): The concentration time-series with the background
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
    ax.plot(concentration, color="gray", alpha=bg_alpha)
    for plume_id, sub_df in plumes.groupby('plume_id'):
        ax.plot(sub_df['concentration'])
    plt.show()

def max_wavelet_level(concentration: pd.Series):
    """
    Identifies the maximum number of levels a signal can be decomposed into
    using wavelets.

    This is a thin wrapper around `pywt.dwt_max_level`.

    Args:
        concentration (pd.Series): The input concentration time-series.

    Returns:
        An integer with the maximum number of levels that can be used in
        `detect_plumes_wavelets`.
    """
    return pywt.dwt_max_level(len(concentration), 'haar')

def detect_plumes_wavelets(concentration: pd.Series,
                           levels: Optional[list[int]] = None,
                           plume_threshold: float =1,
                           plume_starting: float =0.5,
                           plume_buffer: float = 10,
                           interpolate: bool = False,
                           plot: bool = False):
    """
    Detects peaks using a wavelet decomposition.

    This function works differently to the main `detect_plumes`. Rather than
    taking a 2-step approach to identify the background and then remove it leaving
    the plumes behind, the plumes are instead directly identified.

    The concentration time-series is decomposed using a multi-level Haar
    wavelet.
    By reconstructing the signal using a partial set of lower-order Wavelet
    coefficients, the resulting signal has the high-frequency plume content but
    without the slowly-varying background.
    This is a very quick process, so it is recommended to iteratively try
    different levels and plotting the results (with `plot=True`) before
    saving the output plume locations.

    For a gentle introduction to the Haar wavelet, see the first 11 minutes of
    this video:
    https://www.youtube.com/watch?v=c1XL5BeI9_s

    Args:
        - concentration (pd.Series): Concentration time-series, must have a 
            Datetime index.
        - levels (list[int]): Which levels of the decomposition to use when
            reconstructing the peaks
        - plume_threshold (float): The threshold determine whether a signal is a
            plume.
        - plume_starting (float): The threshold from where plumes are determined
            to start.
        - plume_buffer (float): A buffer in seconds applied to plumes, so
            that if they are overlapping they are merged into the same plume.
        - interpolate (bool): Whether to linearly interpolate missing values
            from the raw concentration time-series. Can help in certain conditions.
        - plot (bool): Whether to plot a diagnostic plot of the wavelet
            recomposition with the threshold lines. Used to help identify optimal
            parameter settings.

    Returns:
        A pd.DataFrame where each row corresponds to a timepoint within a plume
        with 3 columns:
          - plume_id: Integer ID identifying the plume this measurement came from
          - concentration: Raw concentration at this timepoint
          - reconstruction: Concentration with the background removed
    """
    # Interpolate missing values - can mess with Wavelets
    if interpolate:
        concentration = concentration.interpolate()
    coefs = pywt.wavedec(concentration, 'haar')
    max_level = max_wavelet_level(concentration)

    # TODO Just pass number levels and let it work out levels for us
    # Wrap levels in list if provided single int
    if type(levels) is int:
        levels = [levels]
    if any(max_level < l < 1 for l in levels):
        raise ValueError("Levels must be between 1 and a maximum level, determined by `max_wavelet_level`")

    # Set unselected levels to zero
    # NB since only ever want the detail and not the approximation (index 0).
    # This has the handy bonus effect of meaning we don't need to explicitly
    # convert between 1-index (User-interface) and 0-index (how they are stored)
    if levels is not None:
        coefs_selected = [np.zeros_like(x) if i not in levels else x for i, x in enumerate(coefs)]

    recon = abs(pywt.waverec(coefs_selected, 'haar'))

    if plot:
        fig, ax = plt.subplots()
        ax.plot((concentration - concentration.mean()).reset_index(drop=True), label=f"Normalised raw signal", alpha=0.5)
        ax.plot(recon, label=f"Reconstructed signal using levels {','.join((str(x) for x in levels))}")
        ax.hlines(plume_threshold, xmin=0, xmax=len(recon), colors='C2',
                   label="plume_threshold")
        ax.hlines(plume_starting, xmin=0, xmax=len(recon), colors='C3',
                   label="plume_starting")
        ax.legend()
        plt.show()

    # Ensure the 2 signals are the same length - not guaranteed!
    con_start = np.max(np.array([0, concentration.size-recon.size]))
    recon_end = recon.size if recon.size <= concentration.size else concentration.size - recon.size
    df = pd.DataFrame({"concentration": concentration[con_start:], "reconstruction": recon[:recon_end]})
    # Rename index so can reliably refer to it later
    df.index.rename("index", inplace=True)

    # Derive useful values for identifying plumes
    df["is_plume"] = df["reconstruction"] > plume_threshold
    df["is_plume_starting"] = df["reconstruction"] > plume_starting
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

    # Join back into time-series and extract concentration
    output = []
    plumes_condensed['plume_id'] = range(1, plumes_condensed.shape[0]+1)
    for row in plumes_condensed.itertuples():
        sub_df = df.loc[row.start : row.end].copy(deep=True)
        sub_df['plume_id'] = row.plume_id
        output.append(sub_df)
    output = pd.concat(output)

    return output[['plume_id', 'concentration']]
