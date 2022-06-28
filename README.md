# ACRUISE Peak ID

Python and R functions for identifying peaks from measurements taken from ACRUISE campaigns.

## Python package

### Installation

To install the Python package, run `pip install git+https://github.com/wacl-york/acruise-peakid#subdirectory=acruisepy`.

### Usage

Below shows an example of loading some CO2 data into Pandas for use with `acruisepy` and subsetting it to a period of interest.
The first column in the CSV is a datetime, which has been set to to the DataFrame's index.

```python
import datetime
import pandas as pd
from acruisepy import peakid

df_co2 = pd.read_csv(
    "C258_FGGA_FAAM.txt",
    names=["conc", "a", "b"],
    parse_dates=True,
    header=1,
    date_parser=lambda x: datetime.datetime.strptime(x, "%d/%m/%Y %H:%M:%S.%f"),
)
start_time = datetime.datetime(2021, 10, 4, 9, 40, 0)
end_time = datetime.datetime(2021, 10, 4, 13, 30, 0)
df_co2 = df_co2.loc[start_time : end_time]
```

The first step is to identify the background level, against which plumes will be compared.
The `peakid.identify_background` function takes in the concentration and several tuning parameters, here the default values are used but it is highly likely you will need to tune them for your dataset. 
Run `help(peakid.identify_background)` to see a description of these parameters.

```python
bg = peakid.identify_background(df_co2['conc'], bg_sd_window=180, bg_sd_threshold=0.5, bg_mean_window=660)
```

To determine the appropriateness of the extracted background, `plot_background` plots the concentration time-series with the background highlighted alongside the limit of what is considered a plume, defined as `plume_sd_threshold` standard deviations about the mean background.
Any plumes will then be determined to start from the point at which they cross `plume_sd_starting` standard deviations above the mean background.
If they were considered to start from the point at which they cross the `plume_sd_threshold` boundary instead then some plume data would be lost.

```python
peakid.plot_background(df_co2['conc'], bg, plume_sd_threshold=4, plume_sd_starting=2)
```

![the extracted background concentration level highlighted on the concentration time-series](images/background_py.png)

Once both a satisfactory background and suitable values for `plume_sd_threshold` and `plume_sd_starting` have been identified, the plumes can be detected using `peakid.detect_plumes` (again use `help(peakid.detect_plumes)` to see full details for what arguments it takes).
It's crucial here that both the concentration and the background have a DatetimeIndex, which should be the case if they are loaded in from CSV as shown in the example above.

The plumes can be visually inspected using the `peakid.plot_plumes` function, adjusting the parameters in the previous step until necessary.

```python
plumes = peakid.detect_plumes(df_co2['conc'], bg, plume_sd_threshold=4, plume_sd_starting=2, plume_buffer=5)
peakid.plot_plumes(df_co2['conc'], plumes)
```

![the extracted plumes highlighted on the concentration time-series](images/plumes_py.png)

Once the plumes have been finalised, the area under the plumes can be calculated.
This is currently done using a trapezoidal approach, see the documentation for `numpy.trapz` for further details.
It's important to subtract the background from the concentration time-series for this function.

```python
co2_areas = peakid.integrate_aup_trapz(df_co2['conc'] - bg, plumes, dx=0.1)
co2_areas
```

```
                    start                     end        area
0 2021-10-04 09:52:56.700 2021-10-04 09:52:58.200   11.551411
0 2021-10-04 09:57:34.900 2021-10-04 09:57:35.200    2.463899
0 2021-10-04 10:01:56.200 2021-10-04 10:01:59.500   18.605635
0 2021-10-04 10:05:13.700 2021-10-04 10:05:15.200    4.229698
0 2021-10-04 10:11:53.900 2021-10-04 10:11:57.200    7.469866
0 2021-10-04 10:20:06.100 2021-10-04 10:20:06.700    6.102518
0 2021-10-04 10:31:23.300 2021-10-04 10:31:23.500    0.650395
0 2021-10-04 10:31:31.800 2021-10-04 10:31:50.900   47.789589
...
```

## R package

### Installation

To install the R package, run `remotes::install_github("https://github.com/wacl-york/acruise-peakid", subdir="acruiseR")`

### Usage

Below shows an example of loading some CO2 data for use with `acruiseR` and subsetting it to a period of interest.
The `UTC_time` column contains the timestamps which are recorded with a sample period of 100ms.
The default R data structure for timestamps (`POSIXct`, which is also used under the hood by the `lubridate` package) cannot accurately store times at millisecond resolution, so instead it is parsed with the `nanotime` package which is accurate to nanoseconds (the same as `pandas`).

```r
library(tidyverse)
library(nanotime)
library(acruiseR)

df_co2 <- read_csv(
    "../data/C258_FGGA_FAAM.txt"
)
df_co2$time_nano <- as.nanotime(df_co2$UTC_time, format="%d/%m/%Y %H:%M:%E3S", tz="UTC")
# Restrict to time period of interest if necessary
study_window <- as.nanoival("+2021-10-04T09:40:00UTC -> 2021-10-04T13:30:00UTC+")
df_co2 <- df_co2[df_co2$time_nano %in% study_window,]
```

The first step is to identify the background level, against which plumes will be compared.
The `identify_background` function takes in the concentration and several tuning parameters, here the default values are used but it is highly likely you will need to tune them for your dataset. 
Run `?identify_background` to see a description of these parameters.

```r
bg <- identify_background(df_co2$CO2_ppm, bg_sd_window=180, bg_sd_threshold=0.5, bg_mean_window=660)
```

To determine the appropriateness of the extracted background, `plot_background` plots the concentration time-series with the background highlighted alongside the limit of what is considered a plume, defined as `plume_sd_threshold` standard deviations about the mean background.
Any plumes will then be determined to start from the point at which they cross `plume_sd_starting` standard deviations above the mean background.
If they were considered to start from the point at which they cross the `plume_sd_threshold` boundary instead then some plume data would be lost.

```r
plot_background(df_co2$CO2_ppm, df_co2$time_nano, bg, plume_sd_threshold=4, plume_sd_starting=2)
```

![the extracted background concentration level highlighted on the concentration time-series](images/background_r.png)

This function returns a `ggplot2` object, so any additional tweaks can be made.

```r
p <- plot_background(df_co2$CO2_ppm, df_co2$time_nano, bg, plume_sd_threshold=4, plume_sd_starting=2)
p + theme(legend.position="bottom")
```

![the same plot as before but with the legend underneath](images/background_r_2.png)

Once both a satisfactory background and suitable values for `plume_sd_threshold` and `plume_sd_starting` have been identified, the plumes can be detected using `detect_plumes` (again use `?peakid.detect_plumes` to see full details for what arguments it takes).

The plumes can be visually inspected using the `plot_plumes` function, adjusting the parameters in the previous step until necessary.
Again, the returned plot can be modified through usual `ggplot2` functions.

```r
plumes <- detect_plumes(df_co2$CO2_ppm, bg, df_co2$time_nano, plume_sd_threshold=4, plume_sd_starting=2, plume_buffer=5)
plot_plumes(df_co2$CO2_ppm, df_co2$time_nano, plumes)
```

![the extracted plumes highlighted on the concentration time-series](images/plumes_r.png)

Once the plumes have been finalised, the area under the plumes can be calculated using a trapezoidal approach.
It's crucial to subtract the background from the concentration time-series for this function.

```r
co2_areas <- integrate_aup_trapz(df_co2$CO2_ppm - bg, df_co2$time_nano, plumes, dx=0.1)
co2_areas
```

```
                    start                     end        area
0 2021-10-04 09:52:56.700 2021-10-04 09:52:58.200   11.551411
0 2021-10-04 09:57:34.900 2021-10-04 09:57:35.200    2.463899
0 2021-10-04 10:01:56.200 2021-10-04 10:01:59.500   18.605635
0 2021-10-04 10:05:13.700 2021-10-04 10:05:15.200    4.229698
0 2021-10-04 10:11:53.900 2021-10-04 10:11:57.200    7.469866
0 2021-10-04 10:20:06.100 2021-10-04 10:20:06.700    6.102518
0 2021-10-04 10:31:23.300 2021-10-04 10:31:23.500    0.650395
0 2021-10-04 10:31:31.800 2021-10-04 10:31:50.900   47.789589
...
```

