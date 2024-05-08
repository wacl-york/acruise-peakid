# ACRUISE Peak ID

R and Python functions for identifying peaks from measurements taken from ACRUISE campaigns.

NB: There are 3 distinct peak detection methods used in this package:

  1. Using wavelets to directly extract the plumes (only in Python) **Recommended**
  2. Using a Generalized Additive Model (GAM) to estimate the background (only in the R package)
  3. Using a rolling averaging approach to estimate the background (in both R and Python packages)

The wavelet method is the preferred choice owing to its speed, limited number of parameters, and direct approach in estimating peaks compared to the 2-step method used by the other two algorithms.
If the Python package isn't an option, the GAM method in the R package works well and only has a single parameter to tune, but it can be rather slow.
The third method isn't recommended anymore owing to its large number of parameters.

Refer to the READMEs within the 2 sub-directories `acruisepy` and `acruiseR` for specific details of the Python and R packages respectively.

## Methodology

The workflow of the peak extraction comprises 2 or 3 steps depending on the algorithm.
The wavelet method simply directly extracts the peaks themselves, while the other two approaches first estimate the background and then extract the peaks via subtraction.
Once the peaks are established the plume concentrations can simply be identified by integrating the area under the curve via trapezoidal approximation.

There are no global optimum parameters, and instead an iterative approach is necessary to identify the best parameter values for a given dataset.
  
The following sections provide a few additional details.

### Background identification

The observed time-series $y(t)$ is assumed to comprise a background concentration $b(t)$ with irregular, sparse plumes $p(t)$ with background noise having constant variance $\epsilon \sim N(0, \sigma)$:

$$y(t) = b(t) + p(t) + \epsilon(t)$$

The background is estimated using either a Generalized Additive Model (R only for now), or a rolling window approach (both R and Python), giving $\hat{b}(t)$.

### Plume identification

Subtracting $\hat{b}(t)$ from $y(t)$ results in the signal comprising the normally distributed noise and the irregular plumes, from which the noise variance is estimated as the standard deviation ($\hat{\sigma}$)

Plumes are then identified as any timepoints where $y(t) \geq \hat{b}(t) + k \hat{\sigma}$, where $k$ is a user-chosen parameter (`plume_sd_threshold`).
This threshold will usually be at least 3, so if plumes were exclusively defined by this criteria then they would exclude some points between $p(t)$ and $k$ standard deviations.
Instead, for each plume identified by the above criteria, it is then expanded to include points up to $j$ standard deviations away from $p(t)$, where $k > j$. 
The default value for $j$ is 1 and is specified in the code as `plume_sd_starting`.

Plumes that are close together can be combined by using the `plume_buffer` parameter, which combines plumes if they are within `plume_buffer` timepoints.

### Integrating peaks

Currently the only method for integrating the area under the peaks is to use the trapezoidal approximation.
