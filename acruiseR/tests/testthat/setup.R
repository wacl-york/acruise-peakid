##### Generates simulated emissions
set.seed(42)
# Parameters
n_obs <- 10000
n_plumes <- 11
mean_bg <- 500
sd_bg <- 1
mean_plume <- 50
sd_plume <- 30
mean_plume_length <- 20
sd_plume_length <- 2

# Simulating the overall time and background is straightforward
time_sim <- seq.POSIXt(from = as.POSIXct("2021-02-18 06:00:00"), by = "1 sec", length.out = n_obs)
bg_sim <- rnorm(n_obs, mean = mean_bg, sd = sd_bg)

# Generate the required number of plumes
plumes <- lapply(1:n_plumes, function(i) {
    duration <- round(abs(rnorm(1, mean = mean_plume_length, sd = sd_plume_length)))
    times <- seq.POSIXt(from = sample(time_sim, 1), length.out = duration, by = "1 sec")
    # Linearly interpolate between background and plume height
    midpoint <- floor(duration / 2)
    concs <- rep(NA, duration)
    concs[1] <- 0
    concs[duration] <- 0
    concs[midpoint] <- rnorm(1, mean_plume, sd_plume)
    interpolated <- approx(concs, xout = which(is.na(concs)))
    concs[interpolated$x] <- interpolated$y
    data.frame(time = times, conc = concs)
})

# Combine
df_sim <- data.frame(
    time = time_sim,
    conc = bg_sim
)
for (i in seq(length(plumes))) {
    plume_times <- df_sim$time %in% plumes[[i]]$time
    df_sim$conc[plume_times] <- df_sim$conc[plume_times] + plumes[[i]]$conc
}
