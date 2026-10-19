knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

# require(devtools) # Install from CRAN if not already installed
# devtools::install_github("ilyamaclean/microclimlearn")
# devtools::install_github("ilyamaclean/micropoint")
# library(microclimlearn)

# # Calculate absorbed radiation
# alb <- 0.23 # albedo
# em <- 0.97  # emissivity
# Rabs_sw <- (1 - alb) * climdata$swdown
# Rabs_lw <- em * climdata$lwdown
# Rabs <- Rabs_sw + Rabs_lw
# 
# # Plot through time
# tme <- as.POSIXct(climdata$obs_time, tz = "UTC")
# plot(Rabs ~ tme, type = "l",
#      xlab = "Month",
#      ylab = expression("Absorbed radiation (W " * m^-2 * ")"))

# # Calculate resistance to heat loss
# zref <- 2        # height above ground of weather data (m)
# veghgt <- 0.12   # vegetation height (m)
# ka <- 0.41       # von Kármán constant
# d <- 0.65 * veghgt  # zero-plane displacement height (m)
# zm <- 0.1 * veghgt  # roughness length for momentum (m)
# zh <- 0.2 * zm      # roughness length for heat (m)
# 
# u2 <- climdata$windspeed
# u2[u2 < 0.5] <- 0.5 # minimum wind speed (m s-1)
# 
# rHa <- (log((zref - d) / zm) * log((zref - d) / zh)) /
#        (ka^2 * u2)  # resistance to heat loss

# # Calculate bulk surface stomatal resistance (rS)
# gsmax <- 0.33 # maximum stomatal conductance of leaves (mol/m^2/s)
# rho <- 43 # approximate molar density of air (mol / m^3)
# # About half of shortwave radiation is PAR, and a joule of PAR is 4.6 umol of photons
# rPAR <- climdata$swdown * 0.5 * 4.6  # Flux density of PAR (umol/m^2/s)
# # The surface conducts about three times as much as its leaves (Kelliher et al. 1995)
# rS <- rho * (rPAR + 300) /  (gsmax * 3 * rPAR) # s/m

# sb <- 5.67e-8 # Stefan–Boltzmann constant W/m^2/K^4
# Rem <- em * sb * (climdata$temp + 273.15) ^ 4 # Emitted radiation (approx) W/m^2
# Rnet <- Rabs - Rem # Net radiation (W/m^2)
# G <- ifelse(Rnet > 0, 0.1 * Rnet, 0.5 * Rnet) # Rate of ground heat storage (W / m^2)

# # Find surface temperature using uniroot method
# Ts1 <- SolveEnergyBalance(climdata$temp, climdata$pres, climdata$relhum,
#                           Rabs, rHa, rS, G, em = 0.97,
#                           method = "uniroot", iters = 1)
# par(mfrow = c(2, 1))
# plot(Ts1 ~ tme, type = "l", col = "firebrick",
#      xlab = "Month", ylab = "Temperature", ylim = c(-10, 65))
# par(new = TRUE)
# plot(climdata$temp ~ tme, type = "l", lwd = 2, col = "gray",
#      ylim = c(-10, 65), xlab = "", ylab = "")
# 
# # Compare with Penman-Monteith solution
# Ts2 <- SolveEnergyBalance(climdata$temp, climdata$pres, climdata$relhum,
#                           Rabs, rHa, rS, G, em = 0.97,
#                           method = "Penman", iters = 1)
# plot(Ts1 ~ Ts2, pch = 15, cex = 0.5,
#      xlim = c(-10, 80), ylim = c(-10, 80),
#      xlab = "Penman-Monteith", ylab = "uniroot")
# abline(a = 0, b = 1, col = "firebrick")
# 
# # Improve Penman-Monteith solution by iterating
# Ts2 <- SolveEnergyBalance(climdata$temp, climdata$pres, climdata$relhum,
#                           Rabs, rHa, rS, G, em = 0.97,
#                           method = "Penman", iters = 4)
# par(new = TRUE)
# plot(Ts1 ~ Ts2, pch = 15, cex = 0.5, col = "steelblue",
#      xlim = c(-10, 80), ylim = c(-10, 80),
#      xlab = "", ylab = "")

# # Calculate wind friction velocity
# uf <- (ka * u2) / log((zref - d) / zm)
# 
# # Select hottest hour
# hr <- which.max(Ts1)
# Ts <- Ts1[hr]              # surface temperature
# Ta <- climdata$temp[hr]    # air temperature
# 
# # Calculate sensible heat flux
# pk <- climdata$pres[hr]
# cp <- cpair(Ta)         # specific heat of air
# ph <- rhoair(Ta, pk)   # molar density of air
# 
# H <- ((ph * cp) / rHa[hr]) * (Ts - Ta)
# 
# # Calculate temperature profile above canopy
# z <- seq(0.12, 2, length.out = 100)
# Tz <- Ts - (H / (ka * ph * cp * uf[hr])) * log((z - d) / zh)
# 
# par(mfrow = c(1, 1))
# plot(z ~ Tz, type = "l", lwd = 2,
#      xlab = "Temperature", ylab = "Height")

# tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
# solp <- sunposition(tme, lat = 49.96807, long = -5.215668)
# 
# plot(solp$zen ~ as.POSIXct(tme), type = "l",
#      xlab = "Month", ylab = "Zenith angle")

# # Select daylight hours on June 21
# sel <- which(tme$mon + 1 == 6 & tme$mday == 21 & solp$zen < 90)
# solp2 <- lapply(solp, `[`, sel)
# 
# # Compute solar coefficient for flat surface
# si <- solarindex(slope = 0, aspect = 180, solp2)
# 
# # Compute extinction coefficient
# k1 <- cank(solp2, x = 1, si)
# k2 <- cank(solp2, x = 0.1, si)
# k3 <- cank(solp2, x = 10, si)
# 
# # Plot extinction coefficient
# plot(k1 ~ solp2$zen, pch = 15, col = "black",
#      xlim = c(0, 90), ylim = c(0, 12),
#      xlab = "Zenith angle", ylab = "Extinction coefficient")
# par(new = TRUE)
# plot(k2 ~ solp2$zen, pch = 15, col = "steelblue",
#      xlim = c(0, 90), ylim = c(0, 12),
#      xlab = "", ylab = "")
# par(new = TRUE)
# plot(k3 ~ solp2$zen, pch = 15, col = "firebrick",
#      xlim = c(0, 90), ylim = c(0, 12),
#      xlab = "", ylab = "")
# 
# legend("topleft", legend = c("x = 1", "x = 0.1", "x = 10"),
#        col = c("black", "steelblue", "firebrick"), pch = 15, bty = "n")

# pai <- 3      # plant area index (m^2 m^-2)
# lref <- 0.4   # leaf reflectance
# ltra <- 0.2   # leaf transmittance
# x <- 1        # leaf inclination coefficient
# gref <- 0.15  # ground reflectance
# slope <- 30   # ground slope (degrees)
# aspect <- 180 # ground aspect (degrees)
# 
# # Compute albedo variables
# albvars <- albedo(pai, lref, ltra, x, gref, slope, aspect, solp,
#                   scattermethod = "Sellers", leafh = FALSE)

# # Direct radiation on a surface facing the sun, from that on flat ground
# Rb0 <- (climdata$swdown - climdata$difrad) / cos(solp$zen * pi / 180)
# 
# # Direct radiation arriving at the sloping surface
# si <- solarindex(slope, aspect, solp)
# Rdir <- Rb0 * si
# 
# # Fraction of the radiation arriving at the slope that is diffuse
# Rtot <- climdata$difrad + Rdir
# difprop <- climdata$difrad / Rtot
# difprop[Rtot == 0] <- 1
# 
# # Calculate albedo
# alb <- difprop * albvars$albd + (1 - difprop) * albvars$albb
# 
# plot(alb ~ as.POSIXct(tme), type = "l", ylim = c(0, 1),
#      xlab = "Month", ylab = "Albedo")

# # Create list of vegetation parameters
# vegvars <- list(h = 0.12,    # vegetation height (m)
#                 pai = 1,     # plant area index
#                 x = 1,       # leaf inclination coefficient
#                 lref = 0.4,  # vegetation reflectance
#                 ltra = 0.2,  # vegetation transmittance
#                 em = 0.97)   # vegetation emissivity
# 
# # Create list of ground parameters
# groundvars <- list(gref = 0.15, # ground reflectance
#                    slope = 30,  # slope (degrees)
#                    aspect = 180,# aspect (degrees)
#                    em = 0.95)   # ground emissivity
# 
# # Calculate absorbed radiation
# Rabs <- calcRabs(climdata, vegvars, groundvars,
#                  lat = 49.96807, long = -5.215668)
# 
# # Plot absorbed radiation through time
# plot(Rabs ~ as.POSIXct(tme), type = "l",
#      xlab = "Month",
#      ylab = expression("Absorbed radiation (W " * m^-2 * ")"))

# # Plot wind profiles for three canopy heights
# uref <- 3            # wind speed at reference height (m s-1)
# zref <- 2            # reference height (m)
# h <- c(0.05, 0.1, 0.5) # vegetation heights (m)
# pai <- 1             # plant area index
# 
# uf <- windfric(uref, zref, h, pai)  # friction velocity
# 
# # Build height grids and wind profiles
# z_list <- lapply(h, function(hi) seq(hi, zref, length.out = 100))
# uz_list <- Map(function(z, ufi, hi)
#                  windprofile_above(z, ufi, hi, pai),
#                z_list, as.list(uf), as.list(h))
# 
# # Bind to matrices for plotting
# Z <- do.call(cbind, z_list)
# UZ <- do.call(cbind, uz_list)
# 
# # Plot wind profiles
# cols <- c("firebrick", "black", "steelblue")
# 
# matplot(UZ, Z, type = "l", lwd = 2, col = cols, lty = 1,
#         xlab = "Wind speed (m s-1)",
#         ylab = "Height (m)",
#         xlim = c(0, 3), ylim = c(0, 2))
# 
# abline(h = h, lty = 2, col = cols)
# 
# legend("topleft",
#        legend = c("h = 0.05 m", "h = 0.1 m", "h = 0.5 m"),
#        col = cols, lwd = 2, bty = "n")

# # Plot temperature profiles for three canopy heights
# tref <- 15        # air temperature at reference height (°C)
# pk <- 101.3       # pressure (kPa)
# H <- 200          # sensible heat flux (W m^-2)
# 
# phcp <- rhoair(tref, pk) * cpair(tref) # volumetric heat capacity of air
# 
# # Resistance to reference height for each canopy
# rHa_ref <- canopyresistance(zref, uf, h, pai)
# 
# # Heat-exchange surface temperature
# th <- tref + (rHa_ref / phcp) * H
# 
# # Resistance and temperature profiles
# rHa_z_list <- Map(function(z, ufi, hi)
#                     canopyresistance(z, ufi, hi, pai),
#                   z_list, as.list(uf), as.list(h))
# 
# tz_list <- Map(function(rHa_i, th_i)
#                  th_i - (rHa_i / phcp) * H,
#                rHa_z_list, as.list(th))
# 
# # Bind for plotting
# TZ <- do.call(cbind, tz_list)
# Z <- do.call(cbind, z_list)
# 
# # Plot profiles
# matplot(TZ, Z, type = "l", lwd = 2, col = cols, lty = 1,
#         xlab = expression("Temperature (" * degree * "C)"),
#         ylab = "Height (m)",
#         xlim = c(15, 25), ylim = c(0, 2))
# 
# abline(h = h, lty = 2, col = cols)
# 
# legend("topright",
#        legend = paste("h =", h, "m"),
#        col = cols, lwd = 2, lty = 1, bty = "n")

# # -------------------------------------------------------------- #
# # ************  Heated surface example ************************* #
# # -------------------------------------------------------------- #
# 
# Rabs <- 1000 # absorbed radiation (W m^-2); high value heats canopy
# z <- seq(0.12, 2, length.out = 100) # heights above a 0.12 m canopy (m)
# 
# # * * Compute wind profile assuming neutral atmospheric conditions * * #
# uf <- windfric(uref = 3, zref = 2, h = 0.12, pai = 1) # friction velocity
# uz1 <- windprofile_above(z, uf, h = 0.12, pai = 1)
# 
# # * * Compute temperature profile assuming neutral conditions * * #
# rHa <- canopyresistance(z = 2, uf, h = 0.12, pai = 1) # aerodynamic resistance
# th <- SolveEnergyBalance(Ta = 15, pk = 101.3, rh = 75, Rabs, rHa, rS = 100,
#                          G = 0, em = 0.97,
#                          method = "uniroot") # canopy temperature
# phcp <- rhoair(15, 101.3) * cpair(15) # volumetric heat capacity of air
# H <- (phcp / rHa) * (th - 15) # sensible heat flux
# rHa <- canopyresistance(z, uf, h = 0.12, pai = 1) # resistance at each height
# tz1 <- th - (rHa / phcp) * H # air temperature profile
# 
# # * * * * * Iterate to derive the Obukhov length * * * * * #
# LL <- Obukhov(15, 101.3, uf, H) # Obukhov length from the neutral solution
# tol <- 1e-9 # convergence tolerance
# error <- tol * 100
# 
# while (error > tol) {
#   olduf <- uf
#   uf <- windfric(uref = 3, zref = 2, h = 0.12, pai = 1, LL) # update friction velocity
#   rHa <- canopyresistance(z = 2, uf, h = 0.12, pai = 1, LL) # update aerodynamic resistance
# 
#   th <- SolveEnergyBalance(Ta = 15, pk = 101.3, rh = 75, Rabs, rHa,
#                            rS = 100, G = 0, em = 0.97,
#                            method = "uniroot") # updated canopy temperature
# 
#   H <- (phcp / rHa) * (th - 15) # updated sensible heat flux
#   LL <- Obukhov(15, 101.3, uf, H) # updated Obukhov length
# 
#   error <- abs(uf - olduf)
# }
# 
# # * * Recalculate profiles including diabatic corrections * * #
# uz2 <- windprofile_above(z, uf, h = 0.12, pai = 1, LL) # wind profile with stability corrections
# rHa <- canopyresistance(z, uf, h = 0.12, pai = 1, LL) # resistance at each height
# tz2 <- th - (rHa / phcp) * H # temperature profile with stability corrections

# # -------------------------------------------------------------- #
# # ************  Cooled surface example ************************* #
# # -------------------------------------------------------------- #
# 
# Rabs <- 250 # absorbed radiation (W m^-2); low value cools canopy
# z <- seq(0.12, 2, length.out = 100) # heights above a 0.12 m canopy (m)
# 
# # * * Compute temperature profile assuming neutral conditions * * #
# uf <- windfric(uref = 3, zref = 2, h = 0.12, pai = 1) # friction velocity
# rHa <- canopyresistance(z = 2, uf, h = 0.12, pai = 1) # aerodynamic resistance
# th <- SolveEnergyBalance(Ta = 15, pk = 101.3, rh = 75, Rabs, rHa, rS = 100,
#                          G = 0, em = 0.97,
#                          method = "uniroot") # canopy temperature
# 
# phcp <- rhoair(15, 101.3) * cpair(15) # volumetric heat capacity of air
# H <- (phcp / rHa) * (th - 15) # sensible heat flux
# 
# rHa <- canopyresistance(z, uf, h = 0.12, pai = 1) # resistance at each height
# uz3 <- uz1 # the neutral wind profile is the same as for the heated surface
# tz3 <- th - (rHa / phcp) * H # air temperature profile
# 
# # * * * * * Iterate to derive the Obukhov length * * * * * #
# LL <- Obukhov(15, 101.3, uf, H) # Obukhov length from the neutral solution
# tol <- 1e-9 # convergence tolerance
# error <- tol * 100
# 
# while (error > tol) {
#   olduf <- uf
#   uf <- windfric(uref = 3, zref = 2, h = 0.12, pai = 1, LL) # update friction velocity
#   rHa <- canopyresistance(z = 2, uf, h = 0.12, pai = 1, LL) # update aerodynamic resistance
# 
#   th <- SolveEnergyBalance(Ta = 15, pk = 101.3, rh = 75, Rabs, rHa,
#                            rS = 100, G = 0, em = 0.97,
#                            method = "uniroot") # updated canopy temperature
# 
#   H <- (phcp / rHa) * (th - 15) # updated sensible heat flux
#   LL <- Obukhov(15, 101.3, uf, H) # updated Obukhov length
# 
#   error <- abs(uf - olduf)
# }
# 
# # * * Recalculate temperature and wind profile * * #
# uz4 <- windprofile_above(z, uf, h = 0.12, pai = 1, LL) # wind profile with stability corrections
# rHa <- canopyresistance(z, uf, h = 0.12, pai = 1, LL) # resistance at each height
# tz4 <- th - (rHa / phcp) * H # temperature profile with stability corrections

# # -------------------------------------------------------------- #
# # ************  Plot both sets of profiles  ******************** #
# # -------------------------------------------------------------- #
# 
# # Compare wind profiles
# plot(z ~ uz1, type = "l", lwd = 2,
#      xlim = c(0, 3), ylim = c(0, 2),
#      xlab = expression("Wind speed (m s"^-1 * ")"),
#      ylab = "Height (m)")
# 
# par(new = TRUE)
# plot(z ~ uz2, type = "l", lwd = 2, col = "firebrick",
#      xlim = c(0, 3), ylim = c(0, 2),
#      xlab = "", ylab = "")
# 
# legend("topleft", legend = c("Neutral", "Diabatic"),
#        col = c("black", "firebrick"), lwd = 2, bty = "n")
# 
# plot(z ~ uz3, type = "l", lwd = 2,
#      xlim = c(0, 3), ylim = c(0, 2),
#      xlab = expression("Wind speed (m s"^-1 * ")"),
#      ylab = "Height (m)")
# 
# par(new = TRUE)
# plot(z ~ uz4, type = "l", lwd = 2, col = "steelblue",
#      xlim = c(0, 3), ylim = c(0, 2),
#      xlab = "", ylab = "")
# 
# legend("topleft", legend = c("Neutral", "Diabatic"),
#        col = c("black", "steelblue"), lwd = 2, bty = "n")
# 
# # Compare temperature profiles
# plot(z ~ tz1, type = "l", lwd = 2,
#      xlim = c(15, 24), ylim = c(0, 2),
#      xlab = expression("Temperature (" * degree * "C)"),
#      ylab = "Height (m)")
# 
# par(new = TRUE)
# plot(z ~ tz2, type = "l", lwd = 2, col = "firebrick",
#      xlim = c(15, 24), ylim = c(0, 2),
#      xlab = "", ylab = "")
# 
# legend("topright", legend = c("Neutral", "Diabatic"),
#        col = c("black", "firebrick"), lwd = 2, bty = "n")
# 
# plot(z ~ tz3, type = "l", lwd = 2,
#      xlim = c(10, 15), ylim = c(0, 2),
#      xlab = expression("Temperature (" * degree * "C)"),
#      ylab = "Height (m)")
# 
# par(new = TRUE)
# plot(z ~ tz4, type = "l", lwd = 2, col = "steelblue",
#      xlim = c(10, 15), ylim = c(0, 2),
#      xlab = "", ylab = "")
# 
# legend("topleft", legend = c("Neutral", "Diabatic"),
#        col = c("black", "steelblue"), lwd = 2, bty = "n")

# # Canopy temperature of a strongly heated surface at a given wind speed
# canopytemp <- function(uref, convection) {
#   phcp <- rhoair(15, 101.3) * cpair(15) # volumetric heat capacity of air
#   H <- 0     # sensible heat flux, none initially
#   LL <- Inf  # neutral conditions initially
#   uf <- 0
#   error <- 1
# 
#   while (error > 1e-9) {
#     olduf <- uf
#     u <- uref
#     # Add the stirring of the air by convection to the measured wind
#     if (convection) u <- effectivewind(uref, H, tc = 15, pk = 101.3)
# 
#     uf <- windfric(u, zref = 2, h = 0.12, pai = 1, LL) # friction velocity
#     rHa <- canopyresistance(z = 2, uf, h = 0.12, pai = 1, LL) # aerodynamic resistance
# 
#     th <- SolveEnergyBalance(Ta = 15, pk = 101.3, rh = 75, Rabs = 1000, rHa,
#                              rS = 100, G = 0, em = 0.97,
#                              method = "Penman", iters = 5) # canopy temperature
# 
#     H <- (phcp / rHa) * (th - 15) # sensible heat flux
#     LL <- Obukhov(15, 101.3, uf, H) # Obukhov length
# 
#     error <- abs(uf - olduf)
#   }
#   return(th)
# }
# 
# # Canopy temperature across a range of wind speeds
# uref <- seq(0.05, 3, length.out = 60) # wind speed at reference height (m/s)
# th1 <- sapply(uref, canopytemp, convection = FALSE) # measured wind only
# th2 <- sapply(uref, canopytemp, convection = TRUE)  # with convection
# 
# plot(th1 ~ uref, type = "l", lwd = 2, ylim = c(20, 90),
#      xlab = "Wind speed (m/s)", ylab = "Canopy temperature")
# lines(th2 ~ uref, lwd = 2, col = "firebrick")
# 
# legend("topright", legend = c("Measured wind", "Effective wind"),
#        col = c("black", "firebrick"), lwd = 2, bty = "n")

# # Wind profile above a forest canopy with and without a roughness sublayer
# uref <- 3   # wind speed at reference height (m/s)
# zref <- 20  # reference height (m)
# h <- 10     # vegetation height (m)
# pai <- 3    # plant area index
# 
# z <- seq(h, zref, length.out = 100) # heights above canopy (m)
# 
# # Standard profile
# uf1 <- windfric(uref, zref, h, pai) # friction velocity
# uz1 <- windprofile_above(z, uf1, h, pai)
# 
# # Profile with a roughness sublayer
# uf2 <- windfric(uref, zref, h, pai, MOST = FALSE) # friction velocity
# uz2 <- windprofile_above(z, uf2, h, pai, MOST = FALSE)
# 
# # Height of the top of the roughness sublayer
# d <- zeroplanedis(h, pai) # zero-plane displacement height
# zw <- d + 2 * (h - d)
# 
# plot(z ~ uz1, type = "l", lwd = 2,
#      xlim = c(1, 3), ylim = c(h, zref),
#      xlab = "Wind speed (m/s)", ylab = "Height (m)")
# lines(z ~ uz2, lwd = 2, col = "firebrick")
# abline(h = zw, lty = 2)
# 
# legend("topleft", legend = c("Standard", "Roughness sublayer"),
#        col = c("black", "firebrick"), lwd = 2, bty = "n")

# # Create list of vegetation parameters
# vegparams <- createplant_inputs("C3")
# 
# # Plot stomatal response to climate inputs
# z <- vegparams$h / 2 # average canopy height
# 
# par(mfrow = c(2, 3))
# 
# # CO2 concentration
# Ca <- 280:450 # CO2 concentration (ppm)
# 
# gs <- stomatalcond_calc(Ca, Rswabs = 500, tair = 20, tleaf = 25,
#                         rh = 75, pk = 101.3, psi_r = -0.005,
#                         vegparams, z, C3 = TRUE)
# 
# plot(gs ~ Ca, type = "l", lwd = 2,
#      xlab = expression(CO[2] ~ "concentration (ppm)"),
#      ylab = "Stomatal conductance",
#      ylim = c(0, 1))
# 
# # Shortwave radiation
# Rswabs <- 0:1000 # radiation absorbed by the leaf (W m^-2), of which 48% is taken to be PAR
# 
# gs <- stomatalcond_calc(Ca = 400, Rswabs, tair = 20, tleaf = 25,
#                         rh = 75, pk = 101.3, psi_r = -0.005,
#                         vegparams, z, C3 = TRUE)
# 
# plot(gs ~ Rswabs, type = "l", lwd = 2,
#      xlab = expression("Absorbed radiation (W " * m^-2 * ")"),
#      ylab = "Stomatal conductance",
#      ylim = c(0, 1))
# 
# # Air temperature
# tair <- (100:300) / 10 # air temperature (°C)
# 
# gs <- stomatalcond_calc(Ca = 400, Rswabs = 500, tair, tleaf = 25,
#                         rh = 75, pk = 101.3, psi_r = -0.005,
#                         vegparams, z, C3 = TRUE)
# 
# plot(gs ~ tair, type = "l", lwd = 2,
#      xlab = expression("Air temperature (" * degree * "C)"),
#      ylab = "Stomatal conductance",
#      ylim = c(0, 1))
# 
# # Leaf temperature
# tleaf <- (100:300) / 10 # leaf temperature (°C)
# 
# gs <- stomatalcond_calc(Ca = 400, Rswabs = 500, tair = 20, tleaf,
#                         rh = 75, pk = 101.3, psi_r = -0.005,
#                         vegparams, z, C3 = TRUE)
# 
# plot(gs ~ tleaf, type = "l", lwd = 2,
#      xlab = expression("Leaf temperature (" * degree * "C)"),
#      ylab = "Stomatal conductance",
#      ylim = c(0, 1))
# 
# # Relative humidity
# rh <- 50:100 # relative humidity (%)
# 
# gs <- stomatalcond_calc(Ca = 400, Rswabs = 500, tair = 20, tleaf = 25,
#                         rh, pk = 101.3, psi_r = -0.005,
#                         vegparams, z, C3 = TRUE)
# 
# plot(gs ~ rh, type = "l", lwd = 2,
#      xlab = "Relative humidity (%)",
#      ylab = "Stomatal conductance",
#      ylim = c(0, 1))
# 
# # Atmospheric pressure
# pk <- c(980:1040) / 10 # pressure (kPa)
# 
# gs <- stomatalcond_calc(Ca = 400, Rswabs = 500, tair = 20, tleaf = 25,
#                         rh = 75, pk, psi_r = -0.005,
#                         vegparams, z, C3 = TRUE)
# 
# plot(gs ~ pk, type = "l", lwd = 2,
#      xlab = "Atmospheric pressure (kPa)",
#      ylab = "Stomatal conductance",
#      ylim = c(0, 1))

# # Check soil parameter table
# ?micropoint::soilparams
# 
# # Select values for clay loam
# Smin <- 0.091 # residual volumetric water content
# Smax <- 0.419 # volumetric water content at saturation
# b <- 5.2      # shape parameter
# psie <- 2.6   # air-entry water potential (J/kg)
# 
# theta <- seq(Smin, Smax, length.out = 100) # volumetric water content
# psi_r <- PsiFromtheta(theta, psie, b, Smax) # soil water potential (MPa)
# 
# gs <- stomatalcond_calc(Ca = 400, Rswabs = 500, tair = 20, tleaf = 25,
#                         rh = 75, pk = 101.3, psi_r,
#                         vegparams, z, C3 = TRUE)
# 
# par(mfrow = c(1, 1))
# 
# plot(gs ~ theta, type = "l", lwd = 2,
#      xlab = expression("Volumetric water content (m"^3 * " m"^-3 * ")"),
#      ylab = "Stomatal conductance",
#      ylim = c(0, 0.5))

# # ===================================================================== #
# # ~~~~~~~~~~ Compute leaf stomatal conductance for C3 grass ~~~~~~~~~~~ #
# # ===================================================================== #
# 
# # Compute sun position at midday on 21 May 2025 in Cornwall, UK
# tme <- as.POSIXlt(0, origin = "2025-05-21 12:00", tz = "UTC")
# solp <- sunposition(tme, 50, -5)
# 
# # Derive vegetation parameters
# vegparams <- createplant_inputs("C3")
# 
# # Compute absorbed PAR radiation
# Rsw <- 800  # shortwave radiation (W m^-2)
# Rdif <- 300 # diffuse radiation (W m^-2)
# 
# om <- vegparams$lrefp + vegparams$ltrap # leaf reflectance + transmittance (PAR)
# 
# k <- cank(solp, x = vegparams$x,
#           cos(solp$zen * pi / 180)) # canopy extinction coefficient
# 
# # Compute absorbed radiation
# # -- (Rsw - Rdif) gives direct beam radiation
# # -- k determines the fraction intercepted by leaves
# # -- Rdif adds diffuse radiation separately
# # -- (1 - om) accounts for absorbed rather than scattered radiation
# Rswabs <- (1 - om) * (k * (Rsw - Rdif) + Rdif)
# 
# # Compute stomatal conductance
# gs <- stomatalcond_calc(Ca = 430, Rswabs, tair = 11, tleaf = 30,
#                         rh = 75, pk = 101.3, psi_r = 0,
#                         vegparams, z = vegparams$h / 2)
# 
# # ===================================================================== #
# # ~~~~~~~~ Compute bulk stomatal conductance for C3 grass ~~~~~~~~~~~~~ #
# # ===================================================================== #
# 
# PAI <- (1:400) / 100 # plant area index
# Lfrac <- 1 # fraction of plant area that is living leaf
# 
# # Compute bulk stomatal resistance
# rS <- 0
# vegp2 <- vegparams
# 
# for (i in 1:length(PAI)) {
#   vegp2$pai <- PAI[i]
# 
#   rS[i] <- bulkstomatalresist_calc(solp, Ca = 430, Rsw, Rdif,
#                                    tair = 11, tcanopy = 30,
#                                    rh = 75, pk = 101.3,
#                                    psi_r = 0, vegp2, Lfrac)
# }
# 
# # Convert to bulk stomatal conductance
# gS <- rhoair(11, 101.3) / rS
# 
# # ===================================================================== #
# # Compute ratio of bulk surface to leaf conductance and plot results
# # ===================================================================== #
# 
# ratio <- gS / gs
# 
# par(mar = c(6, 6, 3, 3))
# 
# plot(ratio ~ PAI, type = "l", lwd = 2,
#      cex.axis = 1.3, cex.lab = 1.3,
#      xlab = "Plant area index",
#      ylab = expression(g[S] / g[s]))
# 
# abline(h = 1, lty = 2, col = "firebrick")

# # ----------------------- Create data inputs ------------------------------- #
# weather <- climdata
# 
# tme <- as.POSIXlt(weather$obs_time, tz = "UTC")
# sel <- which(tme$mon + 1 == 6) # select hours in June
# 
# weather <- weather[sel,]
# 
# groundp <- micropoint::groundparams
# 
# vegp <- createplant_inputs("C3")
# vegp$h <- 0.12  # vegetation height (m)
# vegp$pai <- 0.5 # plant area index
# Lfrac <- 0.6    # fraction of plant area that is living leaf
# 
# Ca <- micropoint::Cafromyear(2017) # atmospheric CO2 concentration (ppm)
# 
# # ------------------- Estimate rate of ground heat storage ----------------- #
# 
# # Calculate diurnal temperature range
# tday <- matrix(weather$temp, nrow = 24)
# 
# tmax <- apply(tday, 2, max, na.rm = TRUE)
# tmin <- apply(tday, 2, min, na.rm = TRUE)
# 
# dtr <- rep(tmax - tmin, each = 24)
# 
# # Simulate sinusoidal daily cycle peaking at 3 pm
# daycycle <- sin(pi / 12 * ((0:719) - 9))
# 
# # Scale by diurnal temperature range
# G <- sqrt(0.5 * dtr) * 25 * daycycle
# 
# # -------- Simulate soil water content and surface humidity ---------------- #
# 
# theta <- (0.1 * sin(pi / 4380 * ((0:8759) + 2000)) + 0.3)[sel] # root zone
# thetas <- (0.15 * sin(pi / 4380 * ((0:8759) + 2000)) + 0.21)[sel] # surface
# 
# soilrh <- with(groundp,
#                soilrelhum(thetas, 15, Psie, b, Smax))
# 
# reqhgt <- 0.2 # target height for air temperature (m)
# 
# # ------------------------- Cycle through each hour ------------------------ #
# 
# Tz <- 0 # initialise temperature vector
# 
# for (hr in 1:720) {
# 
#   # --------------------- Run canopy solver -------------------------------- #
#   can <- solve_wholecanopy(weather, vegp, groundp, hr,
#                            lat = 49.96807, long = -5.215668,
#                            zref = 2, G[hr],
#                            theta[hr], soilrh[hr], Ca, Lfrac = Lfrac)
# 
#   # ---------------- Calculate temperature at reqhgt ----------------------- #
# 
#   # Volumetric heat capacity of air
#   phcp <- rhoair(weather$temp[hr], weather$pres[hr]) *
#           cpair(weather$temp[hr])
# 
#   # Resistance to heat transfer between the canopy and reqhgt
#   rHa <- canopyresistance(reqhgt, can$uf, vegp$h, vegp$pai, can$LL)
# 
#   # Air temperature at reqhgt
#   Tz[hr] <- can$tcanopy - (can$H / phcp) * rHa
# }
# 
# # ----------------------------- Plot results ------------------------------- #
# 
# tmem <- as.POSIXct(tme[sel])
# 
# # Plot reference air temperature
# plot(weather$temp ~ tmem, type = "l",
#      xlab = "Date",
#      ylab = expression("Temperature (" * degree * "C)"),
#      ylim = c(5, 35),
#      lwd = 2,
#      col = rgb(0, 0, 0, 0.5))
# 
# # Plot microclimate air temperature 20 cm above ground
# par(new = TRUE)
# 
# plot(Tz ~ tmem, type = "l",
#      xlab = "Date",
#      ylab = "",
#      ylim = c(5, 35),
#      col = "firebrick")
# 
# legend("topright",
#        legend = c("Reference air temperature",
#                   "Microclimate temperature"),
#        col = c(rgb(0, 0, 0, 0.5), "firebrick"),
#        lwd = 2,
#        bty = "n")

# # Use inbuilt vegetation parameter dataset from micropoint package
# vegp <- micropoint::forestparams
# 
# # Generate vertical profile of plant area index using micropoint::PAIgeometry
# paii <- micropoint::PAIgeometry(PAI = vegp$pai,
#                                 skew = vegp$skew,
#                                 spread = vegp$spread,
#                                 n = 1000)
# 
# z <- ((1:1000) / 1000) * vegp$h # height above ground (m)
# 
# plot(z ~ paii, type = "l", lwd = 2,
#      xlab = "Plant area index of each layer",
#      ylab = "Height (m)",
#      main = paste0("Total PAI: ", round(sum(paii), 2)))

# # Extract data for a sunny hour
# weather <- climdata
# 
# swdown <- weather$swdown[4094] # incoming shortwave radiation (W m^-2)
# difrad <- weather$difrad[4094] # diffuse radiation (W m^-2)
# 
# # Calculate sun position at that hour
# tme <- as.POSIXlt(weather$obs_time[4094], tz = "UTC")
# solp <- sunposition(tme, lat = 49.96807, long = -5.215668)

# # Run two-stream model and plot radiation fluxes
# groundp <- micropoint::groundparams
# 
# tsmod <- twostream(vegp, groundp, paii,
#                    swdown, difrad, solp)
# 
# plot(z ~ tsmod$Rdirdown, type = "l", col = "firebrick", lwd = 2,
#      xlab = expression("Radiation flux (W " * m^-2 * ")"),
#      ylab = "Height (m)",
#      xlim = c(0, 750), ylim = c(0, 10))
# 
# par(new = TRUE)
# plot(z ~ tsmod$Rdifdown, type = "l", col = "steelblue", lwd = 2,
#      xlab = "", ylab = "",
#      xlim = c(0, 750), ylim = c(0, 10))
# 
# par(new = TRUE)
# plot(z ~ tsmod$Rswup, type = "l", col = "black", lwd = 2,
#      xlab = "", ylab = "",
#      xlim = c(0, 750), ylim = c(0, 10))
# 
# par(new = TRUE)
# plot(z ~ tsmod$Rleafabs, type = "l", col = "forestgreen", lwd = 2,
#      xlab = "", ylab = "",
#      xlim = c(0, 750), ylim = c(0, 10))
# 
# legend("bottomright",
#        legend = c("Direct downward",
#                   "Diffuse downward",
#                   "Upward shortwave",
#                   "Leaf absorption"),
#        col = c("firebrick", "steelblue", "black", "forestgreen"),
#        lwd = 2,
#        bty = "n")

# # Create PAI profile using micropoint::PAIgeometry
# paii <- micropoint::PAIgeometry(PAI = 4,
#                                 skew = 7,
#                                 spread = 70,
#                                 n = 100)
# 
# # Generate leaf temperature profile
# z <- (1:100) / 10 # height above ground (m)
# tleaf <- 3 * cos(1 - z / 10) + 13.5 # leaf temperature (°C)
# 
# # Call longwavebelow function
# lwb <- longwavebelow(paii,
#                      lwdown = 300,
#                      tleaf,
#                      tground = 15,
#                      vegem = 0.97)
# 
# # Plot longwave radiation fluxes
# plot(z ~ lwb$Rlwdown, type = "l", col = "firebrick", lwd = 2,
#      xlab = expression("Radiation flux (W " * m^-2 * ")"),
#      ylab = "Height (m)",
#      xlim = c(250, 450), ylim = c(0, 10))
# 
# par(new = TRUE)
# plot(z ~ lwb$Rlwup, type = "l", col = "steelblue", lwd = 2,
#      xlab = "", ylab = "",
#      xlim = c(250, 450), ylim = c(0, 10))
# 
# par(new = TRUE)
# plot(z ~ lwb$RlwLabs, type = "l", col = "forestgreen", lwd = 2,
#      xlab = "", ylab = "",
#      xlim = c(250, 450), ylim = c(0, 10))
# 
# legend("bottomleft",
#        legend = c("Downward longwave",
#                   "Upward longwave",
#                   "Leaf absorption"),
#        col = c("firebrick", "steelblue", "forestgreen"),
#        lwd = 2,
#        bty = "n")

# paii_skew <- micropoint::PAIgeometry(PAI = 4,
#                                      skew = 3,
#                                      spread = 70,
#                                      n = 100) # skewed PAI profile
# 
# paii_unif <- rep(4 / 100, 100) # uniform PAI profile
# 
# # ----------------- Generate wind profiles ------------------------- #
# 
# uh <- 3 # wind speed at canopy top (m s^-1)
# 
# wind_skew <- windprofile_below(hgt = 25, paii_skew, uh)
# wind_unif <- windprofile_below(hgt = 25, paii_unif, uh)
# 
# # ----------------------- Plot results ----------------------------- #
# 
# z <- ((1:100) - 0.5) / 100 * 25 # height of the middle of each layer (m)
# 
# par(mfrow = c(1, 2))
# 
# # Foliage profiles
# plot(z ~ paii_skew, type = "l", lwd = 2, col = "forestgreen",
#      xlab = "Plant area index of each layer",
#      ylab = "Height (m)",
#      xlim = c(0, 0.06))
# 
# par(new = TRUE)
# plot(z ~ paii_unif, type = "l", lwd = 2, col = "steelblue",
#      xlab = "", ylab = "",
#      xlim = c(0, 0.06))
# 
# legend("topright",
#        legend = c("Skewed profile", "Uniform profile"),
#        col = c("forestgreen", "steelblue"),
#        lwd = 2,
#        bty = "n")
# 
# # Wind profiles
# plot(z ~ wind_skew, type = "l", lwd = 2, col = "forestgreen",
#      xlab = expression("Wind speed (m s"^-1 * ")"),
#      ylab = "Height (m)",
#      xlim = c(0, 3))
# 
# par(new = TRUE)
# plot(z ~ wind_unif, type = "l", lwd = 2, col = "steelblue",
#      xlab = "", ylab = "",
#      xlim = c(0, 3))
# 
# legend("bottomright",
#        legend = c("Skewed profile", "Uniform profile"),
#        col = c("forestgreen", "steelblue"),
#        lwd = 2,
#        bty = "n")

# uref <- 3  # wind speed at reference height (m s^-1)
# zref <- 20 # reference height (m)
# h <- 7     # vegetation height (m)
# z <- seq(0, zref, length.out = 500) # heights above ground (m)
# 
# # Dense canopy
# uf_dense <- windfric(uref, zref, h, pai = 3, MOST = FALSE) # friction velocity
# KH_dense <- thermaldiffusivity(z, uf_dense, h, pai = 3)
# 
# # Sparse canopy
# uf_sparse <- windfric(uref, zref, h, pai = 0.01, MOST = FALSE)
# KH_sparse <- thermaldiffusivity(z, uf_sparse, h, pai = 0.01)
# 
# # Bare ground with a roughness length of 0.004 m
# uf_bare <- 0.41 * uref / log(zref / 0.004)
# KH_bare <- thermaldiffusivity(z, uf_bare, h = 0, pai = 0)
# 
# # Plot the three profiles
# plot(z ~ KH_dense, type = "l", lwd = 2, col = "forestgreen",
#      xlab = expression("Thermal diffusivity (" * m^2 ~ s^-1 * ")"),
#      ylab = "Height (m)",
#      xlim = c(0, 2.5))
# lines(z ~ KH_sparse, lwd = 2, col = "darkorange")
# lines(z ~ KH_bare, lwd = 2, col = "black", lty = 2)
# abline(h = h, lty = 3) # canopy height
# 
# legend("bottomright",
#        legend = c("Dense canopy", "Sparse canopy", "Bare ground"),
#        col = c("forestgreen", "darkorange", "black"),
#        lty = c(1, 1, 2), lwd = 2, bty = "n")

# # -------------------------- Derive model inputs -------------------------- #
# groundp <- micropoint::groundparams
# vegp <- createplant_inputs("BDT") # deciduous broadleaf trees
# paii <- micropoint::PAIgeometry(PAI = vegp$pai,
#                                 skew = 7,
#                                 spread = 70,
#                                 n = 100)
# 
# # -------------------- Set temperatures and humidity ---------------------- #
# tleaf <- rep(26, 100) # initial leaf temperature (°C)
# tair <- rep(25, 100)  # air temperature (°C)
# tground <- 25         # ground temperature (°C)
# rh <- rep(80, 100)    # relative humidity (%)
# 
# # ---------------------- Calculate absorbed radiation --------------------- #
# # Run two-stream model
# tme <- as.POSIXlt(0, origin = "2026-05-21 12:00", tz = "UTC")
# solp <- sunposition(tme,
#                     lat = 49.96807,
#                     long = -5.215668)
# tsmod <- twostream(vegp, groundp, paii,
#                    swdown = 839,
#                    difrad = 223,
#                    solp)
# 
# # Run longwave model
# lwb <- longwavebelow(paii,
#                      lwdown = 300,
#                      tleaf,
#                      tground,
#                      vegem = 0.97)
# Rabs <- lwb$RlwLabs + tsmod$Rleafabs # absorbed radiation (W m^-2)
# 
# # ----------------------- Calculate heat resistance ----------------------- #
# uh <- 2 # wind speed at canopy top (m s^-1)
# uz <- windprofile_below(vegp$h, paii, uh) # wind profile within canopy
# 
# dT <- abs(tleaf - tair) # temperature difference - used with low wind speed
# rHa <- leafresistance(tair, dT, uz,
#                       vegp$len, vegp$wid, vegp$x)
# 
# # -------------------- Calculate stomatal resistance ---------------------- #
# # Radiation absorbed by sunlit and shaded leaves, with leaf properties for PAR
# tsmod <- twostream(vegp, groundp, paii,
#                    swdown = 839,
#                    difrad = 223,
#                    solp,
#                    PAR = TRUE)
# 
# z <- ((1:100) - 0.5) / 100 * vegp$h # height of the middle of each layer (m)
# 
# # Sunlit leaf stomatal conductance (mol m^-2 s^-1)
# gs_sun <- stomatalcond_calc(Ca = 430,
#                             tsmod$Rsun,
#                             tair,
#                             tleaf,
#                             rh,
#                             pk = 101.3,
#                             psi_r = -0.5,
#                             vegp,
#                             z)
# 
# # Shaded leaf stomatal conductance (mol m^-2 s^-1)
# gs_shade <- stomatalcond_calc(Ca = 430,
#                               tsmod$Rshade,
#                               tair,
#                               tleaf,
#                               rh,
#                               pk = 101.3,
#                               psi_r = -0.5,
#                               vegp,
#                               z)
# 
# # Combine sunlit and shaded fractions
# gs <- gs_sun * tsmod$sunlitfrac +
#       (1 - tsmod$sunlitfrac) * gs_shade
# 
# ph <- rhoair(tair, pk = 101.3) # air density (mol m^-3)
# rS <- ph / gs # stomatal resistance (s m^-1)
# 
# # ---------------------- Derive leaf temperature -------------------------- #
# tleaf <- SolveEnergyBalance(tair,
#                             pk = 101.3,
#                             rh,
#                             Rabs,
#                             rHa,
#                             rS,
#                             G = 0,
#                             em = 0.97,
#                             method = "Penman",
#                             iters = 4)
# 
# # ---------------------- Plot leaf temperature profile -------------------- #
# par(mfrow = c(1, 1))
# plot(z ~ tleaf, type = "l", lwd = 2, col = "forestgreen",
#      xlab = expression("Leaf temperature (" * degree * "C)"),
#      ylab = "Height (m)",
#      xlim = c(24, 28))

# # ------------------------- Foliage profile ------------------------------- #
# paii <- wythampai$paii # plant area index of each layer
# z <- wythampai$z       # height of the middle of each layer (m)
# hgt <- 24.5            # canopy height (m)
# pai <- sum(paii)       # plant area index of the whole canopy
# 
# # ---------------- Heat released by each canopy layer --------------------- #
# Rnet <- 500 # radiation at canopy top available to heat the canopy (W m^-2)
# k <- 0.5    # extinction coefficient: how strongly foliage attenuates radiation
# 
# # Fraction of that radiation passing the top and the bottom of each layer
# pabove <- rev(cumsum(rev(paii))) - paii # plant area above the top of each layer
# trtop <- exp(-k * pabove)
# trbottom <- exp(-k * (pabove + paii))
# 
# # Half of the radiation absorbed by each layer is released as sensible heat
# Hz <- 0.5 * Rnet * (trtop - trbottom) # sensible heat flux (W m^-2 of ground)
# 
# # ----------------------- Far-field temperature --------------------------- #
# uf <- 0.2       # friction velocity (m s^-1)
# ttop <- 20      # air temperature at canopy top (°C)
# tground <- 19.5 # ground surface temperature (°C)
# 
# # Heights from the ground surface to canopy top (m), closely spaced near the ground
# zall <- c(0, 0.02, 0.05, 0.1, z, hgt)
# tfar <- farfield(Hz, hgt, pai, uf,
#                  top = ttop,
#                  ground = tground,
#                  tair = ttop, # sets the density and heat capacity of air
#                  z = zall)
# 
# # ---------------------------- Plot results ------------------------------- #
# par(mfrow = c(1, 2))
# plot(z ~ paii, type = "l", lwd = 2, col = "forestgreen",
#      xlab = "Plant area index of layer",
#      ylab = "Height (m)")
# plot(zall ~ tfar, type = "l", lwd = 2, col = "firebrick",
#      xlab = expression("Far-field temperature (" * degree * "C)"),
#      ylab = "Height (m)")

# # ----------------------- Near-field temperature -------------------------- #
# tnear <- nearfield(Hz, hgt, pai, uf,
#                    tair = ttop, # sets the density and heat capacity of air
#                    z = zall)
# 
# # --------------------------- Air temperature ----------------------------- #
# tair <- tfar + tnear
# 
# # ---------------------------- Plot results ------------------------------- #
# par(mfrow = c(1, 2))
# plot(zall ~ tnear, type = "l", lwd = 2, col = "steelblue",
#      xlab = expression("Near-field temperature (" * degree * "C)"),
#      ylab = "Height (m)")
# plot(zall ~ tair, type = "l", lwd = 2, col = "black",
#      xlab = expression("Air temperature (" * degree * "C)"),
#      ylab = "Height (m)")
# lines(zall ~ tfar, lwd = 2, col = "firebrick", lty = 2)
# legend("left",
#        legend = c("Far + near", "Far only"),
#        col = c("black", "firebrick"),
#        lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)

# # -------------------------- Derive model inputs -------------------------- #
# weather <- climdata
# groundp <- micropoint::groundparams
# vegp <- createplant_inputs("BDT") # broadleaf deciduous trees
# paii <- wythampai$paii            # plant area index of each layer
# z <- wythampai$z                  # height of the middle of each layer (m)
# n <- length(paii)                 # number of canopy layers
# vegp$h <- 24.5                    # canopy height (m)
# vegp$pai <- sum(paii)             # plant area index of the whole canopy
# Ca <- micropoint::Cafromyear(2017) # atmospheric CO2 concentration (ppm)
# 
# # --------- Run canopy solver to derive conditions at canopy top ---------- #
# lat <- 49.96807
# long <- -5.215668
# hr <- 4094
# can <- solve_wholecanopy(weather, vegp, groundp, hr,
#                          lat, long,
#                          zref = 30,
#                          G = 20,
#                          theta = 0.2,
#                          soilrh = 87,
#                          Ca,
#                          MOST = FALSE,
#                          Lfrac = 0.6)
# 
# # ---------------------- Set reference climate variables ------------------ #
# lwdown <- weather$lwdown[hr] # downward longwave radiation (W m^-2)
# swdown <- weather$swdown[hr] # downward shortwave radiation (W m^-2)
# difrad <- weather$difrad[hr] # diffuse radiation (W m^-2)
# pk <- weather$pres[hr]       # pressure (kPa)
# uf <- can$uf                 # friction velocity (m s^-1)
# LL <- can$LL                 # Obukhov length (m)
# uh <- can$uh                 # wind speed at canopy top (m s^-1)
# th <- can$th                 # air temperature at canopy top (°C)
# 
# # Relative humidity at canopy top (%), from vapour pressure at reference height
# eref <- satvap(weather$temp[hr]) * weather$relhum[hr] / 100
# rhtop <- 100 * eref / satvap(th)
# 
# tground <- th - 1 # ground surface temperature (°C)
# soilrh <- 87      # effective relative humidity of the soil surface (%)
# psi_r <- -0.5     # root-zone water potential (MPa)
# 
# # ---------------- Initial guesses of below-canopy profiles --------------- #
# tleaf <- rep(th + 1, n) # leaf temperature (°C)
# tair <- rep(th, n)      # air temperature (°C)
# rh <- rep(rhtop, n)     # relative humidity (%)
# 
# # --------------------- Calculate shortwave radiation --------------------- #
# tme <- as.POSIXlt(weather$obs_time[hr], tz = "UTC")
# solp <- sunposition(tme, lat, long) # zenith and azimuth angles
# tsmod <- twostream(vegp, groundp, paii,
#                    swdown, difrad, solp)
# # With leaf properties for PAR, for the stomatal calculations
# tspar <- twostream(vegp, groundp, paii,
#                    swdown, difrad, solp,
#                    PAR = TRUE)
# 
# # ------------------------- Calculate wind profile ------------------------ #
# uz <- windprofile_below(vegp$h, paii, uh) # wind speed of each layer (m s^-1)
# 
# # -------------------------- Initialise while loop ------------------------ #
# error <- 1e99
# itr <- 0
# 
# while (error > 1e-3) {
# 
#   # ------------------- Calculate absorbed radiation ---------------------- #
#   lwb <- longwavebelow(paii,
#                        lwdown,
#                        tleaf,
#                        tground,
#                        vegem = 0.97)
#   Rabs <- lwb$RlwLabs + tsmod$Rleafabs # absorbed radiation (W m^-2)
# 
#   # --------------------- Calculate heat resistance ----------------------- #
#   dT <- abs(tleaf - tair) # leaf-air temperature difference (°C)
#   rHa <- leafresistance(tair, dT, uz,
#                         vegp$len, vegp$wid, vegp$x)
# 
#   # ------------------ Calculate stomatal resistance ---------------------- #
#   gs_sun <- stomatalcond_calc(Ca,
#                               tspar$Rsun,
#                               tair,
#                               tleaf,
#                               rh,
#                               pk,
#                               psi_r,
#                               vegp,
#                               z)
# 
#   gs_shade <- stomatalcond_calc(Ca,
#                                 tspar$Rshade,
#                                 tair,
#                                 tleaf,
#                                 rh,
#                                 pk,
#                                 psi_r,
#                                 vegp,
#                                 z)
# 
#   gs <- gs_sun * tspar$sunlitfrac +
#         (1 - tspar$sunlitfrac) * gs_shade
# 
#   ph <- rhoair(tair, pk) # air density (mol m^-3)
#   rS <- ph / gs          # stomatal resistance (s m^-1)
# 
#   # --------------------- Derive leaf temperature ------------------------- #
#   tleaf <- SolveEnergyBalance(tair,
#                               pk,
#                               rh,
#                               Rabs,
#                               rHa,
#                               rS,
#                               G = 0,
#                               em = 0.97,
#                               method = "Penman",
#                               iters = 4)
# 
#   # ------------- Heat and vapour released by each layer ------------------ #
#   # Per unit leaf surface (W m^-2)
#   cp <- cpair(tair)                # heat capacity of air (J mol^-1 K^-1)
#   es <- satvap(tleaf)              # vapour pressure of the leaf (kPa)
#   ea <- satvap(tair) * (rh / 100)  # vapour pressure of the air (kPa)
#   H <- (ph * cp / rHa) * (tleaf - tair)                    # sensible heat
#   L <- (latvap(tleaf) * ph / ((rHa + rS) * pk)) * (es - ea) # latent heat
# 
#   # Per unit ground area: both sides of the leaves, times plant area
#   Hz <- 2 * paii * H
#   Lz <- 2 * paii * L
# 
#   # --------------------- Calculate air temperature ----------------------- #
#   tairn <- temp_below(Hz, vegp$h, vegp$pai, uf,
#                       ttop = th,
#                       tground = tground,
#                       tair = tair,
#                       pk = pk,
#                       LL = LL)
# 
#   # --------------------- Calculate relative humidity --------------------- #
#   rh <- relhum_below(Lz, vegp$h, vegp$pai, uf,
#                      rhtop = rhtop,
#                      ttop = th,
#                      tground = tground,
#                      soilrh = soilrh,
#                      tair = tairn,
#                      pk = pk,
#                      LL = LL)
# 
#   # ------------------------ Check for convergence ------------------------ #
#   error <- max(abs(tairn - tair))
#   tair <- tairn
#   itr <- itr + 1
# }
# 
# # ----------------------------- Plot results ------------------------------ #
# par(mfrow = c(1, 2))
# # Temperature profiles
# plot(z ~ tleaf, type = "l", col = "forestgreen", lwd = 2,
#      xlim = c(24, 32),
#      xlab = expression("Temperature (" * degree * "C)"),
#      ylab = "Height (m)")
# lines(z ~ tair, col = "steelblue", lwd = 2)
# legend("bottomright",
#        legend = c("Leaf", "Air"),
#        col = c("forestgreen", "steelblue"),
#        lwd = 2,
#        bty = "n")
# 
# # Relative humidity profile of the air
# plot(z ~ rh, type = "l", col = "steelblue", lwd = 2,
#      xlim = c(50, 70),
#      xlab = "Relative humidity (%)",
#      ylab = "Height (m)")

# # Composition of a typical clay loam
# Vq <- 0.0568 # volumetric quartz fraction
# Vm <- 0.4818 # volumetric fraction of other minerals
# Vo <- 0.0114 # volumetric organic fraction
# Mc <- 0.34 # mass fraction of clay
# 
# # Volumetric water fraction, from almost dry to saturated
# Vw <- seq(0.02, 0.45, length.out = 100)
# 
# # Thermal conductivity (W m^-1 K^-1)
# k <- thermalConductivity(Vq, Vm, Vo, Vw, Mc,
#                          Tc = 15,    # soil temperature (°C)
#                          pk = 101.3) # atmospheric pressure (kPa)
# 
# # Volumetric heat capacity (MJ m^-3 K^-1)
# Ch <- heatCapacity(Vq, Vm, Vo, Vw,
#                    Tc = 15,
#                    pk = 101.3) / 1e6
# 
# # Plot both against water content
# par(mfrow = c(1, 2))
# plot(k ~ Vw, type = "l", lwd = 2, col = "steelblue",
#      xlab = expression("Volumetric water fraction (m"^3 * " m"^-3 * ")"),
#      ylab = expression("Thermal conductivity (W m"^-1 * " K"^-1 * ")"))
# plot(Ch ~ Vw, type = "l", lwd = 2, col = "firebrick",
#      xlab = expression("Volumetric water fraction (m"^3 * " m"^-3 * ")"),
#      ylab = expression("Volumetric heat capacity (MJ m"^-3 * " K"^-1 * ")"))

# # ----------------------------- Model inputs ------------------------------ #
# weather <- climdata
# tme <- as.POSIXlt(weather$obs_time, tz = "UTC")
# sel <- which(tme$mon + 1 == 6 &
#              tme$mday == 21) # select hours of 21 June
# weather <- weather[sel, ]
# 
# groundp <- micropoint::groundparams # ground parameters
# 
# vegp <- createplant_inputs("C3") # vegetation parameters
# vegp$h <- 0 # vegetation height (m): zero for bare ground
# vegp$pai <- 0 # plant area index: zero for bare ground
# 
# # --------------------- Thermal properties of the soil -------------------- #
# # Clay loam with a volumetric water fraction of 0.3
# k <- thermalConductivity(Vq = 0.0568, Vm = 0.4818, Vo = 0.0114,
#                          Vw = 0.3, Mc = 0.34,
#                          Tc = 15, pk = 101.3) # thermal conductivity (W m^-1 K^-1)
# Ch <- heatCapacity(Vq = 0.0568, Vm = 0.4818, Vo = 0.0114,
#                    Vw = 0.3,
#                    Tc = 15, pk = 101.3) # heat capacity (J m^-3 K^-1)
# 
# tt <- 0:23        # hours of day
# om <- 2 * pi / 24 # angular frequency of the daily cycle (radians per hour)
# ka <- k / Ch      # thermal diffusivity (m^2 s^-1)
# DampDepth <- sqrt(2 * ka * 3600 / om) # damping depth (m)
# 
# # -------------------------- Initialise variables ------------------------- #
# G <- rep(0, 24) # ground heat flux (W m^-2)
# error <- 1e99   # difference from previous iteration
# Told <- 0       # temperature in previous iteration
# itr <- 1        # iteration counter
# 
# cols <- adjustcolor(colorRampPalette(c("firebrick", "steelblue"))(12),
#                     alpha.f = 0.25)
# 
# while (error > 1e-3) {
#   # ------------------------ Cycle through each hour ---------------------- #
#   Ts <- 0 # initialise surface temperature
# 
#   for (hr in 1:24) {
#     can <- solve_wholecanopy(weather, vegp, groundp,
#                              hr,
#                              lat = 49.96807,
#                              long = -5.215668,
#                              zref = 2,
#                              G[hr],
#                              theta = 0.3,
#                              soilrh = 85)
#     Ts[hr] <- can$tcanopy # surface temperature (°C)
#   }
# 
#   # Plot diurnal temperature on each iteration
#   if (itr == 1) {
#     plot(Ts ~ tt, type = "l", col = "firebrick",
#          ylim = c(10, 40),
#          xlab = "Hour",
#          ylab = expression("Temperature (" * degree * "C)"),
#          lwd = 2)
#   } else {
#     lines(Ts ~ tt, col = cols[min(itr, 12)], lwd = 2)
#   }
# 
#   # ---------- Fit a sinusoidal cycle to surface temperature -------------- #
#   m <- lm(Ts ~ sin(om * tt) + cos(om * tt))
#   # Decimal hour at which the fitted cycle peaks
#   hrmx <- as.numeric(((pi / 2 -
#             atan2(coef(m)[3], coef(m)[2])) / om) %% 24)
#   A0 <- (max(Ts) - min(Ts)) / 2 # amplitude of temperature cycle (°C)
# 
#   # ------------------------ Calculate ground flux ------------------------ #
#   # Peaks three hours (a phase of pi / 4) before surface temperature
#   G <- sqrt(2) * A0 * k / DampDepth *
#        sin(om * (tt - hrmx) + 3 * pi / 4)
# 
#   error <- max(abs(Ts - Told)) # convergence error
#   Told <- Ts
#   itr <- itr + 1
# }
# 
# # ---------------------- Plot final converged temperature ----------------- #
# lines(Ts ~ tt, col = "steelblue", lwd = 2)

# # Assumes code above already run
# Tav <- mean(Ts) # average temperature (°C)
# 
# # ---------------- Temperatures at 5 cm depth ---------------------------- #
# z <- -0.05 # depth below surface (m)
# Tz5 <- Tav + A0 * exp(z / DampDepth) *
#        sin(om * (tt - hrmx) + pi / 2 + z / DampDepth)
# 
# # ---------------- Temperatures at 20 cm depth --------------------------- #
# z <- -0.2 # depth below surface (m)
# Tz20 <- Tav + A0 * exp(z / DampDepth) *
#         sin(om * (tt - hrmx) + pi / 2 + z / DampDepth)
# 
# # Assumed sinusoidal cycle at the surface
# Tz0 <- Tav + A0 * sin(om * (tt - hrmx) + pi / 2)
# 
# # ----------------------------- Plot results ----------------------------- #
# plot(Ts ~ tt, type = "l", col = "firebrick",
#      ylim = c(10, 40),
#      xlab = "Hour",
#      ylab = expression("Temperature (" * degree * "C)"),
#      lwd = 2)
# lines(Tz0 ~ tt, col = "firebrick", lwd = 2, lty = 2)
# lines(Tz5 ~ tt, col = "black", lwd = 2)
# lines(Tz20 ~ tt, col = "steelblue", lwd = 2)
# legend("topleft",
#        legend = c("Surface temperature",
#                   "Sinusoidal approximation",
#                   "5 cm depth",
#                   "20 cm depth"),
#        col = c("firebrick", "firebrick", "black", "steelblue"),
#        lty = c(1, 2, 1, 1),
#        lwd = 2,
#        bty = "n")

# # -------------- Generate list of soil model parameters --------------- #
# soilp <- createsoilplist(soiltype = "Clay loam",
#                          nlayers = 15, # number of soil layers
#                          totalDepth = 2, # total soil depth (m)
#                          slope = 0, # slope (degrees)
#                          aspect = 180, # aspect (degrees)
#                          surface_organicmu = 3, # surface organic matter depth
#                          FreeDrain = TRUE) # free drainage lower boundary

# # Initial soil temperature sequence: surface to lower boundary layer
# Te <- seq(18, 11, length.out = 16)
# theta <- rep(0.3, 16) # volumetric soil water fraction
# soilheat <- InitSoilheatmod(soilp, Te, theta)

# # ---------------- Create driving data from inbuilt datasets -------------- #
# weather <- climdata
# groundp <- micropoint::groundparams
# 
# vegp <- createplant_inputs("C3")
# vegp$h <- 0 # vegetation height (m): zero for bare ground
# vegp$pai <- 0 # total plant area index: zero for bare ground
# 
# lat <- 49.96807 # latitude (decimal degrees)
# long <- -5.215668 # longitude (decimal degrees)
# 
# boundaryT <- mean(weather$temp) # lower boundary temperature (°C)
# 
# # ---------------- Select data for first 10 days in June ----------------- #
# tme <- as.POSIXlt(weather$obs_time, tz = "UTC")
# 
# sel <- which(tme$mon + 1 == 6 &
#              tme$mday < 11)
# 
# weather <- weather[sel,]
# tme <- tme[sel]
# 
# # ---------------- Create parameters for soil heat model ----------------- #
# soilp <- createsoilplist(soiltype = "Clay loam",
#                          nlayers = 15, # number of soil layers
#                          totalDepth = 2, # total soil depth (m)
#                          slope = 0, # slope (degrees)
#                          aspect = 180, # aspect (degrees)
#                          surface_organicmu = 3, # surface organic matter depth
#                          FreeDrain = TRUE) # free drainage lower boundary
# 
# Te <- seq(weather$temp[1],
#           boundaryT,
#           length.out = 16) # initial temperature profile (°C)
# 
# theta <- rep(0.3, 16) # volumetric water content
# 
# soilheat <- InitSoilheatmod(soilp, Te, theta) # initialise soil heat model
# 
# # -------------------------- Initialise variables ------------------------- #
# error <- 1e99 # error tolerance for iterative solve
# G <- 0 # ground heat flux (W m^-2)
# 
# Tsurface <- rep(0, 240) # ground surface temperature (°C)
# Gflux <- rep(0, 240) # ground heat flux for each timestep
# 
# for (hr in 1:240) {
# 
#   while (error > 1e-4) {
# 
#      # ------ Absorbed radiation and resistance to heat transfer --------- #
#      can <- solve_wholecanopy(weather, vegp, groundp,
#                               hr,
#                               lat,
#                               long,
#                               zref = 2,
#                               G,
#                               theta = 0.3,
#                               soilrh = 100) # moist soil: air at its surface is saturated
# 
#      # ------------------------- Run soil heat model --------------------- #
#      soilheat <- SoilHeatModel(soilheat,
#                                soilp,
#                                Rabs = can$Rabs, # radiation absorbed by the ground
#                                Tref = weather$temp[hr],
#                                relhum = weather$relhum[hr],
#                                atmPressure = weather$pres[hr],
#                                rHa = can$rHg) # resistance from ground to zref
# 
#      error <- abs(G - soilheat$Gflux) # convergence error
#      G <- soilheat$Gflux
#   }
# 
#   Gflux[hr] <- G
#   Tsurface[hr] <- soilheat$Te[1]
# 
#   soilheat$oldTe <- soilheat$Te # update model state
# 
#   error <- 1e99 # reset convergence criterion
# }
# 
# # ----------------------------- Plot time series -------------------------- #
# par(mfrow = c(2, 2))
# 
# tmec <- as.POSIXct(tme, tz = "UTC")
# 
# # Soil surface temperature
# plot(Tsurface ~ tmec, type = "l", lwd = 2, col = "firebrick",
#      xlab = "Date",
#      ylab = expression("Temperature (" * degree * "C)"))
# 
# # Ground heat flux
# plot(Gflux ~ tmec, type = "l", lwd = 2, col = "steelblue",
#      xlab = "Date",
#      ylab = expression("Ground heat flux (W " * m^-2 * ")"))
# 
# # ---------------- Plot normalised diurnal cycle ------------------------- #
# 
# # Values are scaled between each day's minimum and maximum
# # to emphasise differences in diurnal shape among days
# 
# plotdailycycle <- function(v, clr) {
# 
#   # Convert hourly vector into day × hour matrix
#   vd <- matrix(v, ncol = 24, byrow = TRUE)
# 
#   # Daily minima and maxima
#   vmn <- apply(vd, 1, min)
#   vmx <- apply(vd, 1, max)
# 
#   # Scale values between 0 and 1
#   vdfrac <- (vd - vmn) / (vmx - vmn)
# 
#   # Plot each daily cycle
#   matplot(t(vdfrac), type = "l", lty = 1,
#           col = clr,
#           xlab = "Hour",
#           ylab = "Normalised diurnal cycle")
# }
# 
# # Plot daily cycles for surface temperature and ground heat flux
# plotdailycycle(Tsurface, rgb(0.7, 0.13, 0.13, 0.5))
# plotdailycycle(Gflux, rgb(0.27, 0.51, 0.71, 0.5))

# # ------------ Comparison of Campbell and Van Genuchten --------------------- #
# Campbell_params <- Campbellparams("Clay loam")
# VG_params <- VGparams("Clay loam")
# 
# theta <- seq(0.1, 0.35, length.out = 100) # volumetric water content
# 
# psiw1 <- with(Campbell_params,
#               PsiFromtheta(theta,
#                            psie,
#                            b,
#                            Smax,
#                            MPa = FALSE))
# 
# psiw2 <- with(VG_params,
#               PsiFromthetaVG(theta,
#                              Smin,
#                              Smax,
#                              alpha,
#                              n))
# 
# # ----------------------------- Plot results ------------------------------ #
# par(mar = c(5, 5, 2, 2),
#     mfrow = c(1, 1))
# plot(theta, log10(-psiw1),
#      type = "l",
#      lwd = 2,
#      ylim = c(log10(1e4), log10(1)),
#      yaxt = "n",
#      xlab = expression(theta ~ "(m"^3 * " m"^-3 * ")"),
#      ylab = expression(psi[w] ~ "(J kg"^-1 * ")"),
#      cex.lab = 1.3,
#      cex.axis = 1.3)
# lines(theta, log10(-psiw2),
#       col = "firebrick",
#       lwd = 2)
# exps <- 0:4
# axis(2,
#      at = log10(10^exps),
#      cex.axis = 1.3,
#      labels = paste0("-", 10^exps))
# legend("bottomright",
#        legend = c("Campbell", "van Genuchten"),
#        col = c("black", "firebrick"),
#        lwd = 2,
#        bty = "n",
#        cex = 1.1)

# VG_params1 <- VGparams("Clay")
# VG_params2 <- VGparams("Sand")
# 
# psiw1 <- with(VG_params1,
#               PsiFromthetaVG(theta,
#                              Smin,
#                              Smax,
#                              alpha,
#                              n))
# 
# psiw2 <- with(VG_params2,
#               PsiFromthetaVG(theta,
#                              Smin,
#                              Smax,
#                              alpha,
#                              n))
# 
# # ----------------------------- Plot results ------------------------------ #
# par(mar = c(5, 5, 2, 2))
# plot(theta, log10(-psiw1),
#      type = "l",
#      lwd = 2,
#      ylim = c(log10(1e10), log10(0.1)),
#      yaxt = "n",
#      xlab = expression(theta ~ "(m"^3 * " m"^-3 * ")"),
#      ylab = expression(psi[w] ~ "(J kg"^-1 * ")"),
#      cex.lab = 1.3,
#      cex.axis = 1.3)
# lines(theta, log10(-psiw2),
#       col = "firebrick",
#       lwd = 2)
# exps <- -1:10
# axis(2,
#      at = log10(10^exps),
#      cex.axis = 1.3,
#      labels = paste0("-", 10^exps))
# legend("bottomright",
#        legend = c("Clay", "Sand"),
#        col = c("black", "firebrick"),
#        lwd = 2,
#        bty = "n",
#        cex = 1.1)

# Campbell_params <- Campbellparams("Clay loam")
# VG_params <- VGparams("Clay loam")
# 
# Ks1 <- with(Campbell_params,
#             hydraulicConductivityFromTheta(theta,
#                                            b, # Campbell shape parameter
#                                            Smax, # volumetric water content at saturation
#                                            Ksat)) # saturated hydraulic conductivity
# 
# Ks2 <- with(VG_params,
#             hydraulicConductivityFromThetaVG(theta,
#                                              Smin, # residual water content
#                                              Smax, # volumetric water content at saturation
#                                              n, # van Genuchten shape parameter
#                                              Ko, # conductivity at the matching point
#                                              L)) # pore connectivity parameter
# 
# # ----------------------------- Plot results ------------------------------ #
# par(mar = c(5, 5, 2, 2))
# plot(theta, Ks1,
#      log = "y",
#      type = "l",
#      lwd = 2,
#      ylim = c(1e-14, 1e-5),
#      xlab = expression(theta ~ "(m"^3 * " m"^-3 * ")"),
#      ylab = expression("Conductivity (kg s m"^-3 * ")"),
#      cex.axis = 1.3,
#      cex.lab = 1.3)
# lines(theta, Ks2,
#       col = "firebrick",
#       lwd = 2)
# legend("bottomright",
#        legend = c("Campbell", "van Genuchten"),
#        col = c("black", "firebrick"),
#        lwd = 2,
#        bty = "n",
#        cex = 1.1)

# # ---------------- Load example forcing and parameter datasets ------------ #
# weather <- climdata # meteorological forcing data
# groundp <- micropoint::groundparams # ground surface parameters
# 
# vegp <- createplant_inputs("C3")
# vegp$h <- 0 # vegetation height (m): zero for bare ground
# vegp$pai <- 0 # plant area index: zero for bare ground
# 
# lat <- 49.96807 # latitude (decimal degrees)
# long <- -5.215668 # longitude (decimal degrees)
# 
# boundaryT <- mean(weather$temp) # lower soil boundary temperature (°C)
# 
# # ---------------- Select first 10 days of June -------------------------- #
# tme <- as.POSIXlt(weather$obs_time, tz = "UTC")
# sel <- which(tme$mon + 1 == 6 &
#              tme$mday < 11)
# weather <- weather[sel,]
# tme <- tme[sel]
# 
# # ---------------------- Create soil parameter set ----------------------- #
# soilp <- createsoilplist(soiltype = "Clay loam",
#                          nlayers = 15, # number of soil layers
#                          totalDepth = 2, # total soil depth (m)
#                          slope = 0, # slope (degrees)
#                          aspect = 180, # aspect (degrees)
#                          surface_organicmu = 1, # surface organic matter depth
#                          FreeDrain = TRUE) # free drainage lower boundary
# 
# n <- length(soilp$b) # number of depths: one more than the number of layers
# 
# # ---------------------- Initialise soil profiles ------------------------ #
# Te <- seq(weather$temp[1],
#           boundaryT,
#           length.out = n) # initial temperature profile (°C)
# theta <- seq(0.28,
#              0.36,
#              length.out = n) # initial soil water profile
# soilheat <- InitSoilheatmod(soilp,
#                             Te,
#                             theta) # initialise soil heat model
# soilwater <- InitSoilwatermod(soilp,
#                               soilheat,
#                               root50 = vegp$root50, # rooting depths (m): not
#                               root95 = vegp$root95) # used without transpiration
# 
# # Identify depth closest to 10 cm
# nd <- which.min(abs(soilheat$z - 0.1))
# 
# # ------------------------- Initialise outputs --------------------------- #
# Tsurface <- T10cm <- theta_surface <- theta_10cm <- rep(NA, 240)
# G <- 0 # initial ground heat flux estimate
# error <- 1e99 # convergence criterion
# 
# # --------------------------- Main simulation loop ----------------------- #
# for (hr in 1:240) {
#   itr <- 1
#   # Iteratively solve soil heat and soil water balance
#   while (error > 1e-2 && itr < 20) {
#     # ------------- Relative humidity of air at the soil surface -------- #
#     soilrh <- soilrelhum(soilwater$theta[1], # surface water content
#                          soilheat$Te[1], # surface temperature (°C)
#                          soilp$psie[1], # air-entry water potential
#                          soilp$b[1], # Campbell shape parameter
#                          soilp$Smax[1]) # saturation water content
# 
#     # ------- Absorbed radiation and resistance to heat transfer -------- #
#     can <- solve_wholecanopy(weather,
#                              vegp,
#                              groundp,
#                              hr,
#                              lat,
#                              long,
#                              zref = 2,
#                              G,
#                              theta = soilwater$theta[1],
#                              soilrh = soilrh)
# 
#     # ----------------------- Solve soil heat transport ----------------- #
#     soilheat <- SoilHeatModel(soilheat,
#                               soilp,
#                               Rabs = can$Rabs, # absorbed radiation (W m^-2)
#                               Tref = weather$temp[hr], # air temperature (°C)
#                               relhum = weather$relhum[hr], # relative humidity (%)
#                               atmPressure = weather$pres[hr], # pressure (kPa)
#                               rHa = can$rHg) # resistance from ground to zref (s m^-1)
# 
#     error <- abs(G - soilheat$Gflux) # convergence error
#     G <- soilheat$Gflux
# 
#     # ----------------------- Assemble forcing data --------------------- #
#     forcing <- list(Tair = weather$temp[hr], # air temperature (°C)
#                     relhum = weather$relhum[hr], # relative humidity (%)
#                     rHa = can$rHg, # resistance from ground to zref (s m^-1)
#                     precip = weather$precip[hr], # precipitation (mm h^-1)
#                     Et = 0) # transpiration (mm h^-1): none from bare ground
# 
#     # ----------------------- Solve soil water movement ----------------- #
#     soilwater$Tc <- soilheat$Te
#     soilwater$oldTc <- soilheat$oldTe
# 
#     soilwater <- SoilWaterModel(soilwater,
#                                 soilp,
#                                 forcing)
# 
#     # Update thermal properties using new soil moisture
#     soilheat$wc <- soilwater$theta
#     itr <- itr + 1
#   }
# 
#   # ---------------------------- Store outputs -------------------------- #
#   Tsurface[hr] <- soilheat$Te[1]
#   T10cm[hr] <- soilheat$Te[nd]
# 
#   theta_surface[hr] <- soilwater$theta[1]
#   theta_10cm[hr] <- soilwater$theta[nd]
# 
#   # State at the end of the hour becomes that at the start of the next
#   soilheat$oldTe <- soilheat$Te
#   soilwater$oldvapor <- soilwater$vapor
#   soilwater$oldtheta <- soilwater$theta
# 
#   error <- 1e99 # reset convergence criterion
# }
# 
# # --------------------------- Plot time series --------------------------- #
# par(mfrow = c(2, 1))
# tmec <- as.POSIXct(tme, tz = "UTC")
# # Soil temperature time series
# plot(Tsurface ~ tmec, type = "l",
#      lwd = 2,
#      col = "firebrick",
#      ylim = c(5, 30),
#      xlab = "Date",
#      ylab = expression("Temperature (" * degree * "C)"))
# lines(T10cm ~ tmec, lwd = 2, col = "steelblue")
# legend("topright",
#        legend = c("Surface", "10 cm"),
#        col = c("firebrick", "steelblue"),
#        lwd = 2,
#        bty = "n")
# # Soil moisture time series
# plot(theta_surface ~ tmec, type = "l",
#      lwd = 2,
#      col = "firebrick",
#      ylim = c(0, 0.5),
#      xlab = "Date",
#      ylab = expression(theta ~ "(m"^3 * " m"^-3 * ")"))
# lines(theta_10cm ~ tmec, lwd = 2, col = "steelblue")
# legend("topright",
#        legend = c("Surface", "10 cm"),
#        col = c("firebrick", "steelblue"),
#        lwd = 2,
#        bty = "n")

# snowage <- seq(0, 60 * 24, length.out = 500) # age of the snow (hours)
# alb <- SnowAlbedo(snowage)
# par(mfrow = c(1, 1))
# plot(alb ~ I(snowage / 24), type = "l",
#      col = "steelblue",
#      lwd = 2,
#      ylim = c(0, 1),
#      xlab = "Snow age (days)",
#      ylab = "Albedo")

# # ---------------- Load example forcing data ------------------------------ #
# weather <- climdata # meteorological forcing data
# 
# # ------------------------ Select data for January ----------------------- #
# tme <- as.POSIXlt(weather$obs_time, tz = "UTC")
# sel <- which(tme$mon + 1 == 1)
# weather <- weather[sel,]
# tme <- tme[sel]
# 
# # Reduce air temperature by 9 °C to give snow, and reduce longwave
# # radiation from the sky to that of air colder by the same amount
# temp_new <- weather$temp - 9
# weather$lwdown <- weather$lwdown *
#   ((temp_new + 273.15) / (weather$temp + 273.15))^4
# weather$temp <- temp_new
# 
# # ---------------------- Create soil parameter set ----------------------- #
# soilp <- createsoilplist(soiltype = "Clay loam",
#                          nlayers = 15, # number of soil layers
#                          totalDepth = 2, # total soil depth (m)
#                          slope = 0, # slope (degrees)
#                          aspect = 180, # aspect (degrees)
#                          surface_organicmu = 1, # surface organic matter depth
#                          FreeDrain = TRUE) # free drainage lower boundary
# n <- length(soilp$b) # number of depths: one more than the number of layers
# 
# # --------------------------- Initialise models -------------------------- #
# Te <- seq(0, 4, length.out = n) # initial soil temperature profile (°C)
# theta <- rep(0.3, n) # volumetric soil water content
# soilheat <- InitSoilheatmod(soilp, Te, theta) # initialise soil heat model
# snowmod <- InitSnowModel(snowdepth = 0, # no snow to begin with
#                          snowage = 0,
#                          snowtemp = 0)
# 
# # ------------------------- Initialise outputs --------------------------- #
# nhrs <- nrow(weather)
# snowdepth <- surfacetemp <- soiltemp <- snowfall <- melt <- drainage <- rep(0, nhrs)
# LL <- Inf # Obukhov length (m): neutral air to begin with
# 
# # --------------------- Run model in hourly timesteps -------------------- #
# for (hr in 1:nhrs) {
#   depth <- sum(snowmod$depth) # snow depth (m)
# 
#   # ---------------- Properties of the surface: snow or soil ------------- #
#   if (depth > 0) {
#     alb <- SnowAlbedo(snowmod$age[1]) # albedo, from the age of the surface snow
#     em <- 0.99 # emissivity of snow
#     z0 <- 0.002 # roughness length of snow (m)
#   } else {
#     alb <- soilp$gref # albedo of the soil
#     em <- soilp$groundem # emissivity of the soil
#     z0 <- 0.004 # roughness length of bare ground (m)
#   }
# 
#   # ------------------------ Absorbed radiation ------------------------- #
#   Rabs <- (1 - alb) * weather$swdown[hr] +
#           em * weather$lwdown[hr] # absorbed radiation (W m^-2)
# 
#   # ----------- Resistance to heat transfer and friction velocity ------- #
#   res <- SnowResistance(depth,
#                         hgt = 0, # vegetation height (m)
#                         pai = 0, # plant area index
#                         uref = weather$windspeed[hr], # wind speed (m s^-1)
#                         zref = 2, # reference height (m)
#                         LL = LL, # Obukhov length (m)
#                         z0 = z0)
# 
#   # ------------- Solve coupled snow-soil heat transfer ---------------- #
#   snowsoil <- SnowSoilHeatModel(snowmod,
#                                 soilheat,
#                                 soilp,
#                                 Rabs,
#                                 Tref = weather$temp[hr], # air temperature (°C)
#                                 relhum = weather$relhum[hr], # relative humidity (%)
#                                 atmPressure = weather$pres[hr], # pressure (kPa)
#                                 rHa = res$rHa)
#   snowmod <- snowsoil$snowmod
#   soilheat <- snowsoil$soilheatmod
# 
#   # Surface temperature: that of the soil, or of snow no warmer than 0 °C
#   if (depth > 0) {
#     Tsurface <- min(snowmod$temp[1], 0)
#   } else {
#     Tsurface <- soilheat$Te[1]
#   }
# 
#   # --------------------- Compute atmospheric stability ---------------- #
#   ph <- rhoair(weather$temp[hr],
#                weather$pres[hr]) # air density (mol m^-3)
#   cp <- cpair(weather$temp[hr]) # heat capacity of air
#   H <- (ph * cp / res$rHa) *
#        (Tsurface - weather$temp[hr]) # sensible heat flux (W m^-2)
#   LL <- Obukhov(weather$temp[hr],
#                 weather$pres[hr],
#                 res$uf,
#                 H) # Obukhov length, used in the next hour
# 
#   # ----------------------- Run snow mass balance ---------------------- #
#   forcing <- list(Tair = weather$temp[hr], # air temperature (°C)
#                   relhum = weather$relhum[hr], # relative humidity (%)
#                   pk = weather$pres[hr], # pressure (kPa)
#                   rHa = res$rHa, # resistance to heat transfer (s m^-1)
#                   precip = weather$precip[hr]) # precipitation (mm h^-1)
#   snowmod <- SnowMassBalance(snowmod,
#                              forcing,
#                              snowenv = "Maritime")
# 
#   # -------------------- Regularise snow layer structure ---------------- #
#   snowmod <- RegulariseSnowLayers(snowmod)
# 
#   # ---------------------------- Store outputs -------------------------- #
#   snowdepth[hr] <- sum(snowmod$depth) # snow depth (m)
#   surfacetemp[hr] <- Tsurface # temperature of the snow or soil surface (°C)
#   soiltemp[hr] <- soilheat$Te[1] # temperature of the soil surface (°C)
#   snowfall[hr] <- snowmod$snowfall # snow added (mm)
#   melt[hr] <- snowmod$melt + snowmod$rainmelt # snow melted (mm)
#   drainage[hr] <- snowmod$throughfall # water leaving the base of the snow (mm)
# 
#   # State at the end of the hour becomes that at the start of the next
#   snowmod$oldtemp <- snowmod$temp
#   soilheat$oldTe <- soilheat$Te
# }
# 
# # ----------------------------- Plot results ------------------------------ #
# par(mfrow = c(3, 1))
# tmec <- as.POSIXct(tme, tz = "UTC")
# 
# # Snow depth
# plot(I(snowdepth * 100) ~ tmec, type = "l",
#      col = "steelblue",
#      lwd = 2,
#      xlab = "Date",
#      ylab = "Snow depth (cm)")
# 
# # Temperature of the air, the surface and the soil surface
# plot(weather$temp ~ tmec, type = "l",
#      col = "grey",
#      ylim = c(-16, 6),
#      xlab = "Date",
#      ylab = expression("Temperature (" * degree * "C)"))
# lines(surfacetemp ~ tmec, col = "steelblue")
# lines(soiltemp ~ tmec, col = "firebrick", lwd = 2)
# legend("bottomright",
#        legend = c("Air", "Surface", "Soil surface"),
#        col = c("grey", "steelblue", "firebrick"),
#        lwd = c(1, 1, 2),
#        bty = "n")
# 
# # Water added as snow, melted, and leaving the base of the snow
# plot(cumsum(snowfall) ~ tmec, type = "l",
#      col = "steelblue",
#      lwd = 2,
#      ylim = c(0, 65),
#      xlab = "Date",
#      ylab = "Cumulative water (mm)")
# lines(cumsum(melt) ~ tmec, col = "firebrick", lwd = 2)
# lines(cumsum(drainage) ~ tmec, col = "black", lwd = 2, lty = 2)
# legend("topleft",
#        legend = c("Snowfall", "Melt", "Drainage"),
#        col = c("steelblue", "firebrick", "black"),
#        lty = c(1, 1, 2),
#        lwd = 2,
#        bty = "n")

# # Vary canopy plant area index
# pai_seq <- seq(0.1, 6, length.out = 100) # plant area index
# 
# # Compute snow caught by the canopy (mm SWE)
# interception <- CanopySnowInterception(hgt = 8, # canopy height (m)
#                                        pai = pai_seq,
#                                        uf = 0.3, # friction velocity (m s^-1)
#                                        prec = 5, # snowfall amount (mm SWE)
#                                        tc = -5, # air temperature (°C)
#                                        Li = 0)$intercepted # no snow held to begin with
# 
# # Convert to fraction of snowfall intercepted
# frac_intercepted <- interception / 5
# 
# # ----------------------------- Plot results ------------------------------ #
# par(mfrow = c(1, 1))
# plot(pai_seq, frac_intercepted,
#      type = "l",
#      lwd = 2,
#      ylim = c(0, 1),
#      xlab = "Plant area index (PAI)",
#      ylab = "Interception fraction")

# # Calculate solar zenith and azimuth as function of time of day
# tme <- as.POSIXlt((0:1439) * 60,
#                   origin = "2025-06-21 00:00",
#                   tz = "UTC")
# solp <- sunposition(tme,
#                     lat = 50,
#                     long = 0)
# # Select daylight hours only
# sel <- which(solp$zen < 90)
# tme <- tme[sel]
# solp <- sunposition(tme,
#                     lat = 50,
#                     long = 0)
# 
# # Calculate fraction of direct radiation intercepted by an upright
# # human facing south
# height <- 170 # height (cm)
# width <- 45 # width (cm)
# len <- 25 # body depth (cm)
# svals <- silhouette(solp$zen,
#                     solp$azi,
#                     height,
#                     width,
#                     len,
#                     adir = 180, # facing south (degrees)
#                     atilt = 0.0, # upright posture
#                     position = "fixed")
# 
# # ----------------------------- Plot results ------------------------------ #
# par(mfrow = c(1, 1))
# plot(svals$SolarCoef ~ as.POSIXct(tme),
#      type = "l",
#      col = "firebrick",
#      lwd = 2,
#      ylim = c(0, 0.5),
#      xlab = "Hour",
#      ylab = "Fraction of direct beam radiation intercepted")

# # Animal geometry (cm)
# height <- 12 # height (cm)
# width <- 8 # width (cm)
# len <- 25 # length (cm)
# 
# # Orientation relative to wind direction
# adir <- seq(0, 180, by = 1) # orientation angle (degrees)
# 
# # Characteristic dimension (m)
# d <- chardim(wdir = 0, # wind direction (degrees)
#              height,
#              width,
#              len,
#              adir,
#              atilt = 0, # body tilt angle (degrees)
#              position = "fixed")
# 
# # Boundary layer resistance at fixed conditions
# rHa <- sapply(d, function(dd) {
#   animalrHa(tair = 20, # air temperature (°C)
#             dT = 5, # body-air temperature difference (°C)
#             uz = 1, # wind speed (m s^-1)
#             d = dd) # characteristic dimension (m)
# })
# 
# # ----------------------------- Plot results ------------------------------ #
# par(mfrow = c(2, 1),
#     mar = c(5, 5, 2, 2))
# plot((d * 100) ~ adir,
#      type = "l",
#      lwd = 2,
#      ylim = c(0, 25),
#      xlab = "Orientation relative to wind direction (degrees)",
#      ylab = "Characteristic dimension (cm)")
# plot(rHa ~ adir,
#      type = "l",
#      lwd = 2,
#      ylim = c(0, 150),
#      xlab = "Orientation relative to wind direction (degrees)",
#      ylab = expression("Boundary layer resistance (s m"^-1 * ")"))

# require(micropoint) # ensure micropoint package is installed
# # Animal geometry (cm)
# height <- 12 # height (cm)
# width <- 8 # width (cm)
# len <- 25 # length (cm)
# 
# # Create ectotherm parameter set
# ecto_params <- micropoint::create_ectop(height,
#                                         width,
#                                         len,
#                                         refl = 0.1, # body reflectance
#                                         adir = 0, # orientation angle (degrees)
#                                         position = "randomdir", # random orientation
#                                         type = "reptile") # ectotherm type
# 
# # Assign variable body temperature
# Tbody <- seq(-5, 60, length.out = 100) # body temperature (°C)
# 
# # Compute metabolic rate
# met_rate <- with(ecto_params,
#                  metabolic_rate(volume, # body volume (m^3)
#                                 rho, # body density (kg m^-3)
#                                 Q10, # temperature sensitivity coefficient
#                                 a0, # allometric intercept
#                                 b, # allometric exponent
#                                 Tref, # reference temperature (°C)
#                                 Tbody))
# 
# # ----------------------------- Plot results ------------------------------ #
# par(mfrow = c(1, 1))
# plot(met_rate ~ Tbody,
#      type = "l",
#      col = "firebrick",
#      lwd = 2,
#      xlab = expression("Body temperature (" * degree * "C)"),
#      ylab = expression("Metabolic rate (W)"))

# # -------------------- Create June weather dataset ------------------------- #
# # Reload the weather dataset: the name climdata is reused for model driving variables
# climdata <- microclimlearn::climdata
# weather <- climdata
# 
# # Select data for June
# tme <- as.POSIXlt(weather$obs_time, tz = "UTC")
# sel <- which(tme$mon + 1 == 6)
# weather <- weather[sel,]
# tme <- tme[sel]
# 
# # -------------------- Create vegetation parameter set -------------------- #
# # Bare ground
# vegp <- createplant_inputs("C3")
# vegp$h <- 0     # vegetation height (m): zero for bare ground
# vegp$pai <- 0   # plant area index: zero for bare ground
# 
# # -------------------- Create soil datasets and initialise models --------- #
# groundp <- micropoint::groundparams
# 
# # Lower soil boundary temperature
# boundaryT <- 11.36071
# 
# soilp <- createsoilplist(
#   soiltype = "Clay loam",
#   nlayers = 15,
#   totalDepth = 2,
#   slope = 0,
#   aspect = 180,
#   surface_organicmu = 1,
#   FreeDrain = TRUE
# )
# 
# # Number of soil layers
# n <- length(soilp$b)
# 
# # Initial soil temperature profile
# Te <- seq(weather$temp[1], boundaryT, length.out = n)
# 
# # Initial soil water content profile
# theta <- seq(0.2, soilp$Smax[n], length.out = n)
# 
# # Initialise soil heat model
# soilheat <- InitSoilheatmod(soilp, Te, theta)
# 
# # -------------------- Create ectotherm parameter set --------------------- #
# ecto_params <- micropoint::create_ectop(
#   height = 12,
#   width = 8,
#   len = 25,
#   refl = 0.1,
#   adir = 0,
#   position = "randomdir",
#   type = "reptile"
# )
# 
# # -------------------- Run model in hourly timesteps ---------------------- #
# Tbody <- G <- Tsurface <- 0
# 
# for (hr in 1:720) {
# 
#   # Initialise convergence criterion
#   error <- 1e99
#   itr <- 1
# 
#   # Solve the ground surface and the soil alternately until the ground heat
#   # flux settles
#   while (error > 1e-2 && itr < 20) {
# 
#     # -------------------- Run microclimate model ------------------------ #
# 
#     # Estimate relative humidity at soil surface (percentage)
#     soilrh <- soilrelhum(
#       0.2,              # volumetric water content at the surface
#       soilheat$Te[1],   # soil surface temperature (°C)
#       soilp$psi_e[1],   # air-entry water potential (J/kg)
#       soilp$b[1],       # Campbell shape parameter
#       soilp$Smax[1]     # water content at saturation
#     )
# 
#     # Solve the energy balance of the ground surface
#     can <- solve_wholecanopy(
#       weather,
#       vegp,
#       groundp,
#       hr,
#       49.96807,
#       -5.215668,
#       zref = 2,
#       G,
#       theta = 0.2,
#       soilrh = soilrh
#     )
# 
#     # -------------------- Solve soil heat transport --------------------- #
# 
#     soilheat <- SoilHeatModel(
#       soilheat,
#       soilp,
#       Rabs = can$Rabs,               # radiation absorbed by the ground
#       Tref = weather$temp[hr],
#       relhum = weather$relhum[hr],
#       atmPressure = weather$pres[hr],
#       rHa = can$rHg                  # resistance from ground to zref
#     )
# 
#     # Update convergence criterion
#     error <- abs(G - soilheat$Gflux)
#     G <- soilheat$Gflux
#     itr <- itr + 1
#   }
# 
#   # Move the soil model on to the next time step
#   soilheat$oldTe <- soilheat$Te
# 
#   # --------- Calculate conditions 6 cm above the ground surface --------- #
# 
#   # Resistance to heat transfer between the ground surface and 6 cm
#   r6 <- canopyresistance(0.06, can$uf, vegp$h, vegp$pai, can$LL)
# 
#   # Air temperature at 6 cm: the share of the temperature difference between
#   # ground and reference height that lies above 6 cm
#   Tz <- weather$temp[hr] +
#     (can$tcanopy - weather$temp[hr]) *
#     ((can$rHa - r6) / can$rHa)
# 
#   # Vapour pressures
#   es <- satvap(can$tcanopy)
#   ea <- satvap(weather$temp[hr]) *
#     (weather$relhum[hr] / 100)
# 
#   # Vapour pressure at 6 cm, where vapour from the ground also meets the
#   # resistance of the soil surface
#   ez <- ea +
#     (es - ea) *
#     (can$rHa - r6) / (can$rHa + can$rS)
# 
#   # Relative humidity at 6 cm
#   Rhz <- pmin((ez / satvap(Tz)) * 100, 100)
# 
#   # Wind speed at 6 cm
#   uz <- windprofile_above(0.06, can$uf, vegp$h, vegp$pai, can$LL)
# 
#   # --------- Calculate upward and downward radiation fluxes ------------ #
# 
#   # Solar geometry
#   solp <- sunposition(
#     tme[hr],
#     lat = 49.96807,
#     long = -5.215668
#   )
# 
#   # Surface albedo
#   albs <- albedo(
#     vegp$pai,
#     vegp$lref,
#     vegp$ltra,
#     vegp$x,
#     groundp$gref,
#     groundp$slope,
#     groundp$aspect,
#     solp
#   )
# 
#   # Compute direct and reflected shortwave radiation
#   if (weather$swdown[hr] > 0) {
# 
#     Rdirdown <- (weather$swdown[hr] -
#       weather$difrad[hr]) /
#       cos(solp$zen * pi / 180)
# 
#     Rswup <- albs$albd * weather$difrad[hr] +
#       Rdirdown * cos(solp$zen * pi / 180) *
#       albs$albb
# 
#   } else {
# 
#     Rdirdown <- Rswup <- 0
#   }
# 
#   # ---------------- Calculate ectotherm body temperature --------------- #
# 
#   microcdata <- list(
#     Tair = Tz,
#     Tsurf = can$tcanopy,
#     relhum = Rhz,
#     surfrelhum = soilrh,
#     pk = weather$pres[hr],
#     windspeed = uz,
#     winddir = weather$winddir[hr],
#     Rdirdown = Rdirdown,
#     Rdifdown = weather$difrad[hr],
#     Rswup = Rswup,
#     Rlwdown = weather$lwdown[hr],
#     # longwave emitted by the ground and sky longwave it reflects
#     Rlwup = groundp$em * 5.67e-8 * (can$tcanopy + 273.15)^4 +
#       (1 - groundp$em) * weather$lwdown[hr]
#   )
# 
#   # Solve ectotherm heat balance
#   Tbody[hr] <- Ectothermtemp(
#     ecto_params,
#     tme[hr],
#     microcdata,
#     lat = 49.96807,
#     long = -5.215668
#   )
# }
# 
# # ------------------------- Plot temperatures ----------------------------- #
# tmec <- as.POSIXct(tme)
# 
# # Air temperature
# plot(
#   weather$temp ~ tmec,
#   type = "l",
#   xlab = "Date",
#   ylab = "Temperature (°C)",
#   ylim = c(0, 50),
#   col = rgb(0, 0, 0, 0.5),
#   lwd = 2
# )
# 
# # Ectotherm body temperature
# par(new = TRUE)
# 
# plot(
#   Tbody ~ tmec,
#   type = "l",
#   xlab = "",
#   ylab = "",
#   ylim = c(0, 50),
#   col = "firebrick"
# )
# 
# # Add legend
# legend(
#   "topleft",
#   legend = c("Air temperature", "Body temperature"),
#   col = c(rgb(0, 0, 0, 0.5), "firebrick"),
#   lwd = 2,
#   bty = "n"
# )

# # ----------------------- Create endotherm parameters ----------------------- #
# lion <- create_endop(height = 55,    # body height (cm)
#                      width = 40,     # body width (cm)
#                      len = 170,      # body length (cm)
#                      coat = 0.5,     # coat depth (cm)
#                      type = "mammal",
#                      refl = 0.3,     # reflectance of the coat
#                      adir = 90,      # facing east
#                      position = "fixed")
# lion$Tb <- 38.5 # core temperature of a lion (°C)
# lion$Tcmax <- lion$Tb + 5 # highest core temperature it can survive (°C)

# # -------------------------- Select June weather ---------------------------- #
# # Reload the weather dataset: the name climdata is reused for model driving variables
# climdata <- microclimlearn::climdata
# weather <- climdata
# # Extract time and subset to June
# tme <- as.POSIXlt(weather$obs_time, tz = "UTC")
# sel <- which(tme$mon + 1 == 6)
# weather <- weather[sel, ]
# tme <- tme[sel]
# 
# # -------------------------- Helper functions ------------------------------- #
# # Logit transforms used to adjust bounded variables (0–1)
# logit <- function(x) log(x / (1 - x))
# inv_logit <- function(x) 1 / (1 + exp(-x))
# 
# # -------------------------- Adjust temperature & humidity ------------------ #
# # Increase temperature to represent warmer (e.g. tropical) conditions
# temp_new <- weather$temp + 12
# # Recalculate relative humidity assuming constant vapour pressure
# ea <- satvap(weather$temp) * (weather$relhum / 100)
# relhum_new <- (ea / satvap(temp_new)) * 100
# 
# # -------------------------- Adjust longwave radiation ---------------------- #
# # Derive sky emissivity from observed longwave radiation
# skyem <- weather$lwdown / (5.67e-8 * (weather$temp + 273.15)^4)
# # Reduce emissivity slightly (clearer skies → lower emissivity)
# skyem_new <- inv_logit(logit(skyem) - 0.3)
# 
# # Recalculate downward longwave radiation using updated emissivity and temperature
# lwdown_new <- skyem_new * 5.67e-8 * (temp_new + 273.15)^4
# 
# # Assemble updated weather object
# weather_new <- weather
# weather_new$temp <- temp_new
# weather_new$relhum <- relhum_new
# weather_new$lwdown <- lwdown_new
# weather_new$windspeed <- 0.5 * weather$windspeed
# 
# # -------------------------- Adjust shortwave radiation --------------------- #
# # Compute clear-sky radiation for original conditions
# csr <- clearskyrad(weather, 49.96807, -5.215668)
# 
# # Calculate clear-sky fraction
# csf <- weather$swdown / csr
# 
# # Increase clearness (make conditions sunnier) using logit transform
# csf_new <- inv_logit(logit(csf) + 0.3)
# 
# # Recompute clear-sky radiation for new (warmer) conditions and location
# csr_new <- clearskyrad(weather_new, -1.21867, -5.215668)
# 
# # Apply adjusted clear-sky fraction
# weather_new$swdown <- csf_new * csr_new
# 
# # Ensure nighttime and invalid values are handled
# weather_new$swdown[weather$swdown == 0] <- 0
# weather_new$swdown[is.na(weather_new$swdown)] <- 0
# 
# # -------------------------- Partition into diffuse radiation --------------- #
# # Compute diffuse fraction based on solar geometry and atmospheric conditions
# diffrac <- difprop(weather_new, -1.21867, -5.215668)
# 
# # Calculate diffuse shortwave radiation
# weather_new$difrad <- diffrac * weather_new$swdown

# # -------------------- Ground and soil ---------------------------------------- #
# vegp <- createplant_inputs("C3")
# vegp$h <- 0     # bare ground
# vegp$pai <- 0
# groundp <- micropoint::groundparams
# lat <- -1.21867
# long <- -5.215668
# soilp <- createsoilplist(soiltype = "Clay loam", nlayers = 15, totalDepth = 2,
#                          slope = 0, aspect = 180, surface_organicmu = 1,
#                          FreeDrain = TRUE)
# n <- length(soilp$b)
# Te <- seq(weather_new$temp[1], 23.3, length.out = n) # 23.3 °C at the lower boundary
# # Dry savanna soil: air-dry at the surface, moister at depth
# theta <- seq(0.05, soilp$Smax[n], length.out = n)
# soilheat <- InitSoilheatmod(soilp, Te, theta)
# 
# # ----------------------- Run model in hourly timesteps ---------------------- #
# micro <- list() # conditions around the lion for each hour
# Tz_all <- G <- 0
# for (hr in 1:720) {
#   error <- 1e99
#   itr <- 1
#   while (error > 1e-2 && itr < 20) {
#     soilrh <- soilrelhum(0.05, soilheat$Te[1], soilp$psi_e[1], soilp$b[1], soilp$Smax[1])
#     can <- solve_wholecanopy(weather_new, vegp, groundp, hr, lat, long, zref = 2,
#                              G, theta = 0.05, soilrh = soilrh)
#     soilheat <- SoilHeatModel(soilheat, soilp, Rabs = can$Rabs,
#                               Tref = weather_new$temp[hr],
#                               relhum = weather_new$relhum[hr],
#                               atmPressure = weather_new$pres[hr], rHa = can$rHg)
#     error <- abs(G - soilheat$Gflux)
#     G <- soilheat$Gflux
#     itr <- itr + 1
#   }
#   soilheat$oldTe <- soilheat$Te
# 
#   # ------------- Air temperature, humidity and wind at 30 cm --------------- #
#   r30 <- canopyresistance(0.3, can$uf, vegp$h, vegp$pai, can$LL)
#   Tz <- weather_new$temp[hr] +
#     (can$tcanopy - weather_new$temp[hr]) * (can$rHa - r30) / can$rHa
#   # Vapour evaporating from the ground raises the vapour pressure near the
#   # surface, and its excess over that at 2 m falls off with height as heat does
#   ea <- satvap(weather_new$temp[hr]) * (weather_new$relhum[hr] / 100)
#   ez <- ea + (satvap(can$tcanopy) - ea) * (can$rHa - r30) / (can$rS + can$rHa)
#   uz <- windprofile_above(0.3, can$uf, vegp$h, vegp$pai, can$LL)
# 
#   # ------------------- Direct and reflected shortwave ---------------------- #
#   solp <- sunposition(tme[hr], lat, long)
#   albs <- albedo(vegp$pai, vegp$lref, vegp$ltra, vegp$x, groundp$gref,
#                  groundp$slope, groundp$aspect, solp)
#   if (weather_new$swdown[hr] > 0) {
#     Rdirdown <- (weather_new$swdown[hr] - weather_new$difrad[hr]) /
#       cos(solp$zen * pi / 180)
#     Rswup <- albs$albd * weather_new$difrad[hr] +
#       Rdirdown * cos(solp$zen * pi / 180) * albs$albb
#   } else {
#     Rdirdown <- Rswup <- 0
#   }
# 
#   # ------------- Conditions around the lion, stored for each hour ---------- #
#   microcdata <- list(Tair = Tz,
#                      Tsurf = can$tcanopy, # ground surface temperature
#                      relhum = pmin(100 * ez / satvap(Tz), 100),
#                      pk = weather_new$pres[hr],
#                      windspeed = uz,
#                      winddir = weather_new$winddir[hr],
#                      Rdirdown = Rdirdown,
#                      Rdifdown = weather_new$difrad[hr],
#                      Rswup = Rswup,
#                      Rlwdown = weather_new$lwdown[hr],
#                      Rlwup = groundp$em * 5.67e-8 * (can$tcanopy + 273.15)^4 +
#                        (1 - groundp$em) * weather_new$lwdown[hr])
#   micro[[hr]] <- microcdata
#   Tz_all[hr] <- Tz
# }

# # ----------------------- Lion in sun and in shade -------------------------- #
# E_sun <- E_shade <- Emax <- numeric(720)
# for (hr in 1:720) {
#   shadedata <- micro[[hr]]
#   shadedata$Rdirdown <- 0 # shade: direct beam blocked
#   sun <- Endotherm(lion, tme[hr], micro[[hr]], lat, long, evapunits = "L/day")
#   shade <- Endotherm(lion, tme[hr], shadedata, lat, long, evapunits = "L/day")
#   E_sun[hr] <- sun$Lreq
#   E_shade[hr] <- shade$Lreq
#   Emax[hr] <- shade$Lmax
# }
# 
# # ----------------------------- Plot results ------------------------------ #
# tmec <- as.POSIXct(tme)
# plot(E_sun ~ tmec, type = "l", col = "firebrick", ylim = c(0, 20),
#      xlab = "Date", ylab = expression("Evaporation (L day"^-1 * ")"))
# lines(E_shade ~ tmec, col = "steelblue")
# lines(Emax ~ tmec, lty = 2) # most the lion can evaporate
# legend("topleft", legend = c("Required, sun", "Required, shade", "Maximum"),
#        col = c("firebrick", "steelblue", "black"), lty = c(1, 1, 2), bty = "n", cex = 0.8)

# # ----------------------- Lion in sun and in shade -------------------------- #
# PET_sun <- PET_shade <- Tc_sun <- Tc_shade <- numeric(720)
# for (hr in 1:720) {
#   shadedata <- micro[[hr]]
#   shadedata$Rdirdown <- 0 # shade: direct beam blocked
#   sun <- EndothermPET(lion, tme[hr], micro[[hr]], lat, long)
#   shade <- EndothermPET(lion, tme[hr], shadedata, lat, long)
#   PET_sun[hr] <- sun$PET
#   PET_shade[hr] <- shade$PET
#   Tc_sun[hr] <- sun$Tcore
#   Tc_shade[hr] <- shade$Tcore
# }
# 
# # ----------------------------- Plot results ------------------------------ #
# par(mfrow = c(2, 1))
# plot(PET_sun ~ tmec, type = "l", col = "firebrick", ylim = c(10, 80),
#      xlab = "Date", ylab = expression("Temperature (" * degree * "C)"))
# lines(PET_shade ~ tmec, col = "steelblue")
# lines(Tz_all ~ tmec, col = "darkgrey")
# legend("topleft", legend = c("Equivalent temperature, sun",
#                               "Equivalent temperature, shade", "Air temperature"),
#        col = c("firebrick", "steelblue", "darkgrey"), lwd = 1, bty = "n", cex = 0.8)
# plot(Tc_sun ~ tmec, type = "l", col = "firebrick", ylim = c(15, 60),
#      xlab = "Date", ylab = expression("Core temperature (" * degree * "C)"))
# lines(Tc_shade ~ tmec, col = "steelblue")
# abline(h = c(lion$Tcmin, lion$Tcmax), lty = 2) # survivable limits

# require(devtools) # Install from CRAN if not already installed
# devtools::install_github("ilyamaclean/micropoint")
# library(micropoint)
# vignette("new_model_tutorial", package = "micropoint")

# require(devtools) # Install from CRAN if not already installed
# devtools::install_github("ilyamaclean/microclimf", ref = "dev")
# library(microclimf)

# vignette("running-microclimf", package = "microclimf")

# require(devtools) # Install from CRAN if not already installed
# devtools::install_github("ilyamaclean/mesoclim")
# library(mesoclim)

# remotes::install_github(“dklinges9/mcera5”
# library(mcera5)

# require(devtools) # Install from CRAN if not already installed
# devtools::install_github("ilyamaclean/microclimdata")
# library(microclimdata)

# install.packages("TrenchR")
# library(TrenchR)
