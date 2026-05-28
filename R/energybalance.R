#' Calculate energy balance
#'
#' @param Ta Air temperature (degrees C)
#' @param Ts temperature of surface for which energy balance is computed (degrees C)
#' @param pk atmospheric pressure (kPa)
#' @param rh relative humidity (percentage)
#' @param Rabs flux density of absorbed radiation (W/m^2)
#' @param rHa resistance to heat loss (s/m)
#' @param rS stomatal resistance (s/m)
#' @param G rate of heat storage (W/m^2)
#' @param em emissivity of surface (unitless)
#' @returns energy balance (W/m^2).
#' @export
energybalance <- function(Ta, Ts, pk, rh, Rabs, rHa, rS, G = 0, em = 0.97) {
  # Calculate emitted radiation
  sb <-  5.67 *10 ^-8 # Stefan-Boltzmann constant
  Rem <- em * sb * (Ts + 273.15) ^ 4 # W/m^2
  # Calculate sensible heat flux
  cp <- cpair(Ta)  # specific heat of air
  ph <- rhohair(Ta, pk)   # molar density fo air
  H <- ((ph * cp) / rHa) * (Ts - Ta)
  # Calculate latent heat flux
  la <- latvap(Ts) # Latent heat of vapourization
  rV <- rS + rHa   # total resistance to vapour loss
  es <- satvap(Ts) # vapour pressure of leaf (kPa)
  ea <- satvap(Ta) * (rh / 100) # vapour pressur eof air (kPa)
  L <- ((la * ph) / (rV * pk)) * (es - ea)
  # Calculate energy balance
  Ba <- Rabs - Rem - H - L - G
  return(Ba)
}
#' Calculate temperature of surface using Penman-Monteith equation
#'
#' @param Ta Air temperature (degrees C)
#' @param pk atmospheric pressure (kPa)
#' @param rh relative humidity (percentage)
#' @param Rabs flux density of absorbed radiation (W/m^2)
#' @param rHa resistance to heat loss (s/m)
#' @param rS stomatal resistance (s/m)
#' @param G rate of heat storage (W/m^2)
#' @param em emissivity of surface (unitless)
#' @param iters integer specifying how many times to iterate the equation (4 achieves near exact solution)
#' @returns temperature of surface (deg C)
#' @export
PenmanMonteith <- function(Ta, pk, rh, Rabs, rHa, rS, G, em = 0.97, iters = 4) {
  # Initially set surface temperature to air temperature
  Ts <- Ta
  sb = 5.67e-8 # Stefan-Boltzmann constant
  cp <- cpair(Ta)  # specific heat of air
  ph <- rhohair(Ta, pk)   # molar density fo air
  Rema <- em * sb * (Ta + 273.15)^4
  rV = rHa + rS
  for (i in 1:iters) {
    Te <- (Ts + Ta) / 2.0
    Rer <- 4.0 * em * sb * (Te + 273.15) ^ 3  # radiative conductance
    la <- latvap(Ts) # Latent heat of vapourization
    Da <- satvap(Ta) * (1.0 - rh / 100.0) # Vapour pressure deficit
    De <- satvap(Te + 0.5) - satvap(Te - 0.5) # slope of saturated vapour pressure curve
    Ts <- Ta + ((Rabs - Rema - ((la * ph) / (pk * rV)) * Da - G) /
                  (Rer + ph * (cp / rHa + ((la * De) / (pk * rV)))))
  }
  return (Ts)
}
#' Solve energy balance equation to derive surface temperature
#'
#' This function solves for the surface (or organism) temperature (\eqn{T_s})
#' that satisfies the steady-state energy balance, based on various meteorological
#' and physical parameters. The energy balance can be solved using either the
#' **uniroot method** (slower) or the **Penman-Monteith method** (faster).
#' Users also have the option to iterate the Penman method for a more accurate solution.
#'
#' @param Ta Air temperature (degrees Celsius).
#' @param pk Atmospheric pressure (kPa).
#' @param rh Relative humidity (percentage).
#' @param Rabs Flux density of absorbed radiation (W/m^2).
#' @param rHa Resistance to heat loss (s/m).
#' @param rS Stomatal resistance (s/m).
#' @param G Rate of heat storage (W/m^2).
#' @param em Emissivity of surface (unitless). Default is 0.97.
#' @param method One of `'uniroot'` or `'Penman'` (see details for more information on each method).
#' @param iters Integer specifying the number of iterations for the Penman method. Default is 1.
#'
#' @details
#' This function solves for the surface (or organism) temperature (\eqn{T_s})
#' that satisfies the steady-state energy balance, based on various meteorological
#' and physical parameters. The energy balance can be solved using either the
#' **uniroot method** (slower) or the **Penman-Monteith method** (faster).
#' Users also have the option to iterate the Penman method for a more accurate solution.
#' @returns A numeric vector of surface temperatures (in degrees Celsius) corresponding to each input
#'          air temperature (`Ta`).
#' @export
SolveEnergyBalance <- function(Ta, pk, rh, Rabs, rHa, rS, G, em = 0.97,
                               method = c("uniroot", "Penman"), iters = 1) {
  method <- match.arg(method)
  if (method == "uniroot") {
    Ts <- 0
    # One search interval that’s wide enough for all Ta
    mn <- min(Ta) - 30
    mx <- max(Ta) + 50
    for (hr in 1:length(Ta)) {
        Ts[hr] <- uniroot(
          function(Ts_est) energybalance(Ta[hr], Ts_est, pk[hr], rh[hr], Rabs[hr], rHa[hr], rS[hr], G[hr]),
          interval = c(mn, mx)
        )$root
    }
  } else { # method == "Penman"
    Ts <- PenmanMonteith(Ta, pk, rh, Rabs, rHa, rS, G, em, iters)
  }
  return(Ts)
}
#' Solve energy balance equation for entire canopy returning associated variables
#' @param weather a data.frame of weather variables matching the format of `micropoint::climdata`
#' @param vegp a list of vegetation parameters as returned by [createplant_inputs()]
#' @param groundp a list of ground parameters matching the format of `micropoint::groundparams`
#' @param hr hour in weather for energy balance is to be solved
#' @param lat latitude (decimal degrees)
#' @param long longitude (decimal degrees)
#' @param zref height above ground of weather variables (must be above canopy - see details)
#' @param G Rate of heat storage by ground at `hr` (W/m^2)
#' @param theta soil water fraction in root zone
#' @param soilrh effective soil relative humidity (percentage)
#' @param Ca  Atmospheric CO2 concentration in ppm as returned by `micropoint::Cafromyear()`
#' @param maxiter  Maxiumum number of iterations to achieve convergence
#' @details Meteorological inputs must represent above-canopy conditions,
#' whereas forcing data are often measured 1–2 m above the ground, below the
#' height of surrounding vegetation (e.g. forest canopy). The function
#' `weatherhgt_adjust` in the `micropoint` package  can be used to adjusts
#' weather variables to estimate the values that would occur sufficiently high
#' above the surface, assuming site characteristics consistent with WMO guidelines
#' for weather station siting.
#' @return a list of of the following variables
#' \describe{
#'  \item{tcanopy}{temperature of canopy heat exchange surface (deg C)}
#'  \item{Rabs}{Flux density of radiation absorbed by canopy (W/m^2)}
#'  \item{H}{Flux density of sensible heat exchanged between canopy and `zref` (W/m^2)}
#'  \item{L}{Flux density of latent heat exchanged between canopy and `zref` (W/m^2)}
#'  \item{uf}{wind friction velocity (m/s)}
#'  \item{uh}{wind speed at top of canopy (m/s)}
#'  \item{th}{temperature at top of canopy (deg C)}
#'  \item{rHa}{bulk surface aerodynamic resistance (s/m)}
#'  \item{rSl}{bulk surface stomatal resistance (s/m)}
#'  \item{rS}{bulk surface resistance to vapour including from ground surface (s/m)}
#'  \item{psim}{diabatic influencing factor for momentum}
#'  \item{psih}{diabatic influencing factor for heat}
#'  \item{LL}{Obukhov length}
#'  \item{TL}{Langrangian time-scale}
#' }
#' @rdname solve_wholecanopy
#' @export
solve_wholecanopy <- function(weather, vegp, groundp, hr, lat, long, zref, G, theta, soilrh, Ca = 430, maxiter = 20) {
  if (zref < vegp$h) stop("zref must be above canopy\n")
  # Calculate radiation absorbed by canopy
  Rabs <- calcRabs(weather, vegp, groundp,  lat, long)[hr]
  # Calculate solar position
  tme <- as.POSIXlt(weather$obs_time[hr], tz = "UTC")
  solp <- sunposition(tme, lat, long)
  # extract relevent climate data
  uref <- weather$windspeed[hr]
  Rsw <- weather$swdown[hr]
  Rdif <- weather$difrad[hr]
  Ta <- weather$temp[hr]
  pk <- weather$pres[hr]
  rh <- weather$relhum[hr]
  # initially set diabatic coeffcients to 0
  psih <- psim <- H <- 0
  LL <- 1e10
  # initialize canopy temperature
  tcanopy <- Ta + 0.01 * Rabs
  # Calculate zero-plane displacement height
  d <- zeroplanedis(vegp$h, vegp$pai)
  # Calculate soil water potential in root zone
  psi_r <- PsiFromtheta(theta, groundp$Smax, groundp$Psie, groundp$b)
  error <- 1e99
  itr <- 1
  while (error > 1e-3 && itr < maxiter) {
    # Calculate bulk surface aerodynamic resistance
    zm <- roughlength(vegp$h, vegp$pai, d, psih)  # roughness length for momentum
    uf <- windfric(uref, zref, d, zm, psim) # friction velocity
    rHa <- canopyresistance(zref, d, zm, uf, psih) # resistance to heat loss
    # Calculate bulk surface stomatal resistance
    rSl <- bulkstomatalresist_calc(solp, Ca, Rsw, Rdif, Ta, tcanopy, rh,
                                   pk, psi_r, vegp)
    rSl[rSl > 1e10] <- 1e10
    # Include the ground component  and scale by LAIfrac
    rG <- (groundresistance(vegp$h, vegp$pai, uf, LL, psih, zref) * 100) / soilrh
    rS <- 1 / (vegp$LAIfrac / rSl + 1/rG)
    tcanopy <- PenmanMonteith(Ta, pk, rh, Rabs, rHa, rS, G, vegp$em)
    # cap at dewpoint
    tdew <- dewpoint(Ta, rh)
    if (tcanopy < tdew) tcanopy <- tdew
    # Calculate H
    phcp <- rhohair(Ta, pk)  * cpair(Ta) # Specific heat and density of air
    oldH <- H
    H <-  (phcp / rHa) * (tcanopy - Ta) # Sensible heat flux
    error <- abs(H - oldH)
    # Compute diabatic coefficients
    LL <- Obukhov(Ta, pk, uf, H) # Obukhov length
    psim <- dpsim(zm / LL) - dpsim((zref - d) / LL)
    psih <- dpsih((0.2 * zm) / LL) - dpsih((zref - d) / LL)
    itr <- itr + 1
  }
  # Calculate Latent heat flux
  la <- latvap(tcanopy)
  rV <- rS + rHa
  es <- satvap(tcanopy)
  ea <- satvap(Ta) * (rh / 100)
  L <- ((la * rhohair(Ta, pk)) / (rV * pk)) * (es - ea)
  # Calculate wind speed and temperature at top of canopy
  if (zref == vegp$h) {
    uh <- uref
    th <- weather$temp[hr]
  } else {
    psimh <- dpsim(zm / LL) - dpsim((vegp$h - d) / LL)
    psihh <- dpsih(zm / LL) - dpsim((vegp$h - d) / LL)
    uh <- (uf / 0.41) * (log((vegp$h - d)/ zm) + psimh)
    th <- tcanopy - (H / (0.41 *phcp)) * (log((vegp$h - d)/(0.2 * zm)) + psihh)
  }
  # Calculate diabatic influencing factor
  phih <- dphih((zref - d) / LL)
  # Calculate Langrangian timescale
  a2 <- (phih * 0.41 * (1 - d / vegp$h)) / 1.5625
  TL <- a2 * vegp$h / uf
  # Calculate temperature at top of canopy
  return(list(tcanopy = tcanopy, Rabs = Rabs, H = H, L = L, uf = uf, uh = uh,
              th = th, rHa = rHa, rSl = rSl, rS = rS, psim = psim, psih = psih, LL = LL, TL = TL))
}



