#' @title Snow density parameter sets
#' @description Returns empirical parameter values for a snow density evolution
#' function, based on different climatic snow regimes.
#'
#' These parameterisations are intended for use with [snowdenfun()] to estimate
#' bulk snow density as a function of snow depth and age.
#'
#' @param snowenv Character string specifying the snow environment. Options are:
#' `Alpine`, `Maritime`, `Prairie`, `Tundra` or `Taiga`.
#'
#' @return A numeric vector of length 4 containing density parameters:
#' \describe{
#'   \item{densfun[1]}{Upper asymptotic snow density (g cm^-3)}
#'   \item{densfun[2]}{Initial / minimum snow density (g cm^-3)}
#'   \item{densfun[3]}{Depth-dependent densification coefficient (m^-1)}
#'   \item{densfun[4]}{Age-dependent densification coefficient (day^-1)}
#' }
#' @export
SnowDensParams <- function(snowenv = "Maritime") {
  if (snowenv == "Alpine") {
    densfun <- c(0.5975, 0.2237, 0.0012, 0.0038)
  } else if (snowenv == "Maritime") {
    densfun <- c(0.5979, 0.2578, 0.0010, 0.0038)
  } else if (snowenv == "Prairie") {
    densfun <- c(0.5940, 0.2332, 0.0016, 0.0031)
  } else if (snowenv == "Tundra") {
    densfun <- c(0.3630, 0.2425, 0.0029, 0.0049)
  } else if (snowenv == "Taiga") {
    densfun <- c(0.2170, 0.2170, 0.0000, 0.0000)
  } else {
    stop("snowenv not recognised")
  }
  return(densfun)
}


#' @title Bulk snow density from depth and age
#' @description Computes bulk snow density as a function of snow depth and
#' snow age using an empirical exponential formulation.
#' @param densfun Numeric vector of length 4 as returned by [SnowDensParams()].
#' @param snowdepth Snow depth (m).
#' @param snowage Snow age (hours).
#' @return Numeric value of bulk snow density (kg m^-3).
#' @export
SnowDensFun <- function(densfun, snowdepth, snowage) {
  snowdensity <- ((densfun[1] - densfun[2]) *
                    (1 - exp(-densfun[3] * snowdepth -
                               densfun[4] * snowage / 24)) +
                    densfun[2]) * 1000
  return(snowdensity)
}

#' @title Initialise a multilayer snow model
#' @description Initialises a multilayer snowpack from total snow depth, snow
#' age, and snow temperature.
#' @param snowdepth Total snow depth (m).
#' @param snowage Snow age (hours).
#' @param snowtemp Initial snow temperature (\eqn{^\circ}C). This may be a
#' single value, in which case all layers are assigned the same temperature.
#' If a vector is supplied, it must have length equal to the number of snow
#' layers generated.
#' @param snowenv Character string specifying the snow environment passed to
#' [SnowDensParams()]. Options depend on that function, i.e. `Maritime`,
#' `Alpine`, `Prairie`, `Tundra` or `Taiga`.
#'
#' @return A list containing the snow model state:
#' \describe{
#'   \item{SWE}{Snow water equivalent in each layer (mm)}
#'   \item{depth}{Thickness of each snow layer (m)}
#'   \item{rho}{Bulk density of each layer (kg m^-3)}
#'   \item{temp}{Temperature of each layer (deg C)}
#'   \item{oldtemp}{Temperature of each layer in previous time step (same as temp here)}
#'   \item{age}{Age of each layer (hours)}
#'   \item{melt}{Total melt produced over the timestep (mm)}
#'   \item{snowfall}{Snowfall added over the timestep (mm)}
#'   \item{rainfall}{Rainfall added over the timestep (mm)}
#'   \item{throughfall}{Liquid water leaving the base of the snowpack (mm)}
#'   \item{surfaceWater}{Temporary liquid water stored in the snowpack (mm)}
#' }
#' @export
InitSnowModel <- function(snowdepth, snowage, snowtemp, snowenv = "Maritime") {
  if (snowdepth > 0) {
    dz <- pmin(0.01 * 1.6^(0:100), 0.03)
    dz <- dz[cumsum(dz) <= snowdepth]
    dz <- c(dz, snowdepth - sum(dz))
    # cumulative depth to layer bottoms
    cdz <- cumsum(dz)
    # density profile
    sdp <- SnowDensParams(snowenv)
    rho <- SnowDensFun(sdp, cdz, snowage)
    # layer SWE (mm)
    SWE <- rho * dz
    # layer temperature
    if (length(snowtemp) == 1) {
      snowtemp <- rep(snowtemp, length(SWE))
    } else if (length(snowtemp) != length(SWE)) {
      stop("snowtemp must be length 1 or equal to the number of snow layers")
    }
    # layer age
    snowage <- rep(snowage, length(SWE))
  } else {
    SWE <- numeric(0)
    dz <- numeric(0)
    rho <- numeric(0)
    snowtemp <- numeric(0)
    snowage <- numeric(0)
  }
  list(
    SWE = SWE,             # mm water equivalent in each layer
    depth = dz,            # thickness of each layer (m)
    rho = rho,             # density of each layer (kg m-3)
    oldtemp = snowtemp,    # temperature of each layer in previous timestep (deg C)
    temp = snowtemp,       # temperature of each layer (deg C)
    age = snowage,         # age of each layer (hours)
    melt = 0,              # mm over timestep
    snowfall = 0,          # mm over timestep
    rainfall = 0,          # mm over timestep
    throughfall = 0,       # liquid water reaching ground (mm)
    surfaceWater = 0       # temporary liquid held in snowpack (mm)
  )
}
#' @title Snow albedo as a function of age
#' @description Computes broadband snow albedo as a function of snow age using
#' an empirical logarithmic decay relationship.
#' @param snowage Snow age (hours).
#' @return Numeric value (or vector) of snow albedo (dimensionless, 0 to 1).
#' @export
SnowAlbedo <- function(snowage) {
  alb <- (-9.8740 * log(snowage / 24) + 78.3434) / 100
  alb[alb < 0.1] <- 0.1
  alb[alb > 0.97] <- 0.97
  return(alb)
}

#' @title Aerodynamic resistance above snow
#' @description Computes aerodynamic resistance to heat transfer between a snow
#' surface and the reference height, accounting for partial burial of vegetation
#' by snow where relevant.
#' @param snowdepth Snow depth (m).
#' @param hgt Vegetation height (m).
#' @param pai Plant area index (m^2 m^-2).
#' @param uref Wind speed at the reference height (m s^-1).
#' @param psim Stability correction term for momentum (dimensionless).
#' @param psih Stability correction term for heat (dimensionless).
#' @param zref Reference height above the ground surface (m).
#'
#' @return Aerodynamic resistance to heat transfer (\eqn{s m^{-1}}).
#'
#' @details
#' If snow partially buries the vegetation canopy, vegetation height and plant
#' area index are reduced in proportion to the exposed canopy above the snow
#' surface. Zero-plane displacement height and roughness length for momentum are
#' then recalculated for the reduced canopy using [zeroplanedis()] and
#' [roughlength()].
#' @export
SnowResistance <- function(snowdepth, hgt, pai, uref, psim, psih, zref) {
  if (snowdepth >= zref) {
    stop("zref must be greater than snow depth")
  } else if (hgt > snowdepth) {
    # Adjust veg parameters for presence of snow
    pai <- pai * (hgt - snowdepth) / hgt
    hgt <- hgt - snowdepth
    d <- zeroplanedis(hgt, pai)
    zm <- roughlength(hgt, pai, d, psih)
    zh <- 0.2 * zm # roughness length for heat
    uf <- (0.41 * uref)  / (log((zref - d) / zm) + psim)
    rHa <- (log((zref - d) / zh) + psih) / (0.41 * uf)
  } else {
    # Bare snow surface
    zm <- 0.002 * exp(psih)
    zh <- 0.2 * zm # roughness length for heat
    uf <- (0.41 * uref)  / (log((zref - snowdepth) / zm) + psim)
    rHa <- (log((zref - snowdepth) / zh) + psih) / (0.41 * uf)
  }
  return(rHa)
}

#' @title One-step multilayer snow mass balance
#' @description Updates a multilayer snowpack over one timestep using climate forcing.
#' @param snowmod snow model state as returned by [InitSnowModel()]
#' @param climdata list of climate inputs
#' @param dT timestep (s)
#' @param snow_temp_thresh air temperature at or below which all precipitation
#' is treated as snow (degrees C)
#' @param rain_temp_thresh air temperature at or above which all precipitation
#' is treated as rain (degrees C)
#' @param liquid_water_frac maximum liquid water retained in the snowpack as a
#' fraction of remaining solid snow water equivalent
#' @param snowenv snow environment passed to [SnowDensParams()] and [SnowDensFun()]
#' @return Updated snow model state, including layer SWE, depth, density,
#' temperature, age, melt, snowfall, rainfall, throughfall, surface water,
#' sublimation, deposition and layer melt.
#' @export
SnowMassBalance <- function(snowmod, climdata,
                            dT = 3600,
                            snow_temp_thresh = 0,
                            rain_temp_thresh = 2,
                            liquid_water_frac = 0.05,
                            snowenv = "Maritime") {
  # --- constants ---
  Lf   <- 3.34e5     # latent heat of fusion (J kg-1)
  cice <- 2100       # heat capacity of ice (J kg-1 K-1)
  Rd   <- 287.05     # gas constant for dry air
  eps  <- 0.622      # molecular weight ratio
  # --- forcing ---
  precip <- climdata$precip
  Tair   <- climdata$Tair
  RH     <- climdata$relhum
  rHa    <- climdata$rHa
  pressure <- climdata$pk * 1000
  # --- fresh snow density ---
  sdp <- SnowDensParams(snowenv)
  rho_new_snow <- SnowDensFun(sdp, 0, 0)
  # --- state variables ---
  SWE     <- snowmod$SWE
  depth   <- snowmod$depth
  rho     <- snowmod$rho
  oldtemp <- snowmod$oldtemp
  temp    <- snowmod$temp
  age     <- snowmod$age
  liquid0 <- snowmod$surfaceWater
  # use midpoint temperature for physics
  temp_av <- 0.5 * (temp + oldtemp)
  # ensure empty vectors if no snow
  if (length(SWE) == 0) {
    SWE <- depth <- rho <- oldtemp <- temp <- temp_av <- age <- numeric(0)
  }
  # --- helpers ---
  SpecHumSat <- function(Tc, P) {
    es <- satvap(Tc)
    eps * es / (P - (1 - eps) * es)
  }
  # remove mass progressively from top layers
  remove_from_top <- function(SWEvec, amount) {
    removed <- 0
    if (amount <= 0 || length(SWEvec) == 0) {
      return(list(SWE = SWEvec, removed = removed))
    }
    for (i in seq_along(SWEvec)) {
      if (removed >= amount) break
      take <- min(SWEvec[i], amount - removed)
      SWEvec[i] <- SWEvec[i] - take
      removed <- removed + take
    }
    list(SWE = SWEvec, removed = removed)
  }
  # remove empty layers and compress state vectors
  compress_snowpack <- function(SWE, depth, rho, oldtemp, temp, age) {
    keep <- SWE > 1e-10
    list(
      SWE = SWE[keep],
      depth = depth[keep],
      rho = rho[keep],
      oldtemp = oldtemp[keep],
      temp = temp[keep],
      age = age[keep]
    )
  }
  # --- precipitation partitioning ---
  if (Tair <= snow_temp_thresh) {
    fsnow <- 1
  } else if (Tair >= rain_temp_thresh) {
    fsnow <- 0
  } else {
    fsnow <- (rain_temp_thresh - Tair) /
      (rain_temp_thresh - snow_temp_thresh)
  }
  snowfall <- precip * fsnow
  rainfall <- precip * (1 - fsnow)
  # --- add fresh snow as new surface layer ---
  if (snowfall > 0) {
    Tnew <- min(Tair, 0)
    SWE   <- c(snowfall, SWE)
    depth <- c(snowfall / rho_new_snow, depth)
    rho   <- c(rho_new_snow, rho)
    oldtemp <- c(Tnew, oldtemp)
    temp    <- c(Tnew, temp)
    temp_av <- c(Tnew, temp_av)
    age     <- c(0, age)
  }
  # --- rain-induced melt (wash-off) ---
  rainmelt <- 0
  if (rainfall > 0 && Tair > 0 && sum(SWE) > 0) {
    rainmelt <- 0.0125 * Tair * rainfall
    out <- remove_from_top(SWE, rainmelt)
    SWE <- out$SWE
    rainmelt <- out$removed
  }
  # --- sublimation / deposition ---
  sublimation <- 0
  deposition  <- 0
  if (sum(SWE) > 0 && is.finite(rHa) && rHa > 0) {
    RH_use <- if (RH > 1) RH / 100 else RH
    RH_use <- max(0, min(RH_use, 1))
    Tsurf <- min(temp_av[1], 0)
    qair  <- RH_use * SpecHumSat(Tair, pressure)
    qsurf <- SpecHumSat(Tsurf, pressure)
    rhoair <- pressure / (Rd * (Tair + 273.15))
    sub_flux <- rhoair * (qsurf - qair) / rHa
    sub_mass <- sub_flux * dT
    if (sub_mass > 0) {
      # sublimation removes mass
      out <- remove_from_top(SWE, sub_mass)
      SWE <- out$SWE
      sublimation <- out$removed
    } else if (sub_mass < 0) {
      # deposition adds mass to surface
      deposition <- -sub_mass
      Tnew <- min(Tair, 0)
      if (length(SWE) == 0) {
        SWE   <- deposition
        rho   <- rho_new_snow
        depth <- deposition / rho_new_snow
        oldtemp <- Tnew
        temp    <- Tnew
        temp_av <- Tnew
        age     <- 0
      } else {
        SWE[1] <- SWE[1] + deposition
      }
    }
  }
  # --- internal melt (layers above 0°C) ---
  layermelt <- rep(0, length(SWE))
  if (length(SWE) > 0) {
    for (i in seq_along(SWE)) {
      if (SWE[i] > 0 && temp_av[i] > 0) {
        melt_i <- SWE[i] * cice * temp_av[i] / Lf
        melt_i <- min(melt_i, SWE[i])
        layermelt[i] <- melt_i
        SWE[i] <- SWE[i] - melt_i
        oldtemp[i] <- 0
        temp[i] <- 0
        temp_av[i] <- 0
      }
    }
  }
  melt <- sum(layermelt)
  # --- liquid water storage and drainage ---
  liquid_in <- liquid0 + rainfall + rainmelt + melt
  max_liquid_storage <- liquid_water_frac * sum(SWE)
  retained_liquid <- min(liquid_in, max_liquid_storage)
  throughfall <- liquid_in - retained_liquid
  # --- ageing ---
  if (length(age) > 0) {
    age <- age + dT / 3600
  }
  # --- remove empty layers ---
  packed <- compress_snowpack(SWE, depth, rho, oldtemp, temp, age)
  SWE     <- packed$SWE
  depth   <- packed$depth
  rho     <- packed$rho
  oldtemp <- packed$oldtemp
  temp    <- packed$temp
  age     <- packed$age
  # --- update density profile ---
  if (length(SWE) > 0) {
    sdp <- SnowDensParams(snowenv)
    depth <- SWE / rho
    cdz <- cumsum(depth)
    rho <- SnowDensFun(sdp, cdz, age)
    depth <- SWE / rho
  } else {
    rho <- depth <- oldtemp <- temp <- age <- numeric(0)
    retained_liquid <- 0
  }
  # --- update state ---
  snowmod$SWE <- SWE
  snowmod$depth <- depth
  snowmod$rho <- rho
  snowmod$oldtemp <- oldtemp
  snowmod$temp <- temp
  snowmod$age <- age
  # diagnostics
  snowmod$melt <- melt
  snowmod$snowfall <- snowfall
  snowmod$rainfall <- rainfall
  snowmod$throughfall <- throughfall
  snowmod$surfaceWater <- retained_liquid
  snowmod$rainmelt <- rainmelt
  snowmod$sublimation <- sublimation
  snowmod$deposition <- deposition
  snowmod$layermelt <- layermelt
  snowmod
}

#' @title Regularises snow layer thickness
#' @description Merges overly thin layers and splits overly thick layers to maintain
#' a numerically sensible snow discretisation for the heat model.
#' @param snowmod snow model state as returned by [InitSnowModel()] or updated by model steps
#' @param snowenv snow environment passed to [SnowDensParams()]
#' @return updated \code{snowmod} with adjusted layer structure
#' @export
RegulariseSnowLayers <- function(snowmod, snowenv = "Maritime") {
  # Set min and max thickness
  min_thickness = 0.005
  max_thickness = 0.05
  split_factor = 2
  # Extract form snow model
  SWE <- snowmod$SWE
  depth <- snowmod$depth
  rho <- snowmod$rho
  temp <- snowmod$temp
  oldtemp <- snowmod$oldtemp
  age <- snowmod$age
  if (length(SWE) == 0) return(snowmod)
  # remove empty layers first
  keep <- SWE > 1e-10 & depth > 1e-10
  SWE <- SWE[keep]
  depth <- depth[keep]
  rho <- rho[keep]
  temp <- temp[keep]
  oldtemp <- oldtemp[keep]
  age <- age[keep]
  if (length(SWE) == 0) {
    snowmod$SWE <- numeric(0)
    snowmod$depth <- numeric(0)
    snowmod$rho <- numeric(0)
    snowmod$temp <- numeric(0)
    snowmod$oldtemp <- numeric(0)
    snowmod$age <- numeric(0)
    return(snowmod)
  }
  # --- merge overly thin adjacent layers ---
  i <- 1
  while (i <= length(SWE)) {
    if (depth[i] < min_thickness && length(SWE) > 1) {
      # merge with neighbour giving smaller combined thickness
      if (i == 1) {
        j <- 2
      } else if (i == length(SWE)) {
        j <- i - 1
      } else {
        j <- if ((depth[i - 1] + depth[i]) <= (depth[i] + depth[i + 1])) i - 1 else i + 1
      }
      ii <- sort(c(i, j))
      SWE_new <- sum(SWE[ii])
      depth_new <- sum(depth[ii])
      rho_new <- SWE_new / depth_new
      temp_new <- weighted.mean(temp[ii], SWE[ii])
      oldtemp_new <- weighted.mean(oldtemp[ii], SWE[ii])
      # geometric mean for age, weighted by SWE
      age_pos <- pmax(age[ii], 1e-6)
      age_new <- exp(sum(SWE[ii] * log(age_pos)) / sum(SWE[ii]))
      k <- ii[1]
      SWE[k] <- SWE_new
      depth[k] <- depth_new
      rho[k] <- rho_new
      temp[k] <- temp_new
      oldtemp[k] <- oldtemp_new
      age[k] <- age_new
      drop <- ii[2]
      SWE <- SWE[-drop]
      depth <- depth[-drop]
      rho <- rho[-drop]
      temp <- temp[-drop]
      oldtemp <- oldtemp[-drop]
      age <- age[-drop]
      i <- max(1, k - 1)
    } else {
      i <- i + 1
    }
  }
  # --- split overly thick layers ---
  i <- 1
  while (i <= length(SWE)) {
    if (depth[i] > max_thickness) {
      nsplit <- max(2, ceiling(depth[i] / max_thickness))
      nsplit <- max(nsplit, split_factor)
      SWE_split <- rep(SWE[i] / nsplit, nsplit)
      depth_split <- rep(depth[i] / nsplit, nsplit)
      rho_split <- rep(SWE[i] / depth[i], nsplit)
      temp_split <- rep(temp[i], nsplit)
      oldtemp_split <- rep(oldtemp[i], nsplit)
      age_split <- rep(age[i], nsplit)
      SWE <- append(SWE[-i], SWE_split, after = i - 1)
      depth <- append(depth[-i], depth_split, after = i - 1)
      rho <- append(rho[-i], rho_split, after = i - 1)
      temp <- append(temp[-i], temp_split, after = i - 1)
      oldtemp <- append(oldtemp[-i], oldtemp_split, after = i - 1)
      age <- append(age[-i], age_split, after = i - 1)
      i <- i + nsplit
    } else {
      i <- i + 1
    }
  }
  # --- recompute density profile consistently from depth and age ---
  if (length(SWE) > 0) {
    sdp <- SnowDensParams(snowenv)
    depth <- SWE / rho
    cdz <- cumsum(depth)
    rho <- SnowDensFun(sdp, cdz, age)
    depth <- SWE / rho
  }
  snowmod$SWE <- SWE
  snowmod$depth <- depth
  snowmod$rho <- rho
  snowmod$temp <- temp
  snowmod$oldtemp <- oldtemp
  snowmod$age <- age
  snowmod
}


#' @title Snow thermal conductivity
#' @description Calculates snow thermal conductivity from bulk snow density
#' using an exponential relationship following Djachkova (after Anderson, 2006).
#' @param rho Snow density (kg m-3).
#' @return Snow thermal conductivity (W m-1 K-1).
#' @export
SnowThermalConductivity <- function(rho) {
  rhog <- rho / 1000
  k <- 0.0442 * exp(5.181 * rhog) # convert kg m-3 -> g cm-3
  return(k)
}

#' @title Snow volumetric heat capacity
#' @description Calculates volumetric heat capacity of snow as a mixture of ice,
#' air, and optionally liquid water.
#' @param rho Snow density (kg m-3).
#' @param Tc Snow temperature (deg C).
#' @param pk Air pressure (kPa).
#' @return Snow volumetric heat capacity (J m-3 K-1).
#' @export
SnowHeatCapacity <- function(rho, Tc, pk = 101.3) {
  Vw = 0 # set water content of snow to zero
  rhoi <- 917          # density of ice (kg m-3)
  # volumetric fractions
  Vi <- rho / rhoi
  Va <- pmax(1 - Vi - Vw, 0)
  # volumetric heat capacities (J m-3 K-1)
  Ci <- 1.93e6 + 0.0067e6 * Tc   # ice (temp dependent below 0)
  Ci[Tc >= 0] <- 2.1e6
  Ca <- cpair(Tc) * rhohair(Tc, pk)   # air
  Cw <- 4.18e6                        # water
  # mixture
  Cd <- Vi * Ci + Va * Ca + Vw * Cw
  return(Cd)
}

#' @title Calculates snow surface energy balance
#' @description Calculates the snow surface energy balance (W m-2)
#' @param Rabs absorbed radiation at the surface (W m-2)
#' @param Tref reference air temperature (degrees C)
#' @param Tsurface snow surface temperature (degrees C)
#' @param pk air pressure (kPa)
#' @param relhum air relative humidity (Percentage)
#' @param rHa aerodynamic resistance to heat from snow surface to reference height
#' (s m-1), e.g. as returned by [SnowResistance()]
#' @return snow surface energy balance (W m-2)
#' @export
snowsurfaceEB <- function(Rabs, Tref, Tsurface, pk, relhum, rHa) {
  # net radiation
  sb <- 5.67e-8
  Rnet <- Rabs - 0.99 * sb * (Tsurface + 273.15)^4
  # sensible heat
  cp <- cpair(Tref)
  ph <- rhohair(Tref, pk)
  H <- (ph * cp / rHa) * (Tsurface - Tref)
  # latent heat
  es <- satvap(Tsurface)
  ea <- satvap(Tref) * (relhum / 100)
  la <- latvap(Tsurface)
  L <- (la * ph / (rHa * pk)) * (es - ea)
  Rnet - H - L
}
#' Internal functions used by snow-soil heat solver
.appendSnowSoil <- function(snowmod, soilheatmod, soilp, atmPressure) {
  ns <- length(snowmod$depth)
  nz_soil <- length(soilheatmod$dz)
  # prognostic temperatures from previous time step
  T_old <- c(snowmod$oldtemp, soilheatmod$oldTe)
  # layer thicknesses / control volumes
  dz <- c(snowmod$depth, soilheatmod$dz)
  # thermal conductivity
  k_snow <- if (ns > 0) {
    SnowThermalConductivity(snowmod$rho)
  } else {
    numeric(0)
  }
  k_soil <- thermalConductivity(
    soilp$Vq, soilp$Vm, soilp$Vo,
    soilheatmod$wc, soilp$Mc,
    soilheatmod$Te, atmPressure
  )
  k <- c(k_snow, k_soil)
  # volumetric heat capacity × volume
  C_snow <- if (ns > 0) {
    SnowHeatCapacity(snowmod$rho, snowmod$temp, atmPressure) * snowmod$depth
  } else {
    numeric(0)
  }
  C_soil <- heatCapacity(
    soilp$Vq, soilp$Vm, soilp$Vo,
    soilheatmod$wc, soilheatmod$Te, atmPressure
  ) * soilheatmod$vol
  C <- c(C_snow, C_soil)
  # interface conductance between adjacent node centres
  Gint <- 1 / (0.5 * dz[-length(dz)] / k[-length(k)] +
                 0.5 * dz[-1] / k[-1])
  list(
    ns = ns,
    nz = length(dz),
    T_old = T_old,
    dz = dz,
    C = C,
    k = k,
    Gint = Gint
  )
}
.splitSnowSoil <- function(T_new, snowmod, soilheatmod, Gflux_total) {
  ns <- length(snowmod$depth)
  if (ns > 0) {
    snowmod$temp <- T_new[seq_len(ns)]
  } else {
    snowmod$temp <- numeric(0)
  }
  soilheatmod$Te <- T_new[ns + seq_along(soilheatmod$Te)]
  soilheatmod$Gflux <- Gflux_total
  soilheatmod$iters <- NA_integer_
  snowmod$Gflux <- Gflux_total
  list(snowmod = snowmod, soilheatmod = soilheatmod)
}

#' @title Updates snow and soil temperatures using the coupled snow-soil heat model
#' @description Solves heat conduction through the combined snow and soil profile over one
#' time step, updating layer temperatures by combining surface energy balance from
#' [snowsurfaceEB()] or [soilsurfaceEB()], thermal conductivity from
#' [SnowThermalConductivity()] and [thermalConductivity()], volumetric heat capacity from
#' [SnowHeatCapacity()] and [heatCapacity()], and tridiagonal solution with
#' [ThomasBoundaryCondition()].
#' @param snowmod list of variables provided as input and updated by the snow
#' model as returned by [InitSnowModel()] on first step or by previous model run
#' @param soilheatmod list of variables provided as input and updated by soil
#' heat model as returned by [InitSoilheatmod()] on first step or by previous model run
#' @param soilp list of soil parameters as returned by [createsoilplist()]
#' @param Rabs absorbed radiation at the surface (W m-2)
#' @param Tref reference air temperature (degrees C)
#' @param relhum air relative humidity (Percentage)
#' @param atmPressure air pressure (kPa)
#' @param rHa aerodynamic resistance to heat from the surface to reference height
#' as returned by [SnowResistance()] or [groundresistance()] (s m-1)
#' @param dT time step (s)
#' @param Fact weighting factor for the implicit solution scheme
#' @param maxNrIterations maximum number of nonlinear iterations
#' @param tolerance convergence criterion for maximum change in layer temperature
#' between successive iterations (degrees C)
#' @param relax under-relaxation factor applied to temperature updates
#' @return a list comprising
#' \describe{
#'   \item{snowmod}{updated snow model state}
#'   \item{soilheatmod}{updated soil heat model state}
#' }
#' @export
SnowSoilHeatModel <- function(snowmod, soilheatmod, soilp,
                              Rabs, Tref, relhum, atmPressure, rHa, dT = 3600, Fact = 0.5,
                              maxNrIterations = 100, tolerance = 1e-2, relax = 0.3) {

  # --- combine snow and soil profiles ---
  prof <- .appendSnowSoil(snowmod, soilheatmod, soilp, atmPressure)
  ns <- prof$ns
  nz <- prof$nz
  gg <- 1 - Fact
  T_old <- prof$T_old
  T_new <- T_old
  dz <- prof$dz
  # fixed lower boundary
  boundaryT <- T_old[nz]
  nrIterations <- 0
  maxdT <- Inf
  prev_maxdT <- Inf
  # surface flux damping bookkeeping
  max_qsurface <- 0
  qsurf1 <- NA_real_
  qsurf2 <- NA_real_
  # adaptive relaxation parameters
  relax_now <- relax
  relax_min <- 0.05
  relax_max <- 0.7
  while (maxdT > tolerance && nrIterations < maxNrIterations) {
    # --- surface energy balance ---
    Tav <- 0.5 * (T_old[1] + T_new[1])

    qsurface <- if (ns > 0) {
      snowsurfaceEB(Rabs, Tref, Tav, atmPressure, relhum, rHa)
    } else {
      soilsurfaceEB(soilp, Rabs, Tref, Tav, atmPressure, relhum, rHa,
                    soilheatmod$wc[1])
    }
    # damp nonlinear surface flux in early iterations
    if (nrIterations < 10) {
      if (nrIterations < 2 && abs(qsurface) > abs(max_qsurface)) {
        max_qsurface <- abs(qsurface)
      }
      if (nrIterations > 1 && abs(qsurface) > abs(max_qsurface)) {
        qsurface <- sign(qsurface) * abs(max_qsurface)
      }
      if (nrIterations == 8) qsurf1 <- qsurface
      if (nrIterations == 9) qsurf2 <- qsurface
    } else if (nrIterations == 10 && is.finite(qsurf1) && is.finite(qsurf2)) {
      qsurface <- 0.5 * (qsurf1 + qsurf2)
    }
    # --- update thermal properties ---
    if (ns > 0) {
      snow_T <- T_new[seq_len(ns)]
      k_snow <- SnowThermalConductivity(snowmod$rho)
      C_snow <- SnowHeatCapacity(snowmod$rho, snow_T, atmPressure) * snowmod$depth
    } else {
      k_snow <- numeric(0)
      C_snow <- numeric(0)
    }
    soil_T <- T_new[ns + seq_along(soilheatmod$Te)]
    k_soil <- thermalConductivity(
      soilp$Vq, soilp$Vm, soilp$Vo,
      soilheatmod$wc, soilp$Mc, soil_T, atmPressure
    )
    C_soil <- heatCapacity(
      soilp$Vq, soilp$Vm, soilp$Vo,
      soilheatmod$wc, soil_T, atmPressure
    ) * soilheatmod$vol

    k <- c(k_snow, k_soil)
    C <- c(C_snow, C_soil)
    # --- interface conductance (resistance form) ---
    Gint <- 1 / (0.5 * dz[-length(dz)] / k[-length(k)] +
                   0.5 * dz[-1] / k[-1])
    # --- assemble tridiagonal system ---
    aa <- bb <- cc <- dd <- numeric(nz)
    # top node
    bb[1] <- C[1] / dT + Fact * Gint[1]
    cc[1] <- -Fact * Gint[1]
    dd[1] <- C[1] / dT * T_old[1] + qsurface +
      gg * Gint[1] * (T_old[2] - T_old[1])
    # internal nodes
    if (nz > 2) {
      idx <- 2:(nz - 1)
      aa[idx] <- -Fact * Gint[idx - 1]
      bb[idx] <- C[idx] / dT + Fact * (Gint[idx - 1] + Gint[idx])
      cc[idx] <- -Fact * Gint[idx]
      dd[idx] <- C[idx] / dT * T_old[idx] +
        gg * (Gint[idx - 1] * (T_old[idx - 1] - T_old[idx]) +
                Gint[idx] * (T_old[idx + 1] - T_old[idx]))
    }
    # bottom boundary (Dirichlet)
    aa[nz] <- 0
    bb[nz] <- 1
    cc[nz] <- 0
    dd[nz] <- boundaryT
    # --- solve system ---
    T_prev <- T_new
    T_raw <- ThomasBoundaryCondition(aa, bb, cc, dd, T_new, 1, nz - 1)$x
    # --- adaptive under-relaxation with backtracking ---
    repeat {
      # trial update
      T_trial <- (1 - relax_now) * T_prev + relax_now * T_raw
      # enforce monotonicity constraint
      if (nz > 2) {
        idx <- 2:(nz - 1)
        lo <- pmin(T_trial[idx - 1], T_trial[idx + 1])
        hi <- pmax(T_trial[idx - 1], T_trial[idx + 1])
        T_trial[idx] <- pmax(pmin(T_trial[idx], hi), lo)
      }
      trial_maxdT <- max(abs(T_trial - T_prev))
      # accept step if improvement, otherwise reduce relaxation
      if (trial_maxdT <= prev_maxdT || relax_now <= relax_min) break
      relax_now <- max(relax_min, 0.5 * relax_now)
    }
    T_new <- T_trial
    maxdT <- trial_maxdT
    # if converging well, cautiously increase relaxation
    if (maxdT < prev_maxdT) {
      relax_now <- min(relax_max, relax_now * 1.1)
    }
    prev_maxdT <- maxdT
    nrIterations <- nrIterations + 1
  }
  # --- total heat storage (snow + soil, excluding fixed boundary) ---
  Gflux_total <- sum(C[1:(nz - 1)] *
                       (T_new[1:(nz - 1)] - T_old[1:(nz - 1)]) / dT)
  # --- split back into snow and soil ---
  out <- .splitSnowSoil(T_new, snowmod, soilheatmod, Gflux_total)
  out$soilheatmod$iters <- nrIterations
  return(out)
}

#' @title Calculates slope correction for direct radiation using Ross G
#' @description Returns the fraction of direct beam radiation incident on an inclined
#' surface, analogous to [solarindex()], but with reduced angular dependence
#' controlled by \code{y}
#' @param solp a list of solar zenith and azimuths as returned by [sunposition()]
#' @param y parameter controlling angular response, where \code{y = 0} gives a
#' planar surface and \code{y = 1} gives isotropic response
#' @param slope surface slope (degrees)
#' @param aspect surface aspect (degrees)
#' @return fraction of direct beam radiation incident on the inclined surface
#' @export
rossindex <- function(solp, y = 0.7, slope = 0, aspect = 0) {
  zen <- solp$zen * pi / 180
  azi <- solp$azi * pi / 180
  slope <- slope * pi / 180
  aspect <- aspect * pi / 180
  # cosine of incidence angle on the slope
  mu <- cos(zen) * cos(slope) +
    sin(zen) * sin(slope) * cos(azi - aspect)
  mu[mu < 0] <- 0
  mu[solp$zen > 90] <- 0
  if (length(y) == 1) y <- rep(y, length(mu))
  if (length(y) != length(mu)) {
    stop("y must have length 1 or the same length as solp$zen")
  }
  y <- pmax(0, pmin(1, y))
  G <- numeric(length(mu))
  # y = 0 -> x = Inf -> Lambert cosine law
  i0 <- (y == 0)
  G[i0] <- mu[i0]
  # y = 1 -> x = 1 -> isotropic absorption
  i1 <- (y == 1)
  G[i1] <- 0.5
  # general case
  ig <- !(i0 | i1)
  if (any(ig)) {
    x <- 1 / y[ig]
    ii <- acos(pmin(1, pmax(0, mu[ig])))
    k <- sqrt(x^2 + tan(ii)^2) /
      (x + 1.774 * (x + 1.182)^(-0.733))
    G[ig] <- k * cos(ii)
  }
  G[solp$zen > 90] <- 0
  G
}

#' Calculate canopy interception of snowfall
#'
#' Estimates the amount of snowfall intercepted by vegetation canopy over a single
#' time step, following the formulation in the supplied C++ implementation.
#'
#' The function computes a mean within-canopy wind speed from canopy structure and
#' friction velocity, uses this to estimate the fall angle of snow, and then
#' calculates effective canopy cover and canopy snow storage capacity. Interception
#' is limited by both incoming snowfall and the remaining canopy storage capacity.
#'
#' @param hgt Canopy height (m).
#' @param pai Plant area index (dimensionless).
#' @param uf Friction velocity (m s^-1).
#' @param prec Snowfall during the time step (mm snow water equivalent).
#' @param tc Air temperature (degrees C).
#' @param Li Snow load on the canopy from the previous time step
#'   (mm snow water equivalent).
#' @param Sh Snow load per unit branch area coefficient (kg m^-3).
#'   Default is `6.2`.
#'
#' @return Numeric value giving canopy-intercepted snowfall for the time step
#'   (mm snow water equivalent).
#'
#' @details
#' Very small values of `hgt` and `pai` are bounded below at 0.001 to avoid
#' numerical problems. Fresh snow density is estimated as an exponential function
#' of air temperature, and the maximum canopy snow load is scaled by plant area
#' index. The returned interception is capped so that it cannot exceed the
#' snowfall input `prec`.
#'
#' @export
CanopySnowInterception <- function(hgt, pai, uf, prec, tc, Li, Sh = 6.2) {
  # Avoid zero values
  hgt <- pmax(hgt, 0.001)
  pai <- pmax(pai, 0.001)
  # Calculate mean canopy wind
  Be <- sqrt(0.003 + (0.2 * pai) / 2.0)
  uh <- uf / Be
  a <- pai / hgt
  Lc <- (0.25 * a)^(-1.0)
  Lm <- 2.0 * Be^3.0 * Lc
  k1 <- Be / Lm
  uzm <- (uh / (hgt * k1)) * (1 - exp(-k1 * hgt))
  uzm <- pmax(uzm, uf)
  # Calculate snow interception
  rhos <- 67.92 + 51.25 * exp(tc / 2.59)      # fresh snow density (kg/m^3)
  S <- Sh * (0.26 + 46 / rhos)                # maximum snow load per unit branch area (kg/m^3)
  Lstr <- S * pai                             # maximum canopy snow load (mm SWE)
  Z <- atan(uzm / 0.8)                        # zenith angle of snow fall direction (0.8 m/s = terminal velocity of snow flake)
  kc <- 1.0 / (2.0 * cos(Z))                  # extinction coefficient (from Campbell assuming spherical distribution)
  Cp <- 1.0 - exp(-kc * pai)                  # effective canopy cover perpendicular to direction of snow flake
  k2 <- Cp / Lstr                             # dimensionless proportionality factor
  I1 <- (Lstr - Li) * (1.0 - exp(-k2 * prec)) # intercepted snow load at the start of unloading (mm SWE)
  cis <- I1 * 0.678                           # canopy snow interception (mm SWE)
  cis <- pmin(cis, prec)
  return(cis)
}
