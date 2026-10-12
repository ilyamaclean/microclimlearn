#' @title Snow density parameter sets
#' @description Returns empirical parameter values for a snow density evolution
#' function, based on different climatic snow regimes.
#'
#' These parameterisations are intended for use with [SnowDensFun()] to estimate
#' bulk snow density as a function of snow depth and age.
#'
#' @param snowenv Character string specifying the snow environment. Options are:
#' `Alpine`, `Maritime`, `Prairie`, `Tundra` or `Taiga`.
#'
#' @return A numeric vector of length 4 containing density parameters:
#' \describe{
#'   \item{densfun[1]}{Upper asymptotic snow density (g cm^-3)}
#'   \item{densfun[2]}{Initial / minimum snow density (g cm^-3)}
#'   \item{densfun[3]}{Depth-dependent densification coefficient (cm^-1)}
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
  # Snow is densest where most snow lies above it and where it has had longest
  # to settle: density rises with depth (cm) and age (days) towards an upper limit
  snowdensity <- ((densfun[1] - densfun[2]) *
                    (1 - exp(-densfun[3] * snowdepth * 100 -
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
    # Layers are ordered from the snow surface downwards. They are 1 cm thick at
    # the surface, where temperature changes fastest, and thicken to 3 cm below
    dz <- pmin(0.01 * 1.6^(seq_len(ceiling(snowdepth / 0.01)) - 1), 0.03)
    dz <- dz[cumsum(dz) <= snowdepth]
    # The snow left over forms the lowest layer, or is added to the layer above
    # it when there is too little of it to make a layer of its own
    rest <- snowdepth - sum(dz)
    if (length(dz) > 0 && rest < 0.005) {
      dz[length(dz)] <- dz[length(dz)] + rest
    } else {
      dz <- c(dz, rest)
    }
    # Snow is denser the deeper it lies beneath the surface
    cdz <- cumsum(dz)
    sdp <- SnowDensParams(snowenv)
    rho <- SnowDensFun(sdp, cdz, snowage)
    # Water held by each layer as snow (kg m-2, equivalent to mm of water)
    SWE <- rho * dz
    # Every layer starts at the same temperature unless one is given for each
    if (length(snowtemp) == 1) {
      snowtemp <- rep(snowtemp, length(SWE))
    } else if (length(snowtemp) != length(SWE)) {
      stop("snowtemp must be length 1 or equal to the number of snow layers")
    }
    # All of the snow is taken to have fallen at the same time
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
    oldtemp = snowtemp,    # temperature of each layer at the start of the time step (deg C)
    temp = snowtemp,       # temperature of each layer at the end of the time step (deg C)
    age = snowage,         # age of each layer (hours)
    melt = 0,              # snow melted over the time step (mm)
    snowfall = 0,          # snow added over the time step (mm)
    rainfall = 0,          # rain added over the time step (mm)
    throughfall = 0,       # liquid water leaving the base of the pack (mm)
    surfaceWater = 0       # liquid water held within the pack (mm)
  )
}
#' @title Snow albedo as a function of age
#' @description Computes broadband snow albedo as a function of snow age using
#' an empirical logarithmic decay relationship.
#' @param snowage Snow age (hours).
#' @return Numeric value (or vector) of snow albedo (dimensionless, 0 to 1).
#' @export
SnowAlbedo <- function(snowage) {
  # Snow darkens as it ages: its grains grow and its surface collects dirt.
  # Albedo falls with the logarithm of the days since the snow fell
  alb <- (-9.8740 * log(snowage / 24) + 78.3434) / 100
  # Bounded by the albedo of fresh snow and that of very old, dirty snow
  alb[alb < 0.1] <- 0.1
  alb[alb > 0.85] <- 0.85
  return(alb)
}

#' @title Aerodynamic resistance above snow
#' @description Calculates the resistance to heat transfer between a snow
#' surface and the reference height, and the friction velocity of the wind above
#' the snow. Vegetation is taken to be buried by the snow up to the snow depth,
#' so that only the part standing above the snow surface affects the wind.
#' @param snowdepth Snow depth (m).
#' @param hgt Vegetation height above the ground (m). Zero if there is none.
#' @param pai Plant area index of the vegetation (m^2 m^-2).
#' @param uref Wind speed at the reference height (m s^-1).
#' @param zref Reference height above the ground surface (m). Must be above
#' both the snow and the vegetation.
#' @param LL Obukhov length as returned by [Obukhov()] (m). The default, `Inf`,
#' is neutral air.
#' @param MOST Logical. If `TRUE`, profiles above the vegetation are taken to
#' be logarithmic down to canopy top. If `FALSE`, a roughness sublayer is assumed
#' above it. Has no effect when the vegetation is buried.
#' @param z0 roughness length of the snow surface (m)
#' @return a list of the following:
#' \describe{
#'  \item{rHa}{resistance to heat transfer between the snow surface and `zref` (s m^-1)}
#'  \item{uf}{friction velocity (m s^-1)}
#' }
#' @details
#' Plant area is taken to be spread evenly with height, so the plant area
#' standing above the snow is in proportion to the height standing above it.
#' Where the snow buries the vegetation completely, the surface is treated as
#' bare.
#' @export
SnowResistance <- function(snowdepth, hgt, pai, uref, zref, LL = Inf, MOST = TRUE, z0 = 0.002) {
  if (snowdepth >= zref) stop("zref must be greater than snow depth")
  if (zref < hgt) stop("zref must be above the vegetation")
  # Only the vegetation standing above the snow slows the wind, and heights are
  # measured from the snow surface
  hexp <- max(hgt - snowdepth, 0)
  paiexp <- 0
  if (hgt > 0) paiexp <- pai * hexp / hgt
  zs <- zref - snowdepth
  # Vegetation that is wholly buried leaves a bare surface, with no canopy to
  # form a roughness sublayer
  if (hexp == 0) MOST <- TRUE
  # Friction velocity above the snow and what stands above it, and the
  # resistance of the air between the snow surface and the reference height
  uf <- windfric(uref, zs, hexp, paiexp, LL, MOST, z0)
  rHa <- .groundtoheight(zs, uf, hexp, paiexp, LL, MOST, z0)
  return(list(rHa = rHa, uf = uf))
}

#' @title One-step multilayer snow mass balance
#' @description Updates the water held by a multilayer snowpack over one time
#' step: snow that falls, snow that melts, snow lost or gained as vapour, and
#' liquid water held in the pack or draining from its base.
#' @param snowmod snow model state as returned by [InitSnowModel()], with layer
#' temperatures as updated by [SnowSoilHeatModel()]
#' @param climdata list of driving variables for the time step. Must contain:
#' \describe{
#'   \item{Tair}{air temperature (degrees C)}
#'   \item{relhum}{air relative humidity (Percentage)}
#'   \item{pk}{air pressure (kPa)}
#'   \item{rHa}{resistance to heat transfer between the snow surface and the
#'   reference height as returned by [SnowResistance()] (s m^-1). An infinite
#'   value gives no sublimation or deposition.}
#'   \item{precip}{precipitation over the time step (mm)}
#' }
#' @param dT timestep (s)
#' @param snow_temp_thresh air temperature at or below which all precipitation
#' is treated as snow (degrees C)
#' @param rain_temp_thresh air temperature at or above which all precipitation
#' is treated as rain (degrees C)
#' @param liquid_water_frac maximum liquid water retained in the snowpack as a
#' fraction of remaining solid snow water equivalent
#' @param snowenv snow environment passed to [SnowDensParams()] and [SnowDensFun()]
#' @details The heat model has no melting, so it can leave snow warmer than 0
#' degrees C. Here the heat a layer holds above 0 degrees C melts snow, and the
#' layer is returned to 0 degrees C. Heat left over once the whole pack has
#' melted is returned as `heattoground`. Liquid water passing down through snow
#' colder than 0 degrees C freezes, which warms the snow. `oldtemp` is left as
#' it is; to move on to the next time step, set `oldtemp` to `temp`.
#' @return Updated snow model state. As well as the state of each layer (`SWE`,
#' `depth`, `rho`, `temp`, `oldtemp`, `age`) it holds, for the time step:
#' \describe{
#'   \item{melt}{snow melted by heat held above 0 degrees C (mm)}
#'   \item{layermelt}{the same for each layer present at the start of the step (mm)}
#'   \item{rainmelt}{snow melted by the heat of rain (mm)}
#'   \item{refrozen}{liquid water frozen within the pack (mm)}
#'   \item{snowfall}{snow added (mm)}
#'   \item{rainfall}{rain added (mm)}
#'   \item{sublimation}{snow lost as vapour (mm)}
#'   \item{deposition}{vapour deposited as snow (mm)}
#'   \item{surfaceWater}{liquid water held within the pack at the end of the step (mm)}
#'   \item{throughfall}{liquid water leaving the base of the pack (mm)}
#'   \item{heattoground}{heat left over after all the snow has melted, which
#'   passes to the ground beneath (W m^-2)}
#' }
#' @export
SnowMassBalance <- function(snowmod, climdata,
                            dT = 3600,
                            snow_temp_thresh = 0,
                            rain_temp_thresh = 2,
                            liquid_water_frac = 0.05,
                            snowenv = "Maritime") {
  Lf <- 3.34e5 # latent heat of fusion (J kg-1)
  precip <- climdata$precip
  Tair <- climdata$Tair
  rHa <- climdata$rHa
  SWE <- snowmod$SWE
  depth <- snowmod$depth
  rho <- snowmod$rho
  oldtemp <- snowmod$oldtemp
  temp <- snowmod$temp
  age <- snowmod$age
  sdp <- SnowDensParams(snowenv)
  rho_new_snow <- SnowDensFun(sdp, 0, 0) # density of newly fallen snow (kg m-3)

  # Takes a given amount of water (mm) from the snowpack, starting with the
  # surface layer and working down
  remove_from_top <- function(SWEvec, amount) {
    removed <- 0
    for (i in seq_along(SWEvec)) {
      if (removed >= amount) break
      take <- min(SWEvec[i], amount - removed)
      SWEvec[i] <- SWEvec[i] - take
      removed <- removed + take
    }
    list(SWE = SWEvec, removed = removed)
  }

  # Melt. The heat a layer holds above 0 deg C melts it. Heat left over once a
  # layer has melted away passes to the layer beneath
  layermelt <- numeric(length(SWE))
  spare <- 0 # heat passed down from the layer above (J m-2)
  for (i in seq_along(SWE)) {
    # Heat needed to warm the layer by one degree (J m-2 K-1), and the heat it
    # holds relative to snow at 0 deg C
    C <- SnowHeatCapacity(rho[i], temp[i], climdata$pk) * SWE[i] / rho[i]
    heat <- C * temp[i] + spare
    spare <- 0
    if (heat > 0) {
      layermelt[i] <- min(heat / Lf, SWE[i])
      spare <- heat - layermelt[i] * Lf
      SWE[i] <- SWE[i] - layermelt[i]
      temp[i] <- 0
    } else {
      temp[i] <- heat / C
    }
  }
  melt <- sum(layermelt)
  # Heat that remains when all the snow has melted goes to the ground beneath
  heattoground <- spare / dT

  # Precipitation falls as snow in cold air and as rain in warm air, and as a
  # mixture of the two between
  fsnow <- (rain_temp_thresh - Tair) / (rain_temp_thresh - snow_temp_thresh)
  fsnow <- min(max(fsnow, 0), 1)
  snowfall <- precip * fsnow
  rainfall <- precip * (1 - fsnow)

  # New snow forms a layer at the surface, no warmer than 0 deg C
  if (snowfall > 0) {
    Tnew <- min(Tair, 0)
    SWE <- c(snowfall, SWE)
    depth <- c(snowfall / rho_new_snow, depth)
    rho <- c(rho_new_snow, rho)
    oldtemp <- c(Tnew, oldtemp)
    temp <- c(Tnew, temp)
    age <- c(0, age)
  }

  # Rain gives up heat as it cools to 0 deg C, which melts snow. The factor is
  # the specific heat of water over the latent heat of fusion
  rainmelt <- 0
  if (rainfall > 0 && Tair > 0) {
    out <- remove_from_top(SWE, 0.0125 * Tair * rainfall)
    SWE <- out$SWE
    rainmelt <- out$removed
  }

  # Sublimation and deposition. Vapour passes between the snow surface, where
  # the air is saturated, and the air above, down the gradient in vapour pressure
  sublimation <- 0
  deposition <- 0
  if (sum(SWE) > 0 && is.finite(rHa) && rHa > 0) {
    Tsurf <- min(0.5 * (temp[1] + oldtemp[1]), 0)
    es <- satvap(Tsurf)
    ea <- satvap(Tair) * climdata$relhum / 100
    # Mass of vapour leaving the surface over the time step (kg m-2, or mm)
    vapour <- 0.018015 / (8.314 * (Tair + 273.15) * rHa) * (es - ea) * 1000 * dT
    if (vapour > 0) {
      out <- remove_from_top(SWE, vapour)
      SWE <- out$SWE
      sublimation <- out$removed
    } else {
      deposition <- -vapour
      SWE[1] <- SWE[1] + deposition
    }
  }

  # Refreezing. Rain and melt water trickle down and freeze in snow colder than
  # 0 deg C. The heat released warms the snow, and the ice formed makes it denser
  liquid <- snowmod$surfaceWater + rainfall + rainmelt + melt
  refrozen <- 0
  for (i in seq_along(SWE)) {
    if (liquid <= 0) break
    if (SWE[i] > 0 && temp[i] < 0) {
      C <- SnowHeatCapacity(rho[i], temp[i], climdata$pk) * SWE[i] / rho[i]
      freeze <- min(liquid, -C * temp[i] / Lf)
      temp[i] <- min(temp[i] + freeze * Lf / C, 0)
      thickness <- SWE[i] / rho[i]
      SWE[i] <- SWE[i] + freeze
      rho[i] <- min(SWE[i] / thickness, 917)
      liquid <- liquid - freeze
      refrozen <- refrozen + freeze
    }
  }

  # Liquid water that remains is held within the pack up to a fraction of the
  # snow there; the rest drains from its base
  retained_liquid <- min(liquid, liquid_water_frac * sum(SWE))
  throughfall <- liquid - retained_liquid

  # Snow ages, and layers that have melted or sublimated away are dropped
  age <- age + dT / 3600
  keep <- SWE > 1e-10
  SWE <- SWE[keep]
  rho <- rho[keep]
  oldtemp <- oldtemp[keep]
  temp <- temp[keep]
  age <- age[keep]

  # Snow settles as it ages and as more snow builds up above it: its density is
  # the greater of what it has and what its depth and age give
  rho <- pmax(rho, SnowDensFun(sdp, cumsum(SWE / rho), age))
  depth <- SWE / rho

  snowmod$SWE <- SWE
  snowmod$depth <- depth
  snowmod$rho <- rho
  snowmod$oldtemp <- oldtemp
  snowmod$temp <- temp
  snowmod$age <- age
  snowmod$melt <- melt
  snowmod$snowfall <- snowfall
  snowmod$rainfall <- rainfall
  snowmod$throughfall <- throughfall
  snowmod$surfaceWater <- retained_liquid
  snowmod$rainmelt <- rainmelt
  snowmod$refrozen <- refrozen
  snowmod$heattoground <- heattoground
  snowmod$sublimation <- sublimation
  snowmod$deposition <- deposition
  snowmod$layermelt <- layermelt
  snowmod
}

#' @title Regularises snow layer thickness
#' @description Merges layers that have become too thin, as happens when snow
#' melts or falls lightly, and splits layers that have become too thick, so that
#' the heat model has layers of sensible thickness to work with. Layers thinner
#' than 5 mm are merged with a neighbour and layers thicker than 5 cm are split.
#' @param snowmod snow model state as returned by [InitSnowModel()] or updated by model steps
#' @details The water, depth and heat held by the snowpack are unchanged. A
#' merged layer has the density given by the water and depth of the two layers
#' it is made from, and their mean temperature weighted by the water they hold.
#' A snowpack thinner than 5 mm in all is left as a single layer.
#' @return updated \code{snowmod} with adjusted layer structure
#' @export
RegulariseSnowLayers <- function(snowmod) {
  min_thickness <- 0.005
  max_thickness <- 0.05
  # One row for each layer, from the snow surface down, leaving out any that are empty
  L <- data.frame(SWE = snowmod$SWE, depth = snowmod$depth, rho = snowmod$rho,
                  temp = snowmod$temp, oldtemp = snowmod$oldtemp, age = snowmod$age)
  L <- L[L$SWE > 1e-10 & L$depth > 1e-10, ]

  # Merge each layer that is too thin with the thinner of its neighbours. The
  # merged layer holds the water and depth of both, and the heat of both
  i <- 1
  while (nrow(L) > 1 && i <= nrow(L)) {
    if (L$depth[i] < min_thickness) {
      j <- if (i == 1) 2 else if (i == nrow(L)) i - 1 else if (L$depth[i - 1] <= L$depth[i + 1]) i - 1 else i + 1
      ii <- sort(c(i, j))
      w <- L$SWE[ii]
      L[ii[1], ] <- data.frame(SWE = sum(w), depth = sum(L$depth[ii]), rho = sum(w) / sum(L$depth[ii]),
                               temp = weighted.mean(L$temp[ii], w),
                               oldtemp = weighted.mean(L$oldtemp[ii], w),
                               # mean of the logarithm of age, on which albedo depends
                               age = exp(weighted.mean(log(pmax(L$age[ii], 1e-6)), w)))
      L <- L[-ii[2], ]
      i <- max(1, ii[1] - 1) # the merged layer may still be too thin
    } else {
      i <- i + 1
    }
  }

  # Split each layer that is too thick into equal parts
  nsplit <- ceiling(L$depth / max_thickness - 1e-9)
  L <- L[rep(seq_len(nrow(L)), nsplit), ]
  nsplit <- rep(nsplit, nsplit)
  L$SWE <- L$SWE / nsplit
  L$depth <- L$depth / nsplit

  snowmod$SWE <- L$SWE
  snowmod$depth <- L$depth
  snowmod$rho <- L$rho
  snowmod$temp <- L$temp
  snowmod$oldtemp <- L$oldtemp
  snowmod$age <- L$age
  snowmod
}


#' @title Snow thermal conductivity
#' @description Calculates snow thermal conductivity from bulk snow density
#' using an exponential relationship following Djachkova (after Anderson, 2006).
#' It is no greater than the thermal conductivity of solid ice.
#' @param rho Snow density (kg m-3).
#' @return Snow thermal conductivity (W m-1 K-1).
#' @export
SnowThermalConductivity <- function(rho) {
  rhog <- rho / 1000 # density in g cm-3
  # Denser snow has more ice in contact through which heat can pass, so
  # conductivity rises steeply with density, up to that of solid ice
  k <- pmin(0.0442 * exp(5.181 * rhog), 2.2)
  return(k)
}

#' @title Snow volumetric heat capacity
#' @description Calculates the volumetric heat capacity of snow as a mixture of
#' ice and air.
#' @param rho Snow density (kg m-3).
#' @param Tc Snow temperature (deg C).
#' @param pk Air pressure (kPa).
#' @return Snow volumetric heat capacity (J m-3 K-1).
#' @export
SnowHeatCapacity <- function(rho, Tc, pk = 101.3) {
  # Fractions of the volume of the snow that are ice and air
  Vi <- rho / 917 # 917 kg m-3 is the density of ice
  Va <- pmax(1 - Vi, 0)
  # Heat needed to warm a cubic metre of each by one degree (J m-3 K-1). That
  # of ice falls as it gets colder
  Ci <- 1.93e6 + 0.0067e6 * pmin(Tc, 0)
  Ca <- cpair(Tc) * rhoair(Tc, pk)
  # That of the snow is the sum of these, each weighted by its volume
  Cd <- Vi * Ci + Va * Ca
  return(Cd)
}

#' @title Calculates snow surface energy balance
#' @description Calculates the energy available at the snow surface (W m-2):
#' radiation absorbed, less that emitted and the sensible and latent heat given
#' to the air. It is the heat flux into the snow.
#' @param Rabs absorbed radiation at the surface (W m-2)
#' @param Tref reference air temperature (degrees C)
#' @param Tsurface snow surface temperature (degrees C)
#' @param pk air pressure (kPa)
#' @param relhum air relative humidity (Percentage)
#' @param rHa aerodynamic resistance to heat from snow surface to reference height
#' (s m-1), e.g. as returned by [SnowResistance()]
#' @details A snow surface is no warmer than 0 degrees C. If a higher
#' `Tsurface` is given, the heat and radiation leaving the surface are those of
#' a surface at 0 degrees C.
#' @return heat flux into the snow (W m-2)
#' @export
snowsurfaceEB <- function(Rabs, Tref, Tsurface, pk, relhum, rHa) {
  # A snow surface that gains more heat melts and stays at 0 deg C
  Tsurface <- pmin(Tsurface, 0)
  # Radiation absorbed, less that emitted by the snow surface
  sb <- 5.67e-8
  Rnet <- Rabs - 0.99 * sb * (Tsurface + 273.15)^4
  # Sensible heat given to the air
  cp <- cpair(Tref)
  ph <- rhoair(Tref, pk)
  H <- (ph * cp / rHa) * (Tsurface - Tref)
  # Latent heat of the vapour given to the air. The air at the snow surface is saturated
  es <- satvap(Tsurface)
  ea <- satvap(Tref) * (relhum / 100)
  la <- latvap(Tsurface)
  L <- (la * ph / (rHa * pk)) * (es - ea)
  Rnet - H - L
}
# Joins the snow layers and the soil nodes into one column for the snow-soil
# heat solver, from the snow surface down to the lower boundary of the soil.
# Returns the temperatures at the start of the time step, the heat each node
# stores per degree (C, J m-2 K-1) and the conductance between each node and
# the one beneath it (Gint, W m-2 K-1), at the temperatures given in Tcol
.appendSnowSoil <- function(snowmod, soilheatmod, soilp, atmPressure,
                            Tcol = c(snowmod$temp, soilheatmod$Te)) {
  ns <- length(snowmod$depth)
  n <- soilheatmod$n
  Tsnow <- Tcol[seq_len(ns)]
  Tsoil <- Tcol[ns + seq_len(n + 1)]

  # Soil. Each node stores heat for the soil around it, and heat is conducted
  # from each node to the one beneath it across the gap between them. Close to
  # the surface conductivity changes quickly with depth, so the logarithmic
  # mean of the top two nodes is used
  lambda <- thermalConductivity(soilp$Vq, soilp$Vm, soilp$Vo, soilheatmod$wc, soilp$Mc, Tsoil, atmPressure)
  C <- heatCapacity(soilp$Vq, soilp$Vm, soilp$Vo, soilheatmod$wc, Tsoil, atmPressure) * soilheatmod$vol
  Gint <- lambda[seq_len(n)] / soilheatmod$dz[seq_len(n)]
  if (lambda[1] != lambda[2]) {
    Gint[1] <- ((lambda[1] - lambda[2]) / log(lambda[1] / lambda[2])) / soilheatmod$dz[1]
  }

  # Snow. Each layer stores heat through its whole thickness and has its
  # temperature at its middle. Heat passing from the middle of one layer to the
  # middle of the next crosses half of each. From the middle of the lowest layer
  # to the soil surface it crosses half of that layer
  if (ns > 0) {
    d <- snowmod$depth
    k <- SnowThermalConductivity(snowmod$rho)
    halfR <- 0.5 * d / k # resistance of half of each layer (m2 K W-1)
    Gsnow <- 1 / (halfR + c(halfR[-1], 0))
    C <- c(SnowHeatCapacity(snowmod$rho, Tsnow, atmPressure) * d, C)
    Gint <- c(Gsnow, Gint)
  }
  list(ns = ns, nz = ns + n + 1,
       T_old = c(snowmod$oldtemp, soilheatmod$oldTe),
       C = C, Gint = Gint)
}
# Returns the temperatures of the joined column to the snow layers and the soil
# nodes, with the heat flux into the top of the column (Gsurface) and that into
# the soil at its surface (Gsoil), both in W m-2
.splitSnowSoil <- function(T_new, snowmod, soilheatmod, Gsurface, Gsoil) {
  ns <- length(snowmod$depth)
  snowmod$temp <- T_new[seq_len(ns)]
  soilheatmod$Te <- T_new[ns + seq_along(soilheatmod$Te)]
  snowmod$Gflux <- Gsurface
  soilheatmod$Gflux <- Gsoil
  list(snowmod = snowmod, soilheatmod = soilheatmod)
}

#' @title Updates snow and soil temperatures using the coupled snow-soil heat model
#' @description Solves heat conduction through the snow and the soil beneath it
#' as one column over one time step. The energy balance of the surface, from
#' [snowsurfaceEB()] where there is snow and [soilsurfaceEB()] where there is
#' none, is the upper boundary, and the deepest soil node is held at a fixed
#' temperature. Heat is conducted according to [SnowThermalConductivity()] and
#' [thermalConductivity()] and stored according to [SnowHeatCapacity()] and
#' [heatCapacity()], with temperatures at the end of the time step found using
#' [ThomasBoundaryCondition()]. With no snow it is [SoilHeatModel()].
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
#' as returned by [SnowResistance()] (s m-1)
#' @param dT time step (s)
#' @param maxNrIterations maximum number of iterations
#' @param tolerance convergence criterion for maximum change in temperature
#' between successive iterations (degrees C)
#' @details The model steps forward from the temperatures in `snowmod$oldtemp`
#' and `soilheatmod$oldTe`, which are those at the start of the time step, and
#' returns those at its end in `snowmod$temp` and `soilheatmod$Te`. The old
#' temperatures are left as they are; to move on to the next time step, set them
#' to the new ones.
#'
#' Snow that reaches 0 degrees C is held there while heat is conducted through
#' the column, and the heat it gains beyond that is set aside for melting. That
#' heat is returned as a temperature above 0 degrees C: the temperature the
#' layer would have reached had it not melted. Pass the result to
#' [SnowMassBalance()], which melts snow with the heat held above 0 degrees C
#' and returns the snow to 0 degrees C.
#' @return a list comprising
#' \describe{
#'   \item{snowmod}{updated snow model state. `temp` holds the temperature of
#'   each layer at the end of the time step and `Gflux` the heat flux into the
#'   snow at its surface (W m-2)}
#'   \item{soilheatmod}{updated soil heat model state. `Te` holds the
#'   temperature of each node at the end of the time step, `Gflux` the heat
#'   flux into the soil at its surface (W m-2) and `iters` the number of
#'   iterations used}
#' }
#' @export
SnowSoilHeatModel <- function(snowmod, soilheatmod, soilp,
                              Rabs, Tref, relhum, atmPressure, rHa, dT = 3600,
                              maxNrIterations = 100, tolerance = 1e-3) {
  prof <- .appendSnowSoil(snowmod, soilheatmod, soilp, atmPressure)
  ns <- prof$ns
  nz <- prof$nz
  T_old <- prof$T_old
  # Temperatures at the end of the time step are found by iteration, starting
  # from the latest estimate. The deepest node keeps its temperature
  T_new <- c(snowmod$temp, soilheatmod$Te)
  T_new[nz] <- T_old[nz]
  idx <- seq_len(nz - 1) # the nodes solved

  # Energy available at the surface, which is snow where there is any
  surfaceEB <- function(Tsurface) {
    if (ns > 0) {
      snowsurfaceEB(Rabs, Tref, Tsurface, atmPressure, relhum, rHa)
    } else {
      soilsurfaceEB(soilp, Rabs, Tref, Tsurface, atmPressure, relhum, rHa, soilheatmod$wc[1])
    }
  }

  melting <- rep(FALSE, nz - 1) # snow layers held at 0 deg C
  spare <- numeric(nz - 1)      # heat each of them gains beyond that (W m-2)
  nrIterations <- 0
  maxdT <- Inf
  settled <- FALSE
  while ((maxdT > tolerance || !settled) && nrIterations < maxNrIterations) {
    # Heat stored per degree by each node and conductance between each node and
    # the one beneath it, at the current temperatures. After twenty iterations
    # they are held fixed, which lets the solve settle when a soil node sits at
    # the freezing point, where the heat capacity of soil jumps
    if (nrIterations < 20) prof <- .appendSnowSoil(snowmod, soilheatmod, soilp, atmPressure, T_new)
    C <- prof$C[idx]
    below <- prof$Gint
    above <- c(0, below[-(nz - 1)])

    # Heat flux into the column at its surface, and how it changes as the
    # surface warms: a warmer surface emits more radiation and gives more heat
    # to the air. For snow the change is taken towards colder temperatures,
    # since a snow surface at 0 deg C can cool but cannot warm
    qsurface <- surfaceEB(T_new[1])
    h <- if (ns > 0) -0.01 else 0.01
    dqsurface <- (surfaceEB(T_new[1] + h) - qsurface) / h

    # Energy balance of each node over the time step: the heat it stores equals
    # the heat conducted to it from the nodes above and below, at the
    # temperatures of the end of the time step
    aa <- -above
    cc <- -below
    bb <- C / dT + above + below
    dd <- C / dT * T_old[idx]
    # The surface node also gains the heat flux into the column
    bb[1] <- bb[1] - dqsurface
    dd[1] <- dd[1] + qsurface - dqsurface * T_new[1]
    # The node above the lower boundary exchanges heat with a fixed temperature
    cc[nz - 1] <- 0
    dd[nz - 1] <- dd[nz - 1] + below[nz - 1] * T_new[nz]
    # Melting snow stays at 0 deg C
    aa[melting] <- 0
    cc[melting] <- 0
    bb[melting] <- 1
    dd[melting] <- 0

    Tsolved <- ThomasBoundaryCondition(aa, bb, cc, dd, numeric(nz - 1))$x
    maxdT <- max(abs(Tsolved - T_new[idx]))
    T_new[idx] <- Tsolved

    # Heat that each melting layer gains beyond what holds it at 0 deg C: that
    # conducted to it and, at the surface, the energy available there, less the
    # heat that brought it to 0 deg C
    spare <- above * (c(0, T_new[idx])[idx] - T_new[idx]) + below * (T_new[-1] - T_new[idx]) -
      C / dT * (T_new[idx] - T_old[idx])
    spare[1] <- spare[1] + surfaceEB(T_new[1])
    # Snow starts to melt when it would otherwise pass 0 deg C, and stops when
    # it loses more heat than it gains
    wasmelting <- melting
    melting[seq_len(ns)] <- (melting[seq_len(ns)] & spare[seq_len(ns)] >= 0) |
      (!melting[seq_len(ns)] & T_new[seq_len(ns)] > 0)
    settled <- identical(melting, wasmelting)
    nrIterations <- nrIterations + 1
  }

  # Heat flux into the column at its surface, and that passing into the soil:
  # from the lowest snow layer where there is snow
  Gsurface <- surfaceEB(T_new[1])
  Gsoil <- Gsurface
  if (ns > 0) Gsoil <- prof$Gint[ns] * (T_new[ns] - T_new[ns + 1])
  # The heat set aside for melting is returned as the temperature it would
  # have brought each melting layer to
  T_new[idx][melting] <- (spare * dT / C)[melting]
  out <- .splitSnowSoil(T_new, snowmod, soilheatmod, Gsurface, Gsoil)
  out$soilheatmod$iters <- nrIterations
  return(out)
}
#' @title Canopy snow interception
#' @description Calculates the snow caught by a vegetation canopy during a time
#' step and the snow that falls from the canopy to the ground, following
#' Hedstrom and Pomeroy (1998). The canopy catches a smaller share of the
#' snowfall the more snow it already holds, up to the most it can hold, and the
#' snow it holds falls off gradually.
#' @param hgt Canopy height (m).
#' @param pai Plant area index (dimensionless).
#' @param uf Friction velocity (m s^-1).
#' @param prec Snowfall during the time step (mm snow water equivalent).
#' @param tc Air temperature (degrees C).
#' @param Li Snow load on the canopy at the start of the time step
#'   (mm snow water equivalent).
#' @param Sh Snow load that a unit area of branch can hold (kg m^-2). Hedstrom
#'   and Pomeroy give 6.6 for pine and 5.9 for spruce.
#' @param dT time step (s)
#' @return a list of the following, each in mm snow water equivalent:
#' \describe{
#'  \item{intercepted}{snowfall caught by the canopy during the time step}
#'  \item{unloaded}{snow falling from the canopy to the ground during the time step}
#'  \item{load}{snow load on the canopy at the end of the time step}
#' }
#' @details
#' Wind carries falling snow sideways, so that it meets more of the canopy than
#' snow falling straight down. The mean wind speed within the canopy is taken
#' from [windprofile_below()] for foliage spread evenly with height, and snow is
#' taken to fall at 0.8 m/s. The most snow the canopy can hold rises with plant
#' area and is greater for the light snow that falls in cold air. Snow held by
#' the canopy falls off at a rate that leaves 0.678 of it after a week.
#' Sublimation and melt of the snow held are not included. Very small values of
#' `hgt` and `pai` are bounded below at 0.001.
#' @references Hedstrom NR, Pomeroy JW (1998) Measurements and modelling of snow
#' interception in the boreal forest. Hydrological Processes 12: 1611-1625.
#' @export
CanopySnowInterception <- function(hgt, pai, uf, prec, tc, Li, Sh = 6.2, dT = 3600) {
  hgt <- pmax(hgt, 0.001)
  pai <- pmax(pai, 0.001)
  # Mean wind speed within the canopy: wind at canopy top, slowed by the foliage
  # on the way down
  Be <- pmin(sqrt(0.003 + 0.15 * pai), 0.3) # ratio of friction velocity to wind at canopy top
  uzm <- mapply(function(h, p, u) mean(windprofile_below(h, rep(p / 20, 20), u)), hgt, pai, uf / Be)

  # Most snow the canopy can hold (mm). Branches hold more of the light snow
  # that falls in cold air than of dense, wet snow
  rhos <- 67.92 + 51.25 * exp(tc / 2.59) # density of the falling snow (kg m-3)
  S <- Sh * (0.27 + 46 / rhos)           # snow held per unit area of branch (kg m-2)
  Lstr <- S * pai

  # Share of the falling snow that meets the canopy. Wind makes snow fall at an
  # angle, lengthening its path through the foliage
  Z <- atan(uzm / 0.8)                 # angle of fall from the vertical, for snow falling at 0.8 m/s
  Cp <- 1 - exp(-0.5 * pai / cos(Z))   # foliage facing all directions equally

  # Snow caught. The canopy catches less of each further millimetre of snowfall
  # as it fills
  intercepted <- pmax(Lstr - Li, 0) * (1 - exp(-Cp * prec / Lstr))

  # Snow held by the canopy falls off at a steady fractional rate
  U <- -log(0.678) / (7 * 86400) # per second
  unloaded <- (Li + intercepted) * (1 - exp(-U * dT))
  list(intercepted = intercepted, unloaded = unloaded, load = Li + intercepted - unloaded)
}
