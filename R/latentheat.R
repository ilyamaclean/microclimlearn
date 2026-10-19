#' Calculate saturated vapour pressure
#'
#' @param tc temperature (degrees C)
#' @returns saturated vapour pressure (kPa).
#' @export
satvap <- function(tc) {
  neg <- tc < 0
  es <- rep(NA_real_, length(tc))
  # Saturated vapour pressure over water above freezing
  es[!neg] <- 0.61078 * exp(17.27 * tc[!neg] / (tc[!neg] + 237.3))
  # Saturated vapour pressure over ice below freezing, which is lower than over water
  es[neg] <- 0.61078 * exp(21.875 * tc[neg] / (tc[neg] + 265.5))
  return(es)
}
#' Calculate latent heat of vapourisation
#'
#' Above freezing this is the energy needed to evaporate liquid water. Below
#' freezing, water vapour comes from ice, which takes the energy of melting as
#' well. Between 0 and -1 degrees C the value passes gradually from one to the
#' other.
#' @param Ts surface temperature (degrees C)
#' @returns latent heat of vapourisation (J / mol)
#' @export
latvap <- function(Ts) {
  # Energy needed to evaporate a mole of liquid water, which falls as temperature rises
  vaporisation <- 45068.7 - 42.8428 * Ts
  # Below freezing water vapour comes from ice, which takes the energy of melting as well
  sublimation <- 51078.69 - 4.338 * Ts - 0.06367 * Ts^2
  # Water freezes gradually over the first degree below zero, so the share of
  # vapour that comes from ice rises from none at 0 to all of it at -1
  frozen <- pmin(pmax(-Ts, 0), 1)
  la <- vaporisation + frozen * (sublimation - vaporisation)
  return(la)
}
#' Creates a list of vegetation parameter inputs for a given plant function type
#' @param PFT a plant function type as listed in `PFTtable`
#' @returns a list of vegetation parameters.
#' @export
createplant_inputs <- function(PFT) {

  s <- match(PFT, names(PFTtable))
  if (is.na(s)) stop("PFT not found in PFTtable")

  # Values for the plant functional type, converted from the units of the table
  # to those the functions use
  v <- PFTtable[, s] * PFTtable$multiplier
  names(v) <- PFTtable$varname

  plant_inputs <- list(
    h       = v[["h"]],       # vegetation height (m)
    pai     = v[["pai"]],     # plant area index
    x       = v[["x"]],       # Campbell x parameter
    lref    = v[["lref"]],    # leaf reflectance (shortwave)
    ltra    = v[["ltra"]],    # leaf transmittance (shortwave)
    lrefp   = v[["lrefp"]],   # leaf reflectance (PAR)
    ltrap   = v[["ltrap"]],   # leaf transmittance (PAR)
    Vcmax25 = v[["Vcmx25"]],  # maximum rate of Rubisco carboxylation at 25 deg C (mol / m^2 / s)
    Tup     = v[["Tup"]],     # upper temperature limit of photosynthesis (deg C)
    Tlw     = v[["Tlow"]],    # lower temperature limit of photosynthesis (deg C)
    alpha   = v[["alpha"]],   # quantum efficiency (mol CO2 / mol PAR)
    rpmin   = v[["rpmin"]],   # minimum whole-plant hydraulic resistance (m^2 s MPa / mol)
    fd      = v[["fd"]],      # dark respiration as a fraction of Vcmax
    psi50   = v[["psi50"]],   # water potential at which half of hydraulic conductance is lost (MPa)
    apsi    = v[["apsi"]],    # shape of the loss of hydraulic conductance as water potential falls
    len     = v[["len"]],     # leaf length (m)
    wid     = v[["wid"]],     # leaf width (m)
    em      = v[["vegem"]],   # vegetation emissivity
    pTAW    = v[["pTAW"]],    # fraction of available soil water that can be used before stress
    root50  = v[["root50"]],  # depth above which half of all roots lie (m)
    root95  = v[["root95"]]   # depth above which 95 percent of all roots lie (m)
  )
  return(plant_inputs)
}
# Shape of the loss of hydraulic conductance as water potential falls, from the
# water potential at which half of conductance is lost
.psi50toa <- function(psi50){
  stem_Slope <- 65.15*(-psi50)^(-1.25)
  a <- -4*stem_Slope/100*psi50
  return(a)
}
# Rate of assimilation where two limiting rates meet: a little below the lesser
# of the two, passing smoothly from one to the other (Collatz et al. 1991)
.colimited <- function(a, b, beta) {
  ((a + b) - sqrt(pmax((a + b)^2 - 4 * beta * a * b, 0))) / (2 * beta)
}
#' Calculate leaf stomatal conductance
#' @param Ca atmospheric carbon dioxide concentration (ppm)
#' @param Rswabs flux density of radiation absorbed by leaf (W/m^2), of which 48 percent is taken to be PAR
#' @param tair temperature of air (deg C)
#' @param tleaf temperature of leaf (deg C)
#' @param rh relative humidity of air (percentage)
#' @param pk atmospheric pressure (kPa)
#' @param psi_r mean water potential in root zone (MPa)
#' @param plant_inputs a list of vegetation parameters as returned by [createplant_inputs()]
#' @param z height above ground
#' @param C3 optional logical indicating whether vegetation has C3 or C4 photosynthetic pathway
#' @returns stomatal conductance to water vapour (mol / m^2 / s)
#' @export
stomatalcond_calc <- function(Ca, Rswabs, tair, tleaf, rh, pk, psi_r, plant_inputs, z, C3 = TRUE) {
  # Photosynthetic capacity at leaf temperature. It doubles with every 10 degrees
  # of warming, and falls away towards the lower and upper temperature limits
  Vcmax <- plant_inputs$Vcmax25 * 2^(0.1 * (tleaf - 25)) /
    ((1 + exp(0.3 * (tleaf - plant_inputs$Tup))) * (1 + exp(0.3 * (plant_inputs$Tlw - tleaf))))
  Rd <- plant_inputs$fd * Vcmax # dark respiration

  # CO2 concentration at which photosynthesis balances photorespiration (Pa)
  P <- pk * 1000 # air pressure (Pa)
  Oa <- 0.2095 * P # partial pressure of oxygen (Pa)
  photocomp <- Oa / (2 * 2600 * 0.57^(0.1 * (tleaf - 25)))

  ca <- Ca * pk / 1000 # partial pressure of CO2 in the air (Pa)
  IPAR <- Rswabs * (0.48 / 0.219) * 1e-6 # absorbed PAR (mol photons / m^2 / s)
  light <- plant_inputs$alpha * IPAR # rate at which light can drive assimilation

  if (C3) {
    # Affinity of Rubisco for CO2, allowing for competition from oxygen (Pa)
    Kc <- 30 * 2.1^(0.1 * (tleaf - 25))
    Ko <- 30000 * 1.2^(0.1 * (tleaf - 25))
    Km <- Kc * (1 + Oa / Ko)
    We <- 0.5 * Vcmax # rate limited by the transport of the products of photosynthesis
    # CO2 concentration inside the leaf above which more CO2 adds nothing: where
    # the rate limited by Rubisco meets the co-limited rate of light and transport
    Wcol <- .colimited(We, light, 0.93)
    cicol <- (-Vcmax * photocomp - Km * Wcol) / (Wcol - Vcmax)
    # Net assimilation at a given CO2 concentration inside the leaf: the
    # co-limited rate of the rates limited by Rubisco, light and transport
    assimilation <- function(ci) {
      Wc <- pmax(Vcmax * (ci - photocomp) / (ci + Km), 0)
      Wl <- pmax(light * (ci - photocomp) / (ci + 2 * photocomp), 0)
      .colimited(.colimited(Wc, Wl, 0.83), We, 0.93) - Rd
    }
  } else {
    # CO2 concentration inside the leaf above which more CO2 adds nothing: where
    # the rate limited by CO2 meets the co-limited rate of Rubisco and light
    k <- 2e4
    Wcol <- .colimited(Vcmax, light, 0.83)
    cicol <- (Wcol * P) / (k * Vcmax)
    # Net assimilation at a given CO2 concentration inside the leaf: the
    # co-limited rate of that and the rate limited by CO2
    assimilation <- function(ci) {
      .colimited(Wcol, k * Vcmax * ci / P, 0.93) - Rd
    }
  }

  # Gain in assimilation per unit mole fraction of CO2 from opening the stomata, from
  # the point at which CO2 stops limiting to CO2 inside the leaf equal to that of the air
  Aca <- assimilation(ca)
  Acol <- assimilation(cicol)
  dadc <- P * (Aca - Acol) / (ca - cicol)

  # Water potential of the leaf before dawn: that of the root zone less the pull
  # of gravity on the water column below the leaf (MPa)
  rhow <- 1000 * (1 - ((tleaf + 288.9414) * (tleaf - 3.9863)^2) / (508929.2 * (tleaf + 68.12963)))
  psi_pd <- psi_r - z * 9.81 * rhow * 1e-6

  # Fraction of maximum hydraulic conductance that remains at this water
  # potential, and how steeply it is lost if water potential falls further
  psi_50 <- plant_inputs$psi50
  apsi <- plant_inputs$apsi
  if (apsi < 0) apsi <- .psi50toa(psi_50)
  K_psi_pd <- 1 / (1 + (psi_pd / psi_50)^apsi)
  psi_50f <- (psi_pd + psi_50) / 2
  K_50f <- 1 / (1 + (psi_50f / psi_50)^apsi)
  dKdpKi <- ((K_psi_pd - K_50f) / (psi_pd - psi_50f)) * (1 / K_psi_pd)

  # Resistance of the plant to the flow of water from root to leaf, which rises
  # as the xylem loses conductance
  rp <- plant_inputs$rpmin / K_psi_pd

  # Cost in lost water supply of the carbon gained by opening the stomata, which
  # rises as the air dries. A leaf at or below the dew point has no deficit
  DD <- pmax(satvap(tleaf) - satvap(tair) * (rh / 100), 1e-9) / pk # vapour pressure deficit (mol / mol)
  zeta <- 2 / (dKdpKi * rp * 1.6 * DD)

  # Conductance at which the gain in carbon balances the loss of water supply
  mu <- pmax(1 + (4 * zeta) / dadc, 1)
  gs <- 0.5 * dadc * (sqrt(mu) - 1)
  # Where more CO2 brings no gain, conductance just brings the CO2 concentration
  # inside the leaf to the point at which it stops limiting
  flat <- rep_len(dadc <= 0, length(gs))
  gs[flat] <- rep_len(pmax(Acol * P / (ca - cicol), 0), length(gs))[flat]

  # Stomata have an upper limit to conductance and are closed in the dark
  gsmax <- plant_inputs$Vcmax25 * 1e4
  gs <- pmin(gs, gsmax)
  gs[rep_len(Rswabs <= 0, length(gs))] <- 0

  # Water vapour diffuses through the stomata 1.6 times as fast as CO2
  return(1.6 * gs)
}
#' Calculate bulk surface stomatal resistance
#' @param solp a list of solar positions as returned by [sunposition()]
#' @param Ca atmospheric carbon dioxide concentration (ppm)
#' @param Rsw flux density of shortwave radiation on the horizontal (W/m^2)
#' @param Rdif flux density of diffuse radiation on the horizontal (W/m^2)
#' @param tair temperature of air (deg C)
#' @param tcanopy temperature of canopy (deg C)
#' @param rh relative humidity of air (percentage)
#' @param pk atmospheric pressure (kPa)
#' @param psi_r mean water potential in root zone (MPa)
#' @param plant_inputs a list of vegetation parameters as returned by [createplant_inputs()]
#' @param Lfrac fraction of plant area that is living leaf (default 0.8)
#' @param C3 optional logical indicating whether vegetation has C3 or C4 photosynthetic pathway
#' @returns bulk surface stomatal resistance (s/m)
#' @export
bulkstomatalresist_calc <- function(solp, Ca, Rsw, Rdif, tair, tcanopy, rh, pk, psi_r, plant_inputs, Lfrac = 0.8, C3 = TRUE) {
  # Extract parameters from plant inputs
  z <- plant_inputs$h / 2
  x <- plant_inputs$x
  lrefp <- plant_inputs$lrefp
  ltrap <- plant_inputs$ltrap
  PAI <-  plant_inputs$pai
  om <- lrefp + ltrap
  # Calculate canopy extinction coefficient
  zen <- solp$zen
  azi <- solp$azi
  zen[zen > 90] <- 90
  if (x == 1) {
    k <- 1 / (2 * cos(zen * pi/180))
  } else if (x == 0) {
    k <- tan(zen * pi/180)
  } else if (is.infinite(x)) {
    k<-rep(1, length(zen))
  } else {
    k <- sqrt((x^2 + (tan(zen * pi/180)^2)))/(x + 1.774 * (x + 1.182)^(-0.733))
  }
  # Sunlit leaf area
  L_sun <- (1 - exp(-k* PAI * Lfrac)) / k
  L_shade <- (PAI * Lfrac) - L_sun
  # Absorbed - sun and shaded parts
  Rshade_abs <- Rdif * ((1 - exp(-PAI)) / PAI) * (1 - om)
  Rsun_abs <- (Rsw - Rdif) * k  * (1 - om) + Rshade_abs
  # Calculate stomatal conductance
  gs_sun <- stomatalcond_calc(Ca, Rsun_abs, tair, tcanopy, rh, pk, psi_r, plant_inputs, z, C3)
  gs_shade <- stomatalcond_calc(Ca, Rshade_abs, tair, tcanopy, rh, pk, psi_r, plant_inputs, z, C3)
  Gs <- gs_sun * L_sun + gs_shade * L_shade
  # Convert to bulk surface resistance
  ph <- rhoair(tair, pk)
  rS <- ph / Gs
  return(rS)
}
#' Dew-point temperature
#'
#' Calculate dew-point temperature (degrees C) from air temperature (degrees C) and
#' relative humidity (\%). Uses Magnus equations with coefficients for
#' water and ice, applying the ice formulation when dew point is below 0 degrees C.
#'
#' @param tc Air temperature (degrees C).
#' @param rh Relative humidity (\%).
#'
#' @return Dew-point temperature (degrees C).
#'
#' @seealso [satvap()]
#' @export
dewpoint <- function(tc, rh) {
  ea <- satvap(tc) * rh / 100 # vapour pressure of the air (kPa)
  # Temperature at which this vapour pressure saturates the air, over water
  gamma <- log(ea / 0.61078)
  td <- (237.3 * gamma) / (17.27 - gamma)
  # Below freezing the air saturates over ice
  ice <- td < 0
  td[ice] <- (265.5 * gamma[ice]) / (21.875 - gamma[ice])
  return(td)
}
#' Calculate relative humidity within a canopy
#'
#' Relative humidity at any height within a canopy, from the latent heat
#' released by the foliage of each layer and the humidity of the air at canopy
#' top and of the ground surface. Vapour pressure is calculated as the sum of
#' the far field, as returned by [farfield()], and the near field, as returned by
#' [nearfield()], and converted to relative humidity at the temperature of the air.
#' @param Lz vector of the flux of latent heat from the foliage of each canopy
#' layer, per unit ground area (W/m^2), ordered from canopy bottom to top
#' @param hgt height of canopy (m)
#' @param pai plant area index of the canopy as a whole
#' @param uf friction velocity (m/s) as returned by [windfric()] or [solve_wholecanopy()]
#' @param rhtop relative humidity of the air at canopy top (percentage)
#' @param ttop air temperature at canopy top (deg C)
#' @param tground temperature of the ground surface (deg C)
#' @param soilrh effective relative humidity of the soil surface (percentage)
#' @param tair air temperature within the canopy (deg C), as returned by
#' [temp_below()]. A single value or one for each height at which values are
#' returned.
#' @param pk atmospheric pressure (kPa)
#' @param LL Obukhov length (m) as returned by [Obukhov()]. Infinite under neutral conditions.
#' @param z optional vector of heights above ground (m) at which values are
#' returned, between the ground surface (0) and canopy top (`hgt`). The default
#' is the middle of each canopy layer.
#' @returns a vector of relative humidities (percentage) at each height `z`
#' @seealso [temp_below()] for air temperature
#' @export
relhum_below <- function(Lz, hgt, pai, uf, rhtop, ttop, tground, soilrh, tair, pk = 101.3, LL = Inf,
                         z = NULL) {
  # Vapour pressure of the air at canopy top and at the ground surface (kPa)
  etop <- satvap(ttop) * rhtop / 100
  eground <- satvap(tground) * soilrh / 100
  # Vapour mixed gradually through the canopy, and vapour still concentrated
  # close to the foliage that released it
  far <- farfield(Lz, hgt, pai, uf, etop, eground, tair, pk, LL, vapour = TRUE, z = z)
  near <- nearfield(Lz, hgt, pai, uf, tair, pk, LL, vapour = TRUE, z = z)
  # As a percentage of what air at that temperature can hold
  rh <- 100 * (far + near) / satvap(tair)
  rh <- pmin(pmax(rh, 0), 100)
  return(rh)
}
