#' @title Calculates the projected silhouette area of an ellipsoidal body
#' @description Computes the surface area, volume, projected silhouette area,
#' and solar interception coefficient of an ellipsoid representing an organism
#' or object with specified height, width, and length. The silhouette area is
#' calculated for a specified solar position and body orientation, or under
#' simplified assumptions about orientation.
#' @param zen solar zenith angle (degrees)
#' @param azi solar azimuth angle (degrees)
#' @param height body height along the vertical axis (cm)
#' @param width body width along the lateral axis (cm)
#' @param len body length along the longitudinal axis (cm)
#' @param adir body azimuth direction (degrees). Only used when
#' \code{position = "fixed"}. Default is 0.
#' @param atilt body tilt from vertical (degrees). Only used when
#' \code{position = "fixed"} or \code{position = "randomdir"}. Default is 0.
#' @param position assumption about body orientation relative to the sun. Must
#' be one of:
#' \describe{
#'   \item{\code{"fixed"}}{body orientation is fixed according to
#'   \code{adir} and \code{atilt}}
#'   \item{\code{"max"}}{body orientation is assumed to maximise projected
#'   silhouette area}
#'   \item{\code{"min"}}{body orientation is assumed to minimise projected
#'   silhouette area}
#'   \item{\code{"randomdir"}}{body azimuth direction is assumed random, but
#'   tilt is fixed by \code{atilt}}
#'   \item{\code{"random"}}{both body azimuth direction and tilt are assumed
#'   random}
#' }
#' @return a list comprising
#' \describe{
#'   \item{SurfArea}{total ellipsoid surface area (m^2), approximated using
#'   Knud Thomsen's formula}
#'   \item{Volume}{ellipsoid volume (m^3)}
#'   \item{silA}{projected silhouette area normal to the solar beam (m^2)}
#'   \item{SolarCoef}{ratio of projected silhouette area to total surface area
#'   (-)}
#' }
#' @details The three body dimensions are treated as the full axis lengths of an
#' ellipsoid and are converted internally from centimetres to metres before
#' calculation. Surface area is approximated because no exact closed-form
#' solution exists for a general triaxial ellipsoid.
#' @export
silhouette <- function(zen, azi, height, width, len, adir = 0, atilt = 0.0, position = "fixed") {

  # Compute semi-axis (converted from cm to m)
  AA <- len / 200
  BB <- width / 200
  CC <- height / 200
  # Convert directions to radians
  zenr <- zen * pi / 180
  azir <- azi * pi / 180
  adirr <- adir * pi / 180
  atiltr <- atilt * pi / 180
  # Calculate surface area using Knud Thomsens formula for ellipsoid
  P <- 1.6075 # Approximation constant
  Area <- 4 * pi * (((AA^P * BB^P) + (AA^P * CC^P) + (BB^P * CC^P)) / 3)^(1 / P)
  Volume <- (4 / 3) * pi * AA * BB * CC
  if (position == "max") { # assume solar radiation maximised
    silA <- pi * max(c(BB * CC, AA * CC, AA * BB))
  }
  else if (position == "min") { # assume solar radiation minimised
    silA <- pi * min(c(BB * CC, AA * CC, AA * BB))
  }
  else if (position == "randomdir") { # direction (horizontal rotation) assumed random
    phi <- acos(cos(zenr) * cos(-atiltr) + sin(zenr) * sin(-atiltr) * cos(azir))
    M <- sin(phi)
    N <- cos(phi)
    silA <- pi * sqrt(0.5 * M^2 * (BB^2 * CC^2 + CC^2 * AA^2) + N^2 * AA^2 * BB^2)
  }
  else if (position == "random") { # direction and tilt assumed random
    silA <- pi * sqrt((AA^2 * BB^2 + AA^2 * CC^2 + BB^2 * CC^2) / 3)
  }
  else if (position == "fixed") {
    theta <- (adirr - azir)
    phi <- acos(cos(zenr) * cos(-atiltr) + sin(zenr) * sin(-atiltr) * cos(theta))
    L <- cos(theta) * sin(phi)
    M <- sin(theta) * sin(phi)
    N <- cos(phi)
    silA <- pi * sqrt(L^2 * BB^2 * CC^2 +
                        M^2 * CC^2 * AA^2 +
                        N^2 * AA^2 * BB^2)
  }
  else stop("Position not recongised")
  out <- list(SurfArea = Area, Volume = Volume, silA = silA, SolarCoef = silA / Area)
  return (out)
}
#' @title Calculates a characteristic dimension based on projected area
#' @description Computes a characteristic dimension of an ellipsoidal body as
#' the ratio of volume to projected silhouette area in the direction of wind
#' flow. The projected area is obtained from [silhouette()], with body
#' orientation specified directly or inferred from simple assumptions.
#' @param wdir wind direction (degrees)
#' @param height body height along the vertical axis (cm)
#' @param width body width along the lateral axis (cm)
#' @param len body length along the longitudinal axis (cm)
#' @param adir body azimuth direction (degrees). Only used when
#' \code{position = "fixed"}. Default is 0.
#' @param atilt body tilt from vertical (degrees). Only used when
#' \code{position = "fixed"}. Default is 0.
#' @param position assumption about body orientation. Must be one of:
#' \item{\code{"fixed"}}{body orientation is fixed according to
#'   \code{adir} and \code{atilt}}
#'   \item{\code{"max"}}{body orientation is assumed to maximise projected
#'   silhouette area}
#'   \item{\code{"min"}}{body orientation is assumed to minimise projected
#'   silhouette area}
#'   \item{\code{"randomdir"}}{body azimuth direction is assumed random, but
#'   tilt is fixed by \code{atilt}}
#'   \item{\code{"random"}}{both body azimuth direction and tilt are assumed
#'   random}
#' @return a numeric value giving the characteristic dimension (m), defined as
#' volume divided by projected silhouette area
#' @details The characteristic dimension is calculated as:
#' \deqn{d = V / A_p}
#' where \eqn{V} is ellipsoid volume and \eqn{A_p} is projected area normal to
#' the flow direction. The projected area is always evaluated for a horizontal
#' beam (i.e. \eqn{90^\circ} zenith angle) aligned with the wind direction.
#' Solar angles (\code{zen}, \code{azi}) are only used to determine body
#' orientation when \code{position} is set to \code{"max"} or \code{"min"}.
#' @export
chardim <- function(wdir, height, width, len,
                    adir = 0.0, atilt = 0.0, position = "fixed") {
  sa <- silhouette(
    zen = 90,
    azi = wdir,
    height = height,
    width = width,
    len = len,
    adir = adir,
    atilt = atilt,
    position = position
  )
  d <- sa$Volume / sa$silA
  return(d)
}
#' @title Calculates boundary layer resistance to heat transfer for an animal
#' @description Computes the aerodynamic resistance to heat transfer between an
#' animal and the surrounding air, accounting for both forced and free
#' convection using standard dimensionless relationships.
#' @param tair air temperature (degrees C)
#' @param dT initial estimate of temperature difference between the animal surface and air
#' (degrees C)
#' @param uz wind speed at the animal scale (m s^-1)
#' @param d characteristic dimension of the animal (m), typically defined as
#' volume divided by projected area
#' @param rHmax maximum allowed boundary layer resistance (s m^-1). Default is
#' 300.
#' @return a numeric value giving the boundary layer resistance to heat transfer
#' (s m^-1)
#' @details Thermal diffusivity and kinematic viscosity are calculated as
#' functions of air temperature. The Nusselt number for forced convection is
#' estimated using empirical relationships across Reynolds number regimes,
#' while free convection is estimated from Grashof and Prandtl numbers. The two
#' contributions are combined using a cubic mean to represent mixed convection.
#' @export
animalrHa <- function(tair, dT, uz, d, rHmax = 300.0) {
  Tk <- tair + 273.15
  # Thermal diffusivity
  Kh <- (1.6667e-10 * Tk^2 + 2.9935e-8 * Tk - 1.7128e-6)
  # Kinematic viscosity
  v <- 1.326e-5 * (Tk / 273.15)^1.5 * (393.55 / (Tk + 120.0))
  # Reynolds number
  Re <- (uz * d) / v
  # Prandtl number
  Pr <- v / Kh
  # Nusselt number (forced convection)
  if (Re > 2e5) {
    Nuf <- 0.37 * Re^0.6 * Pr^(1/3) + 9.08
  } else if (Re > 1000) {
    Nuf <- 2.0 + (0.48 * Re^0.6 - 11.31) * Pr^(1/3)
  } else {
    Nuf <- 2.0 + 0.6 * sqrt(Re) * Pr^(1/3)
  }
  # Grashof number (g = 9.81 m s^-2)
  Gr <- (9.81 * d^3 * dT) / (Tk * v^2)
  # Nusselt number (free convection)
  Nun <- (0.825 + (0.387 * (Gr * Pr)^(1/6)) /
            (1 + (0.492 / Pr)^(9/16))^(8/27))^2
  # Mixed convection
  Nu <- (Nuf^3 + Nun^3)^(1/3)
  # Boundary layer resistance
  rHa <- d / (Kh * Nu)
  if (rHa > rHmax) rHa <- rHmax
  return(rHa)
}

#' Calculates metabolic heat production as a function of body mass and temperature.
#'
#' @param volume Animal volume (m^3)
#' @param rho Animal density (kg m^-3)
#' @param Q10 Temperature sensitivity coefficient
#' @param a0 Normalisation constant
#' @param b Mass-scaling exponent
#' @param Tref Reference temperature (°C)
#' @param Tbody Body temperature (°C)
#'
#' @return Numeric vector of metabolic rate (W)
#' @export
metabolic_rate <- function(volume, rho, Q10, a0, b, Tref, Tbody) {
  m <- volume * rho
  a0 * m^b * Q10^((Tbody - Tref) / 10)
}

#' Calculates ectotherm body temperature over a single timestep
#'
#' @param ecto_params A list of ectotherm parameters as returned by
#' [micropoint::create_ectop()]
#' @param tme POSIXlt object indicating date and time in UTC
#' @param microcdata A list of microclimate variables with the following entries:
#'   \describe{
#'     \item{Tair}{Air temperature (degrees C).}
#'     \item{Tsurf}{Surface temperature (degrees C).}
#'     \item{Rdirdown}{Direct shortwave radiation flux perpendicular to the solar beam (W m^-2).}
#'     \item{Rdifdown}{Diffuse downward shortwave radiation flux (W m^-2).}
#'     \item{Rswup}{Upward shortwave radiation flux (W m^-2).}
#'     \item{Rlwdown}{Downward longwave radiation flux (W m^-2).}
#'     \item{Rlwup}{Upward longwave radiation flux (W m^-2).}
#'     \item{winddir}{Wind direction (degrees from north).}
#'     \item{windspeed}{Wind speed at animal height (m s^-1).}
#'     \item{relhum}{Air relative humidity (Percentage).}
#'     \item{pk}{Atmospheric pressure (kPa).}
#'     \item{surfrelhum}{Relative humidity at the surface (Percentage).}
#'   }
#' @param lat Latitude (decimal degrees).
#' @param long Longitude (decimal degrees).
#' @param tolerance Convergence tolerance used when iterating (degrees C).
#' @param maxiter Maximum number of iterations.
#'
#' @return Body temperature (deg C)
#' @export
Ectothermtemp <- function(ecto_params, tme, microcdata, lat, long, tolerance = 1e-3, maxiter = 100) {

  # ------------------- Set constants ---------------------------------------- #
  RgasC <- 8.314 # Universal gas constant
  sb <- 5.67e-8 # Stefan-Boltzman constant
  D0 <- 2.26e-5 # water diiffusivity and standard pressure and temperature
  P0 <- 101.325 # standard pressure
  T0 <- 273.15 # standard temperature
  # ---------------- Calculate absorbed radiation ---------------------------- #
  # Calculate solar position
  solp <- sunposition(tme, lat, long)
  # Calculate absorbed radiation
  sil <- with(ecto_params, silhouette(solp$zen, solp$azi, height, width, len,
                                      adir, atilt, position))
  sc <-  with(sil, silA / SurfArea)
  if (ecto_params$confrac < 0.5) { # Only lower portion in contact with surface
    Rswabs <- (sc * microcdata$Rdirdown + 0.5 * microcdata$Rdifdown +
                 0.5 * (1 - ecto_params$confrac)  * microcdata$Rswup) *
      (1 - ecto_params$refl)
    Rlwabs <- 0.5 * (microcdata$Rlwdown + (1 - ecto_params$confrac)  *
                       microcdata$Rlwup) * ecto_params$em
  } else {  # Upper portion in contact with surface too
    upcon <-  ecto_params$confrac - 0.5 # upper portion in contact with surface
    Rswabs <- (sc * microcdata$Rdirdown + 0.5 * microcdata$Rdifdown) *
      (1 - ecto_params$refl) * (1 - upcon)
    Rlwabs <- (0.5 * microcdata$Rlwdown) * (1 - upcon) * ecto_params$em
  }
  Rabs <- Rswabs + Rlwabs
  # ---------------- Calculate Characteristic dimension (m) ------------------ #
  svals <- with(ecto_params, silhouette(zen = 90, azi = microcdata$winddir,
                                        height, width, len, adir, atilt, position))
  d <- svals$Volume / svals$silA
  # ------------------------------- Initialize values ------------------------ #
  error <- 1e99
  iter <- 1
  Tbody <- microcdata$Tair + 2
  while (error > tolerance && iter < maxiter) {
    Te <- (Tbody + microcdata$Tair) / 2
    Tf <- (Tbody + microcdata$Tsurf) / 2
    # -------------- Compute boundary layer resistance ----------------------- #
    dT <- abs(Tbody - microcdata$Tair)
    rHa <- animalrHa(microcdata$Tair, dT, microcdata$windspeed, d)
    # ------------------------- Compute metaboloc rate ----------------------- #
    M <- with(ecto_params, metabolic_rate(volume, rho, Q10, a0, b, Tref, Tbody)
              / area)
    # ----------------------- Re-compute body tmeperature -------------------- #
    cp <- cpair(microcdata$Tair)
    ph <- rhohair(microcdata$Tair, microcdata$pk)
    cd <- 1 - ecto_params$confrac
    ea <- satvap(microcdata$Tair)
    Da <- ea * (1 - microcdata$relhum / 100) # Vapour pressure deficit
    Dac <- (microcdata$surfrelhum / 100) * satvap(Tf) - ea #  contact Vapour pressure deficit
    # latent heat exchange through a contact interface
    keff <- with(ecto_params, k / (0.5 * height))
    # Calculate mus
    mr <- cd * ecto_params$em * 5.67e-8 # emmited ratiation mu
    mrc <- with(ecto_params, confrac * em * 5.67e-8) # emmited ratiation mu
    mc <- keff * ecto_params$confrac # conductance mu
    mv <- ph * cp * cd / rHa
    la <- latvap(Tbody)
    # Total vapour resistance for exposed surface
    rv <- rHa + ecto_params$rc
    ml <- (ph * la * cd) / (rv * microcdata$pk)
    Dw <- D0 * (P0 / microcdata$pk) *
      ((Tf + 273.15) / T0)^1.75 # water diffusivity
    rv_contact <- ecto_params$rc + ((0.5 * ecto_params$height *
                         (microcdata$Tsurf + 273.15) * RgasC) / (Dw * 1000))
    mlc <- (la * ecto_params$confrac) / rv_contact
    DeV <- satvap(Te + 0.5) - satvap(Te - 0.5) # Slope of the saturated vapour pressure curve (Lv)
    DeC <- satvap(Tf + 0.5) - satvap(Tf - 0.5) # Slope of the saturated vapour pressure curve (Lc)
    # Calculate various contributions to body temperature
    top <- (Rabs + M - mr * (microcdata$Tair + 273.15)^4 -
              mrc * (microcdata$Tair + 273.15)^4 -
              mc * dT - ml * Da - mlc * Dac)
    btm <- 4 * mr * (Te + 273.15)^3 + 4 * mrc * (Tf + 273.15)^3 +
      mc + mv + ml * DeV + (microcdata$surfrelhum / 100) * mlc * DeC
    newTbody <- microcdata$Tair + top / btm
    error <- abs(newTbody - Tbody)
    error
    Tbody <- newTbody
    iter <- iter + 1
  }
  return (Tbody)
}


#' Mammal heat transfer resistance
#'
#' Calculates total resistance to heat loss for a mammal, combining boundary
#' layer, coat, and vascular resistances under warm and cold conditions.
#'
#' @param height Body height (cm)
#' @param width Body width (cm)
#' @param len Body length (cm)
#' @param coatthick Coat thickness (cm)
#' @param tbody Body temperature (degrees C)
#' @param tair Air temperature (degrees C)
#' @param adir body azimuth direction (degrees)
#' @param atilt body tilt from vertical (degrees)
#' @param position assumption about body orientation. See [chardim()]
#' @param windspeed Wind speed at animal height (m s^-1)
#' @param wdir wind direction (degrees)
#' @return A list with:
#'   \describe{
#'     \item{rHa_bound}{boundary layer resistance (s m^-1)}
#'     \item{rHa_warm}{Total resistance (s m^-1) under vasodilated (warm) conditions}
#'     \item{rHa_cold}{Total resistance (s m^-1) under vasoconstricted (cold) conditions}
#'   }
#'
#' @details
#' Total resistance is calculated as the sum of boundary layer resistance,
#' coat resistance, and a fixed vascular resistance representing either
#' vasodilation or vasoconstriction. Coat conductance decreases with coat
#' thickness and increases with wind speed.
#'
#' @export
mammal_resistance <- function(height, width, len, coatthick, tbody, tair,
                              adir = 0.0, atilt = 0.0, position = "fixed",
                              windspeed = 0, wdir = 0) {
  # boundary layer resistance
  d <- chardim(wdir, height, width, len, adir, atilt, position)
  dT <- abs(tbody - tair)
  rHa <- animalrHa(tair, dT, windspeed, d)
  # coat resistance
  kmm0 <- 2.7058 * coatthick^(-0.803)
  kmm <-  kmm0 * (1 + 0.1 * windspeed)
  rHcoat <- 1000 / kmm
  # vascular conductance
  rHvas_cold <- 110
  rHvas_warm <- 30
  return(list(rHa_bound = rHa,
              rHa_warm = rHa + rHcoat + rHvas_warm,
              rHa_cold = rHa + rHcoat + rHvas_cold))
}

#' Bird heat transfer resistance
#'
#' Calculates total resistance to heat loss for a bird, combining boundary
#' layer resistance with plumage conductance under minimum and maximum
#' insulation states.
#'
#' @param height Body height (cm)
#' @param width Body width (cm)
#' @param len Body length (cm)
#' @param rho Body density (kg m^-3). Default is 730.
#' @param tbody Body temperature (degrees C)
#' @param tair Air temperature (degrees C)
#' @param adir body azimuth direction (degrees)
#' @param atilt body tilt from vertical (degrees)
#' @param position assumption about body orientation. See [chardim()]
#' @param windspeed Wind speed at animal height (m s^-1)
#' @param wdir wind direction (degrees)
#'
#' @return A list with:
#'   \describe{
#'     \item{rHa_bound}{boundary layer resistance (s m^-1)}
#'     \item{rHa_warm}{Total resistance (s m^-1) under reduced insulation (warm conditions)}
#'     \item{rHa_cold}{Total resistance (s m^-1) under increased insulation (cold conditions)}
#'   }
#'
#' @details
#' Boundary layer resistance is calculated from body geometry and wind speed.
#' Plumage conductance is estimated from body mass using an allometric
#' relationship for minimum conductance, with a higher conductance representing
#' reduced insulation. Conductance is converted to resistance and adjusted for
#' wind speed before being combined with boundary layer resistance.
#'
#' @export
bird_resistance <- function(height, width, len, rho = 730, tbody, tair,
                            adir, atilt, position, windspeed, wdir) {
  # boundary layer resistance
  d <- chardim(wdir, height, width, len, adir, atilt, position)
  dT <- abs(tbody - tair)
  rHa <- animalrHa(20, 5, windspeed, d)
  # calculate body mass
  vol <- (4 / 3) * pi * (height / 200) * (width / 200) * (len / 200)
  mass <- rho * vol
  g_min <- 0.06 * mass^(-0.15)
  g_max <- g_min * 1.75
  k0_min <- g_min * 8.314 * ((tair + 273.15) / 101325)
  k0_max <- g_max * 8.314 * ((tair + 273.15) / 101325)
  k_min <- k0_min * (1 + 0.05 * windspeed)
  k_max <- k0_max * (1 + 0.05 * windspeed)
  return(list(rHa_bound = rHa,
              rHa_warm = rHa + 1/ k_min,
              rHa_cold = rHa + 1/ k_max))
}

#' Endotherm energy balance
#'
#' Calculates net energy balance for an endotherm at a single time step from
#' body traits and local microclimate.
#'
#' @param endo_params A list of endotherm parameters with the following entries:
#'   \describe{
#'     \item{height}{Body height (cm).}
#'     \item{width}{Body width (cm).}
#'     \item{len}{Body length (cm).}
#'     \item{refl}{Shortwave reflectance of the body surface (unitless).}
#'     \item{em}{Longwave emissivity of the body surface (unitless).}
#'     \item{confrac}{Fraction of body surface in conductive contact with the substrate (unitless).}
#'     \item{adir}{Body orientation relative to north (degrees).}
#'     \item{atilt}{Body tilt angle relative to the horizontal (degrees).}
#'     \item{position}{Descriptor of body position (character).}
#'     \item{Tbody}{Body temperature (degrees C).}
#'     \item{type}{Endotherm type; one of \code{"mammal"} or \code{"bird"}.}
#'     \item{basmet}{Basal metabolic rate (W).}
#'     \item{coatthick}{Coat thickness (cm; mammals only).}
#'     \item{rho}{Body density (kg m^-3; birds only).}
#'   }
#' @param tme POSIXlt object indicating date and time in UTC.
#' @param microcdata A list of microclimate variables with the following entries:
#'   \describe{
#'     \item{Tair}{Air temperature (degrees C).}
#'     \item{Tsurf}{Surface temperature (degrees C).}
#'     \item{Rdirdown}{Direct shortwave radiation flux perpendicular to the solar beam (W m^-2).}
#'     \item{Rdifdown}{Diffuse downward shortwave radiation flux (W m^-2).}
#'     \item{Rswup}{Upward shortwave radiation flux (W m^-2).}
#'     \item{Rlwdown}{Downward longwave radiation flux (W m^-2).}
#'     \item{Rlwup}{Upward longwave radiation flux (W m^-2).}
#'     \item{winddir}{Wind direction (degrees from north).}
#'     \item{windspeed}{Wind speed at animal height (m s^-1).}
#'     \item{pk}{Atmospheric pressure (kPa).}
#'   }
#' @param lat Latitude (decimal degrees).
#' @param long Longitude (decimal degrees).
#'
#' @return Net energy balance (W), with positive values indicating heat gain
#'   and negative values indicating heat loss.
#'
#' @details
#' Absorbed shortwave and longwave radiation are calculated from body geometry,
#' orientation and surface properties. Emitted longwave radiation is determined
#' from body temperature. Sensible heat exchange is calculated using resistances
#' appropriate for mammals or birds, representing different insulation states.
#' Basal metabolic heat production is added and the resulting energy balance is
#' returned for the whole animal.
#'
#' @export
Endothermstress <- function(endo_params, tme, microcdata, lat, long) {
  # ------------------- Set constants ---------------------------------------- #
  sb <- 5.67e-8 # Stefan-Boltzman constant
  # ---------------- Calculate absorbed radiation ---------------------------- #
  # Calculate solar position
  solp <- sunposition(tme, lat, long)
  # Calculate absorbed radiation
  sil <- with(endo_params, silhouette(solp$zen, solp$azi, height, width, len,
                                      adir, atilt, position))
  sc <-  with(sil, silA / SurfArea)
  if (endo_params$confrac < 0.5) { # Only lower portion in contact with surface
    Rswabs <- (sc * microcdata$Rdirdown + 0.5 * microcdata$Rdifdown +
                 0.5 * (1 - endo_params$confrac)  * microcdata$Rswup) *
      (1 - endo_params$refl)
    Rlwabs <- 0.5 * (microcdata$Rlwdown + (1 - endo_params$confrac)  *
                       microcdata$Rlwup) * endo_params$em
  } else {  # Upper portion in contact with surface too
    upcon <-  endo_params$confrac - 0.5 # upper portion in contact with surface
    Rswabs <- (sc * microcdata$Rdirdown + 0.5 * microcdata$Rdifdown) *
      (1 - endo_params$refl) * (1 - upcon)
    Rlwabs <- (0.5 * microcdata$Rlwdown) * (1 - upcon) * endo_params$em
  }
  Rabs <- Rswabs + Rlwabs
  Rabs <- Rswabs + Rlwabs
  # ----------------- Calculate emitted radiation ---------------------------- #
  Rem <- with(endo_params, em * sb * (Tbody + 273.15)^4)
  # ----------------- Calculate sensible heat flux --------------------------- #
  tair <- microcdata$Tair
  windspeed <- microcdata$windspeed
  wdir <-   microcdata$winddir
  if (endo_params$type == "mammal") {
    rHas <- with(endo_params, mammal_resistance(height, width, len, coatthick,
            Tbody, tair, adir, atilt, position, windspeed, wdir))
  } else if (endo_params$type == "bird") {
    rHas <- with(endo_params, bird_resistance(height, width, len, rho,
            Tbody, tair, adir, atilt, position, windspeed, wdir))
  } else stop("Endotherm type not recognised")
  ph <- rhohair(tair, microcdata$pk)
  cp <- cpair(tair)
  H_warm <- (ph * cp / rHas$rHa_warm) * (endo_params$Tbody - tair)
  H_cold <- (ph * cp / rHas$rHa_cold) * (endo_params$Tbody - tair)
  # ----------------- Calculate contact heat flux --------------------------- #
  Hc_warm <- (ph * cp / (rHas$rHa_warm - rHas$rHa_bound)) *
    (endo_params$Tbody - microcdata$Tsurf)
  Hc_cold <- (ph * cp / (rHas$rHa_cold - rHas$rHa_bound)) *
    (endo_params$Tbody - microcdata$Tsurf)
  HH_warm <- (1 - endo_params$confrac) * H_warm + endo_params$confrac * Hc_warm
  HH_cold <- (1 - endo_params$confrac) * H_cold + endo_params$confrac * Hc_cold
  M <- endo_params$basmet / sil$SurfArea
  # --------------------- Calculate energy balance --------------------------- #
  EB_cold <- Rabs - Rem - HH_cold + M # W/m^2
  EB_warm <- Rabs - Rem - HH_warm + M # W/m^2
  if (EB_warm > 0) {
    EB <- EB_warm
  } else {
    EB <- EB_cold
  }
  EB <- EB * sil$SurfArea
  return(EB)
}
