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
#' @param atilt body tilt (degrees): the angle by which the long axis of the body
#' is raised from horizontal at its front, which tilts its vertical axis from
#' vertical by the same angle. Only used when \code{position = "fixed"} or
#' \code{position = "randomdir"}. Default is 0.
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
#' solution exists for a general triaxial ellipsoid. With random orientation,
#' silhouette area is a quarter of surface area, as it is for any convex body.
#' @export
silhouette <- function(zen, azi, height, width, len, adir = 0, atilt = 0.0, position = "fixed") {
  # Semi-axes of the ellipsoid (m): along the length, across and up
  AA <- len / 200
  BB <- width / 200
  CC <- height / 200

  # Surface area, from Knud Thomsen's approximation, and volume
  P <- 1.6075
  Area <- 4 * pi * (((AA^P * BB^P) + (AA^P * CC^P) + (BB^P * CC^P)) / 3)^(1 / P)
  Volume <- (4 / 3) * pi * AA * BB * CC

  # Area of the shadow the body casts on a plane facing the sun, for a body
  # heading towards adir with its front raised by atilt. The sun's direction is
  # resolved along the length, width and height of the body, and each component
  # is weighted by the cross-section of the ellipsoid seen from that axis
  shadow <- function(adir) {
    z <- zen * pi / 180
    theta <- (adir - azi) * pi / 180
    t <- atilt * pi / 180
    L <- sin(z) * cos(t) * cos(theta) + cos(z) * sin(t)
    M <- sin(z) * sin(theta)
    N <- cos(z) * cos(t) - sin(z) * sin(t) * cos(theta)
    pi * sqrt(L^2 * BB^2 * CC^2 + M^2 * CC^2 * AA^2 + N^2 * AA^2 * BB^2)
  }

  if (position == "max") { # body side-on to the sun
    silA <- pi * max(c(BB * CC, AA * CC, AA * BB))
  } else if (position == "min") { # body end-on to the sun
    silA <- pi * min(c(BB * CC, AA * CC, AA * BB))
  } else if (position == "fixed") {
    silA <- shadow(adir)
  } else if (position == "randomdir") { # shadow averaged over all headings
    s <- sapply(0:359, shadow)
    silA <- rowMeans(matrix(s, ncol = 360))
  } else if (position == "random") { # shadow averaged over all orientations
    silA <- Area / 4
  } else stop("Position not recognised")
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
#' @param atilt body tilt (degrees): the angle by which the long axis of the body
#' is raised from horizontal at its front. Only used when
#' \code{position = "fixed"} or \code{position = "randomdir"}. Default is 0.
#' @param position assumption about body orientation. Must be one of:
#' \describe{
#'   \item{\code{"fixed"}}{body orientation is fixed according to
#'   \code{adir} and \code{atilt}}
#'   \item{\code{"max"}}{body turned to present its largest cross-section to
#'   the sun}
#'   \item{\code{"min"}}{body turned to present its smallest cross-section to
#'   the sun}
#'   \item{\code{"randomdir"}}{body azimuth direction is assumed random, but
#'   tilt is fixed by \code{atilt}}
#'   \item{\code{"random"}}{both body azimuth direction and tilt are assumed
#'   random}
#' }
#' @param zen solar zenith angle (degrees). Only used, and then required, when
#' \code{position} is \code{"max"} or \code{"min"}.
#' @param azi solar azimuth angle (degrees). Only used, and then required, when
#' \code{position} is \code{"max"} or \code{"min"}.
#' @return a numeric value giving the characteristic dimension (m), defined as
#' volume divided by projected area
#' @details The characteristic dimension is calculated as:
#' \deqn{d = V / A_p}
#' where \eqn{V} is ellipsoid volume and \eqn{A_p} is the area the body presents
#' to a horizontal wind from direction \code{wdir}. For an ellipsoid this is
#' two-thirds of the length of the body along the wind. When \code{position} is
#' \code{"max"} or \code{"min"}, the body is first turned towards the sun, and
#' the area it then presents to the wind is used. Turning is by heading and
#' tilt only, so when the width is the axis that should point at the sun, the
#' body lies across the sun's direction with its long axis upright.
#' @export
chardim <- function(wdir, height, width, len,
                    adir = 0.0, atilt = 0.0, position = "fixed", zen = NULL, azi = NULL) {
  # A body that maximises or minimises the sunlight it intercepts points its
  # shortest or longest axis at the sun. Its orientation is set that way, and
  # its size along the wind follows from that orientation
  if (position %in% c("max", "min")) {
    if (is.null(zen) || is.null(azi)) stop("zen and azi are needed when position is max or min")
    dims <- c(len, width, height)
    axis <- if (position == "max") which.min(dims) else which.max(dims)
    if (axis == 1) { # long axis towards the sun
      adir <- azi
      atilt <- 90 - zen
    } else if (axis == 2) { # side towards the sun's azimuth, long axis upright
      adir <- azi + 90
      atilt <- 90
    } else { # back towards the sun
      adir <- azi
      atilt <- -zen
    }
    position <- "fixed"
  }
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
#' (degrees C). Its sign does not matter.
#' @param uz wind speed at the animal scale (m s^-1)
#' @param d characteristic dimension of the animal (m): volume divided by the
#' area it presents to the wind, as returned by [chardim()]
#' @param rHmax maximum allowed boundary layer resistance (s m^-1). Default is
#' 300.
#' @return a numeric value giving the boundary layer resistance to heat transfer
#' (s m^-1)
#' @details Thermal diffusivity and kinematic viscosity are calculated as
#' functions of air temperature. The animal is treated as a spheroid, the size
#' of which is the length of its body along the wind (1.5 times `d`). The Nusselt
#' numbers are those for a sphere: Whitaker (1972) for forced convection and
#' Churchill (1983) for free convection. The two contributions are combined
#' using a cubic mean to represent mixed convection.
#' @references Churchill SW (1983) Free convection around immersed bodies. In:
#' Schlunder EU (ed) Heat Exchanger Design Handbook, Section 2.5.7. Hemisphere,
#' New York.
#'
#' Whitaker S (1972) Forced convection heat transfer correlations for
#' flow in pipes, past flat plates, single cylinders, single spheres, and for flow
#' in packed beds and tube bundles. AIChE Journal 18: 361-371.
#' @export
animalrHa <- function(tair, dT, uz, d, rHmax = 300.0) {
  Tk <- tair + 273.15
  # Properties of air that set how readily it carries heat and how readily it
  # flows: thermal diffusivity and kinematic viscosity (m^2/s)
  Kh <- (1.6667e-10 * Tk^2 + 2.9935e-8 * Tk - 1.7128e-6)
  v <- 1.326e-5 * (Tk / 273.15)^2.5 * (393.15 / (Tk + 120.0))
  # Size of the animal as the wind meets it (m): the length of its body along
  # the wind. For a spheroid this is 1.5 times its volume over the area it
  # presents to the wind, and for a sphere it is the diameter
  D <- 1.5 * d

  # Heat carried away by the wind (forced convection). The Reynolds number
  # measures the strength of the flow past the animal, and the Nusselt number how
  # much faster heat is lost than through still air
  Re <- (uz * D) / v
  Pr <- v / Kh # Prandtl number: how readily air flows relative to how readily it carries heat
  # Heat is lost by conduction into still air (the 2), to smooth flow over the
  # front of the body (rising with the square root of the Reynolds number) and
  # to the wake behind it (rising faster)
  Nuf <- 2 + (0.4 * sqrt(Re) + 0.06 * Re^(2/3)) * Pr^0.4

  # Heat carried away by air that rises or sinks because the animal warms or
  # cools it (free convection), which matters when there is little wind. The
  # Grashof number measures the strength of that buoyancy. As for forced
  # convection, the animal is treated as a sphere, which in still air loses heat
  # by conduction at least
  Gr <- (9.81 * D^3 * abs(dT)) / (Tk * v^2)
  Nun <- 2 + 0.589 * (Gr * Pr)^(1/4) / (1 + (0.469 / Pr)^(9/16))^(4/9)

  # The two act together, the stronger dominating
  Nu <- (Nuf^3 + Nun^3)^(1/3)
  # Resistance of the layer of air around the animal to the loss of heat
  rHa <- pmin(D / (Kh * Nu), rHmax)
  return(rHa)
}

#' Calculates metabolic heat production as a function of body mass and temperature.
#'
#' @param volume Animal volume (m^3)
#' @param rho Animal density (kg m^-3)
#' @param Q10 Temperature sensitivity coefficient
#' @param a0 Normalisation constant (W kg^-b)
#' @param b Mass-scaling exponent
#' @param Tref Reference temperature (°C)
#' @param Tbody Body temperature (°C)
#'
#' @return Numeric vector of metabolic rate (W)
#' @export
metabolic_rate <- function(volume, rho, Q10, a0, b, Tref, Tbody) {
  # Metabolic rate rises less than in proportion to body mass, and by a factor
  # of Q10 for every 10 degrees the body is warmer than the reference temperature
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
  sb <- 5.67e-8 # Stefan-Boltzmann constant (W m^-2 K^-4)

  # Radiation absorbed per unit body surface: the direct beam on the silhouette
  # the body presents to the sun, and half the diffuse, reflected and longwave
  # radiation from above and below, since half the body faces each way
  solp <- sunposition(tme, lat, long)
  sil <- with(ecto_params, silhouette(solp$zen, solp$azi, height, width, len,
                                      adir, atilt, position))
  sc <- sil$silA / sil$SurfArea
  Rsw <- sc * microcdata$Rdirdown + 0.5 * microcdata$Rdifdown + 0.5 * microcdata$Rswup
  Rabs <- (1 - ecto_params$refl) * Rsw +
    0.5 * ecto_params$em * (microcdata$Rlwdown + microcdata$Rlwup)

  # Size of the body as the wind meets it, which sets its boundary layer
  d <- with(ecto_params, chardim(microcdata$winddir, height, width, len, adir, atilt,
                                 position, solp$zen, solp$azi))

  # Heat exchange that does not depend on body temperature. The share of the body
  # touching the ground exchanges heat with it by conduction, and water vapour by
  # diffusion, across half the height of the body
  Tair <- microcdata$Tair
  Tsurf <- microcdata$Tsurf
  pk <- microcdata$pk
  cf <- ecto_params$confrac
  cd <- 1 - cf
  cp <- cpair(Tair)
  ph <- rhoair(Tair, pk)
  halfheight <- 0.5 * ecto_params$height / 100 # m
  mc <- cf * ecto_params$k / halfheight # conductance to the ground (W m^-2 K^-1)
  mr <- cd * ecto_params$em * sb # longwave emission from the exposed share
  mrc <- cf * ecto_params$em * sb # longwave emission from the share touching the ground
  # Moist skin is saturated at body temperature. Vapour pressure differences are
  # taken at air temperature, with the change at body temperature added below
  es <- satvap(Tair)
  Da <- es - es * microcdata$relhum / 100 # skin against the air
  Dac <- es - (microcdata$surfrelhum / 100) * satvap(Tsurf) # skin against the ground

  # Body temperature from the energy balance, linearised about air temperature.
  # Boundary layer resistance, metabolic rate and the linearisation depend on
  # body temperature, so it is found by iteration
  error <- 1e99
  iter <- 1
  Tbody <- Tair + 2
  while (error > tolerance && iter < maxiter) {
    Te <- (Tbody + Tair) / 2
    Tf <- (Tbody + Tsurf) / 2
    rHa <- animalrHa(Tair, Tbody - Tair, microcdata$windspeed, d)
    M <- with(ecto_params, metabolic_rate(volume, rho, Q10, a0, b, Tref, Tbody) / area)
    # Sensible and latent heat from the exposed share, with vapour leaving the
    # skin through its own resistance and that of the boundary layer
    la <- latvap(Tbody)
    mv <- ph * cp * cd / rHa
    ml <- (ph * la * cd) / ((rHa + ecto_params$rc) * pk)
    # Latent heat through the share touching the ground, by diffusion through the
    # air in the pores between body and ground
    Dw <- 2.26e-5 * (101.325 / pk) * ((Tf + 273.15) / 273.15)^1.75 # diffusivity of water vapour (m^2/s)
    mlc <- (ph * la * cf) / ((halfheight / Dw + ecto_params$rc) * pk)
    DeV <- satvap(Te + 0.5) - satvap(Te - 0.5) # slope of the saturated vapour pressure curve
    # Heat gained less heat lost at air temperature, over the rate at which losses
    # rise as the body warms
    top <- Rabs + M - (mr + mrc) * (Tair + 273.15)^4 + mc * (Tsurf - Tair) - ml * Da - mlc * Dac
    btm <- 4 * mr * (Te + 273.15)^3 + 4 * mrc * (Tf + 273.15)^3 + mc + mv + (ml + mlc) * DeV
    newTbody <- Tair + top / btm
    error <- abs(newTbody - Tbody)
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
  sb <- 5.67e-8 # Stefan-Boltzmann constant
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
  ph <- rhoair(tair, microcdata$pk)
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
#' Physiologically equivalent temperature for an endotherm
#'
#' @description Computes the physiologically equivalent temperature (PET) for
#' an endotherm at a single time step: the air temperature in a standard
#' reference environment (calm wind, no additional radiation) that would
#' produce the same body temperature and heat loss as the animal has in the
#' environment being considered.
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
#' @param uref  reference (standardising) wind speed used to express the
#' @param iters integer specifying how many times to iterate the equation (4 achieves near exact solution)
#'
#' @return Physiologically equivalent temperature (degrees C)
#' @export
EndothermPET <- function(endo_params, tme, microcdata, lat, long, uref = 0.1, iters = 4) {
  sb <- 5.67e-8
  # ---------------- Absorbed and emitted radiation (as Endothermstress) ----- #
  solp <- sunposition(tme, lat, long)
  sil <- with(endo_params, silhouette(solp$zen, solp$azi, height, width, len,
                                     adir, atilt, position))
  sc <- with(sil, silA / SurfArea)
  if (endo_params$confrac < 0.5) {
    Rswabs <- (sc * microcdata$Rdirdown + 0.5 * microcdata$Rdifdown +
                 0.5 * (1 - endo_params$confrac) * microcdata$Rswup) *
      (1 - endo_params$refl)
    Rlwabs <- 0.5 * (microcdata$Rlwdown + (1 - endo_params$confrac) *
                       microcdata$Rlwup) * endo_params$em
  } else {
    upcon <- endo_params$confrac - 0.5
    Rswabs <- (sc * microcdata$Rdirdown + 0.5 * microcdata$Rdifdown) *
      (1 - endo_params$refl) * (1 - upcon)
    Rlwabs <- (0.5 * microcdata$Rlwdown) * (1 - upcon) * endo_params$em
  }
  Rabs <- Rswabs + Rlwabs
  Rem <- with(endo_params, em * sb * (Tbody + 273.15)^4)
  # ---------------- Actual sensible + contact heat exchange ------------------ #
  tair <- microcdata$Tair
  windspeed <- microcdata$windspeed
  wdir <- microcdata$winddir
  if (endo_params$type == "mammal") {
    rHas <- with(endo_params, mammal_resistance(height, width, len, coatthick,
                                                Tbody, tair, adir, atilt, position, windspeed, wdir))
    rHas_ref <- with(endo_params, mammal_resistance(height, width, len, coatthick,
                                                    Tbody, tair, adir, atilt, position, uref, wdir))
  } else if (endo_params$type == "bird") {
    rHas <- with(endo_params, bird_resistance(height, width, len, rho,
                                              Tbody, tair, adir, atilt, position, windspeed, wdir))
    rHas_ref <- with(endo_params, bird_resistance(height, width, len, rho,
                                                  Tbody, tair, adir, atilt, position, uref, wdir))
  } else stop("Endotherm type not recognised")
  ph <- rhoair(tair, microcdata$pk)
  cp <- cpair(tair)
  H_warm <- (ph * cp / rHas$rHa_warm) * (endo_params$Tbody - tair)
  H_cold <- (ph * cp / rHas$rHa_cold) * (endo_params$Tbody - tair)
  Hc_warm <- (ph * cp / (rHas$rHa_warm - rHas$rHa_bound)) *
    (endo_params$Tbody - microcdata$Tsurf)
  Hc_cold <- (ph * cp / (rHas$rHa_cold - rHas$rHa_bound)) *
    (endo_params$Tbody - microcdata$Tsurf)
  HH_warm <- (1 - endo_params$confrac) * H_warm + endo_params$confrac * Hc_warm
  HH_cold <- (1 - endo_params$confrac) * H_cold + endo_params$confrac * Hc_cold
  # ---------------- Pick regulatory state exactly as Endothermstress -------- #
  M <- endo_params$basmet / sil$SurfArea
  EB_warm <- Rabs - Rem - HH_warm + M
  if (EB_warm > 0) {
    Renv <- Rabs - Rem - HH_warm
    rHa_ref <- rHas_ref$rHa_warm
  } else {
    Renv <- Rabs - Rem - HH_cold
    rHa_ref <- rHas_ref$rHa_cold
  }
  # ---------- Solve for equivalent temperature (Penman-style linearisation) - #
  # Reference environment = blackbody enclosure at Te with convective exchange
  # at the reference wind speed; both must together reproduce Renv:
  #   em*sb*((Te+273.15)^4-(Tbody+273.15)^4) + (ph*cp/rHa_ref)*(Te-Tbody) = Renv
  # Solved by linearising the radiative term and iterating (cf. the "iters"
  # argument on SolveEnergyBalance), rather than a general root solver.
  gconv <- (ph * cp) / rHa_ref
  Te <- endo_params$Tbody
  for (i in 1:iters) {
    grad <- 4 * endo_params$em * sb * (Te + 273.15)^3
    resid <- endo_params$em * sb * ((Te + 273.15)^4 - (endo_params$Tbody + 273.15)^4) +
      gconv * (Te - endo_params$Tbody) - Renv
    Te <- Te - resid / (grad + gconv)
  }
  return(Te)
}
