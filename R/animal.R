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

  # Radiation absorbed per unit body surface
  solp <- sunposition(tme, lat, long)
  Rabs <- animalRabs(ecto_params, solp, microcdata)

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

#' @title Radiation absorbed by an animal
#' @description Calculates the shortwave and longwave radiation absorbed per unit
#' body surface by an animal treated as an ellipsoid.
#' @param params a list of animal parameters containing \code{height},
#' \code{width}, \code{len} (cm), \code{adir}, \code{atilt} (degrees),
#' \code{position} (see [silhouette()]), \code{refl} (shortwave reflectance) and
#' \code{em} (longwave emissivity)
#' @param solp a list with the solar zenith (\code{zen}) and azimuth (\code{azi})
#' angles (degrees), as returned by [sunposition()]
#' @param microcdata a list of radiation fluxes at the animal (W m^-2):
#' \code{Rdirdown} (direct, perpendicular to the solar beam), \code{Rdifdown}
#' (diffuse downward), \code{Rswup} (upward shortwave), \code{Rlwdown} and
#' \code{Rlwup} (downward and upward longwave)
#' @return radiation absorbed per unit body surface (W m^-2)
#' @details The direct beam falls on the silhouette the body presents to the sun.
#' Half the body faces upward and half downward, so it receives half the diffuse
#' and longwave radiation from the sky and half that reflected or emitted from below.
#' @export
animalRabs <- function(params, solp, microcdata) {
  # Share of the body surface that the direct beam strikes
  sil <- with(params, silhouette(solp$zen, solp$azi, height, width, len, adir, atilt, position))
  sc <- sil$silA / sil$SurfArea
  # Shortwave absorbed from the direct beam and from diffuse and reflected light,
  # and longwave from the sky and the ground
  Rsw <- sc * microcdata$Rdirdown + 0.5 * microcdata$Rdifdown + 0.5 * microcdata$Rswup
  Rabs <- (1 - params$refl) * Rsw + 0.5 * params$em * (microcdata$Rlwdown + microcdata$Rlwup)
  return(Rabs)
}
#' @title Operative temperature
#' @description Calculates the operative temperature of an animal surface: the
#' temperature of uniform surroundings with which the surface would exchange the
#' same heat by radiation and convection as in the actual environment.
#' @param Rabs radiation absorbed per unit surface (W m^-2), as returned by
#' [animalRabs()]
#' @param tair air temperature (degrees C)
#' @param rHa boundary layer resistance to heat (s m^-1), as returned by [animalrHa()]
#' @param em longwave emissivity of the surface
#' @param pk atmospheric pressure (kPa)
#' @return a list with:
#' \describe{
#'   \item{Te}{operative temperature (degrees C)}
#'   \item{re}{resistance of the surroundings to heat (s m^-1): boundary layer
#'   and radiative resistance in parallel}
#' }
#' @details Longwave emission is linearised about air temperature, so radiation
#' acts as a resistance in parallel with the boundary layer. Operative temperature
#' equals air temperature when the surface absorbs as much radiation as it would
#' emit at air temperature, as in shade with surroundings at air temperature.
#' @export
operativetemp <- function(Rabs, tair, rHa, em, pk) {
  sb <- 5.67e-8
  rcp <- rhoair(tair, pk) * cpair(tair) # volumetric heat capacity of air (J m^-3 K^-1)
  # Radiative exchange behaves as a resistance alongside convection, smaller the
  # warmer the surroundings
  rR <- rcp / (4 * em * sb * (tair + 273.15)^3)
  re <- 1 / (1 / rHa + 1 / rR)
  # Air temperature raised by the radiation absorbed beyond what the surface
  # would emit at air temperature
  Te <- tair + re * (Rabs - em * sb * (tair + 273.15)^4) / rcp
  return(list(Te = Te, re = re))
}
#' @title Resistance of a coat to heat
#' @description Calculates the resistance of a layer of fur or plumage to the flow
#' of heat from skin to coat surface.
#' @param depth coat depth (cm)
#' @param kcoat effective thermal conductivity of the coat (W m^-1 K^-1),
#' including radiation between fibres. Typically 0.04 to 0.08.
#' @param tair air temperature (degrees C)
#' @param pk atmospheric pressure (kPa)
#' @return resistance to heat (s m^-1)
#' @export
coatresistance <- function(depth, kcoat, tair, pk) {
  # Conduction across the coat, expressed in the same units as a boundary layer
  rcp <- rhoair(tair, pk) * cpair(tair)
  rc <- rcp * (depth / 100) / kcoat
  return(rc)
}
#' @title Resistance of body tissue to heat
#' @description Calculates the resistance to the flow of heat from the body core
#' to the skin, for a body in which heat is produced throughout.
#' @param volume body volume (m^3)
#' @param area body surface area (m^2)
#' @param kflesh effective thermal conductivity of the flesh (W m^-1 K^-1). Blood
#' flow to the skin raises it: about 0.5 with vessels near the skin constricted,
#' several times that with them dilated.
#' @param tair air temperature (degrees C)
#' @param pk atmospheric pressure (kPa)
#' @return resistance to heat (s m^-1)
#' @details In a sphere of radius R producing heat uniformly, the temperature
#' drop from centre to surface is that across a resistance R / 2k per unit area.
#' For an ellipsoid, R is taken as three times volume over surface area, which is
#' the radius for a sphere.
#' @export
tissueresistance <- function(volume, area, kflesh, tair, pk) {
  rcp <- rhoair(tair, pk) * cpair(tair)
  # Radius of the equivalent sphere (m), and conduction from centre to surface
  # with heat produced throughout
  R <- 3 * volume / area
  rt <- rcp * R / (2 * kflesh)
  return(rt)
}
#' @title Endotherm heat balance
#' @description Calculates what a bird or mammal must do to hold its core at body
#' temperature at a single time step, and the physiologically equivalent
#' temperature of its surroundings.
#' @param endo_params a list of endotherm parameters as returned by [create_endop()]
#' @param tme POSIXlt object indicating date and time in UTC
#' @param microcdata a list of microclimate variables with the following entries:
#'   \describe{
#'     \item{Tair}{Air temperature (degrees C).}
#'     \item{relhum}{Air relative humidity (percentage).}
#'     \item{pk}{Atmospheric pressure (kPa).}
#'     \item{windspeed}{Wind speed at animal height (m s^-1).}
#'     \item{winddir}{Wind direction (degrees from north).}
#'     \item{Rdirdown}{Direct shortwave radiation flux perpendicular to the solar beam (W m^-2).}
#'     \item{Rdifdown}{Diffuse downward shortwave radiation flux (W m^-2).}
#'     \item{Rswup}{Upward shortwave radiation flux (W m^-2).}
#'     \item{Rlwdown}{Downward longwave radiation flux (W m^-2).}
#'     \item{Rlwup}{Upward longwave radiation flux (W m^-2).}
#'   }
#' @param lat latitude (decimal degrees)
#' @param long longitude (decimal degrees)
#' @param uref wind speed of the reference setting for the physiologically
#' equivalent temperature (m s^-1). Default 0.1.
#' @param evapunits units in which evaporation is returned: \code{"W"} (heat
#' lost by evaporation, the default), \code{"L/day"} (water evaporated, litres
#' per day) or \code{"percent"} (water evaporated as a percentage of body mass
#' per day)
#' @return a list with:
#' \describe{
#'   \item{regime}{\code{"cold"}, \code{"thermoneutral"} or \code{"hot"}}
#'   \item{Mreq}{metabolic rate required to hold body temperature (W)}
#'   \item{Mbasal}{basal metabolic rate (W)}
#'   \item{Ereq}{evaporation required to hold body temperature, from breathing
#'   and skin together, in \code{evapunits}}
#'   \item{Emax}{most the animal can evaporate, panting at its fastest and, if it
#'   sweats, with its skin fully wet, in \code{evapunits}}
#'   \item{Tskin}{skin temperature (degrees C)}
#'   \item{Tsurf}{coat surface temperature (degrees C)}
#'   \item{Te}{operative temperature (degrees C)}
#'   \item{PET}{physiologically equivalent temperature (degrees C)}
#' }
#' @details The animal is a core held at body temperature, losing heat through
#' three resistances in series: tissue between core and skin, the coat, and the
#' surroundings (boundary layer and radiation in parallel). The animal regulates
#' the first through blood flow to the skin. Some heat is also lost by breathing,
#' as inhaled air leaves saturated at body temperature; the air breathed rises
#' with metabolic rate, 20 percent of the oxygen in it being used (as in
#' NicheMapR). If heat lost with blood flow at its lowest still exceeds what basal
#' metabolism can supply, the animal is cold and must raise its metabolic rate. If
#' heat lost with blood flow at its highest is less than basal metabolism, it is
#' hot and must evaporate the rest. Between the two it is thermoneutral. Panting
#' raises the rate of breathing up to \code{pantmax} times; an animal that sweats
#' can wet its whole skin, and vapour then passes through the coat and boundary
#' layer, which are taken to resist vapour as they resist heat.
#'
#' The physiologically equivalent temperature is the air temperature of a
#' reference setting (shade, surroundings at air temperature, wind \code{uref}) in
#' which the animal would have the same core and skin temperature, and so the same
#' heat flux from the skin, as in the actual conditions. Evaporation is the same in
#' both settings, so humidity does not enter it.
#' @export
Endotherm <- function(endo_params, tme, microcdata, lat, long, uref = 0.1, evapunits = "W") {
  p <- endo_params
  Tair <- microcdata$Tair
  pk <- microcdata$pk
  rcp <- rhoair(Tair, pk) * cpair(Tair)

  # Body size, the radiation the body absorbs and its size as the wind meets it
  solp <- sunposition(tme, lat, long)
  geo <- silhouette(solp$zen, solp$azi, p$height, p$width, p$len, p$adir, p$atilt, p$position)
  A <- geo$SurfArea
  Rabs <- animalRabs(p, solp, microcdata)
  d <- chardim(microcdata$winddir, p$height, p$width, p$len, p$adir, p$atilt, p$position,
               solp$zen, solp$azi)

  # Heat lost by breathing, as a share of metabolic rate: air breathed per joule
  # of metabolism (20.1 kJ per litre of oxygen, 20 percent of the oxygen used)
  # leaves saturated at body temperature, carrying the latent heat of the water
  vapdens <- function(e, tc) e * 1000 * 0.018015 / (8.314 * (tc + 273.15)) # kg m^-3
  ea <- satvap(Tair) * microcdata$relhum / 100
  airperJ <- 1 / (20.1e6 * 0.2095 * 0.2) # m^3 of air per J
  lakg <- latvap(p$Tb) / 0.018015 # latent heat (J kg^-1)
  fresp <- airperJ * lakg * max(vapdens(satvap(p$Tb), p$Tb) - vapdens(ea, Tair), 0)

  # Resistances of coat and of tissue at lowest and highest blood flow, and basal
  # metabolism per unit body surface
  rc <- coatresistance(p$coat, p$kcoat, Tair, pk)
  rtmax <- tissueresistance(geo$Volume, A, p$kmin, Tair, pk)
  rtmin <- tissueresistance(geo$Volume, A, p$kmax, Tair, pk)
  Mb <- p$basmet / A

  # Heat balance. The boundary layer depends on the temperature of the coat
  # surface, so the balance is repeated until that settles
  Tsurf <- Tair
  for (i in 1:10) {
    rHa <- animalrHa(Tair, Tsurf - Tair, microcdata$windspeed, d)
    ot <- operativetemp(Rabs, Tair, rHa, p$em, pk)
    Te <- ot$Te
    re <- ot$re
    # Dry heat loss for a given tissue resistance (W m^-2)
    Q <- function(rt) rcp * (p$Tb - Te) / (rt + rc + re)
    if (Q(rtmax) >= Mb * (1 - fresp)) {
      # Cold: blood flow at its lowest, metabolism makes up the loss
      regime <- "cold"
      rt <- rtmax
      Mreq <- Q(rt) / (1 - fresp)
    } else if (Q(rtmin) >= Mb * (1 - fresp)) {
      # Thermoneutral: blood flow set so that basal metabolism balances the loss
      regime <- "thermoneutral"
      rt <- rcp * (p$Tb - Te) / (Mb * (1 - fresp)) - rc - re
      Mreq <- Mb
    } else {
      # Hot: blood flow at its highest, evaporation removes the rest
      regime <- "hot"
      rt <- rtmin
      Mreq <- Mb
    }
    Ereq <- Mreq - Q(rt)
    # Skin and coat surface temperatures, from the share of the drop from core to
    # surroundings that lies outside each
    Tsurf <- Te + (p$Tb - Te) * re / (rt + rc + re)
    Tskin <- Te + (p$Tb - Te) * (rc + re) / (rt + rc + re)
  }

  # Most the animal can evaporate (W m^-2): breathing at its fastest, and vapour
  # from the wetted share of the skin through coat and boundary layer
  wet <- if (p$sweat) 1 else p$wetfrac
  Eskin <- wet * rhoair(Tair, pk) * latvap(Tskin) * max(satvap(Tskin) - ea, 0) / ((rc + rHa) * pk)
  Emax <- p$pantmax * fresp * Mreq + Eskin

  # Physiologically equivalent temperature: the still-air temperature that drives
  # the same heat flux from skin at the same temperature through coat and
  # reference boundary layer
  Qdry <- Q(rt)
  PET <- Tair
  Tsr <- Tsurf
  for (i in 1:10) {
    rHr <- animalrHa(PET, Tsr - PET, uref, d)
    rer <- operativetemp(p$em * 5.67e-8 * (PET + 273.15)^4, PET, rHr, p$em, pk)$re
    PET <- Tskin - Qdry * (rc + rer) / rcp
    Tsr <- PET + Qdry * rer / rcp
  }

  # Evaporation for the whole animal, in the units asked for
  E <- c(Ereq, Emax) * A
  if (evapunits == "L/day") {
    E <- E / lakg * 86400
  } else if (evapunits == "percent") {
    E <- 100 * E / lakg * 86400 / (geo$Volume * p$rho)
  } else if (evapunits != "W") stop("evapunits must be W, L/day or percent")

  out <- list(regime = regime, Mreq = Mreq * A, Mbasal = p$basmet, Ereq = E[1], Emax = E[2],
              Tskin = Tskin, Tsurf = Tsurf, Te = Te, PET = PET)
  return(out)
}
#' @title Endotherm parameters
#' @description Creates a list of parameters for [Endotherm()] for a bird or
#' mammal from its size and coat. Physiological values are typical ones from the
#' literature, so that only things that can be seen or measured need be given.
#' Any value can be changed in the returned list.
#' @param height body height (cm)
#' @param width body width (cm)
#' @param len body length (cm)
#' @param coat depth of fur or plumage (cm)
#' @param type one of \code{"mammal"} or \code{"bird"}
#' @param sweat logical: can the animal sweat? Default \code{FALSE}.
#' @param refl shortwave reflectance of the coat. Default 0.2.
#' @param adir body azimuth direction (degrees). Default 0.
#' @param atilt body tilt (degrees). See [silhouette()]. Default 0.
#' @param position assumption about body orientation. See [silhouette()].
#' Default \code{"randomdir"}.
#' @return a list of parameters:
#' \describe{
#'   \item{height, width, len, coat, sweat, adir, atilt, position, refl}{as given}
#'   \item{em}{longwave emissivity of the coat: 0.97}
#'   \item{rho}{body density (kg m^-3): 1000}
#'   \item{Tb}{core body temperature (degrees C): 37 for mammals, as in NicheMapR;
#'   41 for birds, the mean for the active phase across 724 species (Prinzinger et
#'   al. 1991)}
#'   \item{basmet}{basal metabolic rate (W), from body mass by Kleiber's relation,
#'   70 kcal per day times mass in kg to the power 0.75}
#'   \item{kmin, kmax}{effective flesh conductivity with least and most blood flow
#'   to the skin (W m^-1 K^-1): 0.9 and 2.8, NicheMapR's values}
#'   \item{kcoat}{effective conductivity of fur or plumage (W m^-1 K^-1): 0.05,
#'   measured for cattle hair coats (Maia et al. 2009). Scholander et al. (1950)
#'   found insulation roughly in proportion to coat depth across arctic and
#'   tropical mammals, so one value serves.}
#'   \item{wetfrac}{share of the skin that is wet without sweating: 0.005, as in
#'   NicheMapR}
#'   \item{pantmax}{most by which panting raises the rate of breathing: 5, as in
#'   NicheMapR}
#' }
#' @references Kearney MR, Briscoe NJ, Mathewson PD, Porter WP (2021) NicheMapR:
#' an R package for biophysical modelling: the endotherm model. Ecography 44:
#' 1595-1605.
#'
#' Kleiber M (1947) Body size and metabolic rate. Physiological Reviews 27:
#' 511-541.
#'
#' Maia ASC, da Silva RG, de Souza Junior JBF, da Silva RB, Domingos HGT (2009)
#' Effective thermal conductivity of the hair coat of Holstein cows in a
#' tropical environment. Revista Brasileira de Zootecnia 38: 2218-2223.
#'
#' Prinzinger R, Pressmar A, Schleucher E (1991) Body temperature in birds.
#' Comparative Biochemistry and Physiology 99A: 499-506.
#'
#' Scholander PF, Walters V, Hock R, Irving L (1950) Body insulation of some
#' arctic and tropical mammals and birds. Biological Bulletin 99: 225-236.
#' @export
create_endop <- function(height, width, len, coat, type = "mammal", sweat = FALSE, refl = 0.2,
                         adir = 0, atilt = 0, position = "randomdir") {
  if (!type %in% c("mammal", "bird")) stop("type must be mammal or bird")
  # Body mass from the volume of the ellipsoid, and basal metabolism from mass
  # (70 kcal per day is 3.39 W)
  volume <- (4 / 3) * pi * (height / 200) * (width / 200) * (len / 200)
  basmet <- 3.39 * (volume * 1000)^0.75
  out <- list(height = height, width = width, len = len, coat = coat, sweat = sweat,
              adir = adir, atilt = atilt, position = position, refl = refl, em = 0.97,
              rho = 1000, Tb = if (type == "bird") 41 else 37, basmet = basmet,
              kmin = 0.9, kmax = 2.8, kcoat = 0.05, wetfrac = 0.005, pantmax = 5)
  return(out)
}
