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
#' @param Tref Reference temperature (degrees C)
#' @param Tbody Body temperature (degrees C)
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
  # diffusion, across half the height of the body. Only the rest of the body
  # exchanges radiation
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
    top <- Rabs + M - mr * (Tair + 273.15)^4 + mc * (Tsurf - Tair) - ml * Da - mlc * Dac
    btm <- 4 * mr * (Te + 273.15)^3 + mc + mv + (ml + mlc) * DeV
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
#' \code{position} (see [silhouette()]), \code{refl} (shortwave reflectance),
#' \code{em} (longwave emissivity) and \code{confrac} (share of the body surface
#' in contact with the ground)
#' @param solp a list with the solar zenith (\code{zen}) and azimuth (\code{azi})
#' angles (degrees), as returned by [sunposition()]
#' @param microcdata a list of radiation fluxes at the animal (W m^-2):
#' \code{Rdirdown} (direct, perpendicular to the solar beam), \code{Rdifdown}
#' (diffuse downward), \code{Rswup} (upward shortwave), \code{Rlwdown} and
#' \code{Rlwup} (downward and upward longwave)
#' @return radiation absorbed per unit body surface (W m^-2)
#' @details The direct beam falls on the silhouette the body presents to the sun.
#' Half the body faces upward and receives diffuse and longwave radiation from the
#' sky. Half faces downward and receives radiation reflected or emitted by the
#' ground, except where it touches the ground, which receives none.
#' @export
animalRabs <- function(params, solp, microcdata) {
  # Share of the body surface that the direct beam strikes
  sil <- with(params, silhouette(solp$zen, solp$azi, height, width, len, adir, atilt, position))
  sc <- sil$silA / sil$SurfArea
  # Share of the body surface facing the ground without touching it
  down <- max(0.5 - params$confrac, 0)
  # Shortwave absorbed from the direct beam, diffuse skylight and light reflected
  # by the ground, and longwave from the sky and the ground
  Rsw <- sc * microcdata$Rdirdown + 0.5 * microcdata$Rdifdown + down * microcdata$Rswup
  Rabs <- (1 - params$refl) * Rsw + params$em * (0.5 * microcdata$Rlwdown + down * microcdata$Rlwup)
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
# Heat lost by evaporation (W m^-2) into air at temperature ta and vapour pressure
# e, for core and skin at temperatures Tc and Tsk and metabolic rate M. Vapour
# leaves the skin through coat and boundary layer, and breathing removes it as
# inhaled air leaves saturated at core temperature (20.1 kJ per litre of oxygen,
# 20 percent of the oxygen used). The evaporation the animal seeks, Ed, comes from
# sweat or from panting, each limited by what the air can take up
.endoevap <- function(p, Tc, Tsk, Ed, ta, e, rHa, rc, M, pk) {
  vapdens <- function(e, tc) e * 1000 * 0.018015 / (8.314 * (tc + 273.15)) # kg m^-3
  airperJ <- 1 / (20.1e6 * 0.2095 * 0.2) # m^3 of air per J
  Eskinmax <- (1 - p$confrac) * rhoair(ta, pk) * latvap(Tsk) * max(satvap(Tsk) - e, 0) /
    ((rc + rHa) * pk)
  Eresp0 <- M * airperJ * latvap(Tc) / 0.018015 *
    max(vapdens(satvap(Tc), Tc) - vapdens(e, ta), 0)
  if (p$sweat) {
    Eskin <- min(max(Ed, p$wetfrac * Eskinmax), Eskinmax)
    Eresp <- Eresp0
  } else {
    Eskin <- p$wetfrac * Eskinmax
    Eresp <- min(Eresp0 + Ed, p$pantmax * Eresp0)
  }
  c(skin = Eskin, resp = Eresp)
}
#' @title Endotherm heat balance
#' @description Calculates what a bird or mammal must do to hold its core at body
#' temperature at a single time step: the metabolic heat it must produce when it
#' is cold, or the heat it must lose by evaporation when it is hot.
#' @param endo_params a list of endotherm parameters as returned by [create_endop()]
#' @param tme POSIXlt object indicating date and time in UTC
#' @param microcdata a list of microclimate variables with the following entries:
#'   \describe{
#'     \item{Tair}{Air temperature (degrees C).}
#'     \item{Tsurf}{Temperature of the ground surface the animal stands or lies on (degrees C).}
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
#' @param evapunits units in which evaporation is returned: \code{"W"} (heat
#' lost by evaporation, the default), \code{"L/day"} (water evaporated, litres
#' per day) or \code{"percent"} (water evaporated as a percentage of body mass
#' per day)
#' @return a list with:
#' \describe{
#'   \item{regime}{\code{"cold"}, \code{"thermoneutral"} or \code{"hot"}}
#'   \item{Mreq}{metabolic rate required to hold body temperature (W)}
#'   \item{Mbasal}{basal metabolic rate (W)}
#'   \item{Mmax}{the highest metabolic rate the animal can sustain (W)}
#'   \item{Lreq}{heat that must be lost by evaporation to hold body temperature,
#'   in \code{evapunits}}
#'   \item{Lmax}{the most the animal can evaporate, in \code{evapunits}}
#'   \item{Tskin}{skin temperature (degrees C)}
#'   \item{Tsurf}{coat surface temperature (degrees C)}
#'   \item{Te}{operative temperature (degrees C)}
#' }
#' @details The animal is a core held at body temperature. Heat passes through
#' the tissue between core and skin, which the animal regulates through blood flow
#' to the skin, and then by two routes side by side: over most of the body through
#' the coat to the surroundings (boundary layer and radiation in parallel), and
#' where the body touches the ground through the layer between skin and ground to
#' the ground surface. If heat lost with blood flow at its lowest still exceeds
#' basal metabolism, the animal is cold and must raise its metabolic rate. If heat
#' lost with blood flow at its highest is less than basal metabolism, it is hot
#' and must lose the rest by evaporation. Between the two it is thermoneutral.
#'
#' The highest metabolic rate is \code{metmax} times basal. The most the animal
#' can evaporate is from panting at its fastest (\code{pantmax} times the resting
#' rate of breathing) and from its skin, fully wet if it sweats; it depends on air
#' humidity. A warning is given when the animal would need more than either, as it
#' then cannot hold its body temperature.
#' @seealso [EndothermPET()] for how hot or cold the surroundings feel to the
#' animal
#' @export
Endotherm <- function(endo_params, tme, microcdata, lat, long, evapunits = "W") {
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

  # Resistances of coat, of the layer between skin and ground where the body
  # touches it, and of tissue at lowest and highest blood flow, and basal
  # metabolism per unit body surface
  rc <- coatresistance(p$coat, p$kcoat, Tair, pk)
  rg <- coatresistance(p$contactdepth, p$kcontact, Tair, pk)
  rtmax <- tissueresistance(geo$Volume, A, p$kmin, Tair, pk)
  rtmin <- tissueresistance(geo$Volume, A, p$kmax, Tair, pk)
  Mb <- p$basmet / A
  cf <- p$confrac

  # Heat balance. The boundary layer depends on the temperature of the coat
  # surface, so the balance is repeated until that settles
  Tsurf <- Tair
  for (i in 1:10) {
    # Operative temperature of the exposed coat, from the radiation it absorbs
    rHa <- animalrHa(Tair, Tsurf - Tair, microcdata$windspeed, d)
    ot <- operativetemp(Rabs / (1 - cf), Tair, rHa, p$em, pk)
    Te <- ot$Te
    re <- ot$re
    # The routes from skin through the coat and through the contact with the
    # ground act together as one resistance to surroundings at a temperature
    # between the operative temperature and that of the ground
    ro <- 1 / ((1 - cf) / (rc + re) + cf / rg)
    Tenv <- ro * ((1 - cf) * Te / (rc + re) + cf * microcdata$Tsurf / rg)
    # Dry heat loss for a given tissue resistance (W m^-2)
    Q <- function(rt) rcp * (p$Tb - Tenv) / (rt + ro)
    if (Q(rtmax) >= Mb) {
      # Cold: blood flow at its lowest, metabolism makes up the loss
      regime <- "cold"
      rt <- rtmax
      Mreq <- Q(rt)
    } else if (Q(rtmin) >= Mb) {
      # Thermoneutral: blood flow set so that basal metabolism balances the loss
      regime <- "thermoneutral"
      rt <- rcp * (p$Tb - Tenv) / Mb - ro
      Mreq <- Mb
    } else {
      # Hot: blood flow at its highest, evaporation removes the rest
      regime <- "hot"
      rt <- rtmin
      Mreq <- Mb
    }
    Lreq <- Mreq - Q(rt)
    # Skin temperature, from the share of the drop from core to surroundings that
    # lies outside the tissue, and coat surface temperature, from the share of the
    # drop from skin to operative temperature that lies outside the coat
    Tskin <- Tenv + (p$Tb - Tenv) * ro / (rt + ro)
    Tsurf <- Te + (Tskin - Te) * re / (rc + re)
  }

  # The animal's limits: its highest sustainable metabolic rate, and the most it
  # can evaporate into the surrounding air at basal metabolism by panting at its
  # fastest and, if it sweats, from fully wet skin
  Mmax <- p$metmax * p$basmet
  ea <- satvap(Tair) * microcdata$relhum / 100
  Lmax <- sum(.endoevap(p, p$Tb, Tskin, Inf, Tair, ea, rHa, rc, Mb, pk))
  if (Mreq * A > Mmax) {
    warning("metabolic rate required exceeds the animal's maximum: it cannot hold its body temperature")
  }
  if (Lreq > Lmax) {
    warning("evaporation required exceeds the most the animal can achieve: it cannot hold its body temperature")
  }

  # Evaporation for the whole animal, in the units asked for
  L <- c(Lreq, Lmax) * A
  lakg <- latvap(p$Tb) / 0.018015 # latent heat (J kg^-1)
  if (evapunits == "L/day") {
    L <- L / lakg * 86400
  } else if (evapunits == "percent") {
    L <- 100 * L / lakg * 86400 / (geo$Volume * p$rho)
  } else if (evapunits != "W") stop("evapunits must be W, L/day or percent")

  out <- list(regime = regime, Mreq = Mreq * A, Mbasal = p$basmet, Mmax = Mmax,
              Lreq = L[1], Lmax = L[2], Tskin = Tskin, Tsurf = Tsurf, Te = Te)
  return(out)
}
#' @title Physiologically equivalent temperature of an endotherm
#' @description Calculates how hot or cold the surroundings feel to a bird or
#' mammal at a single time step, as a physiologically equivalent temperature
#' (PET): the air temperature of a reference setting in which the animal would
#' have the same core and skin temperatures as in its actual surroundings.
#' @param endo_params a list of endotherm parameters as returned by [create_endop()]
#' @param tme POSIXlt object indicating date and time in UTC
#' @param microcdata a list of microclimate variables as for [Endotherm()]
#' @param lat latitude (decimal degrees)
#' @param long longitude (decimal degrees)
#' @param metmult metabolic rate of the animal, as a multiple of its basal rate.
#' Default 2.
#' @param uref wind speed of the reference setting (m s^-1). Default 0.1.
#' @param eref vapour pressure of the reference setting (kPa). Default 1.2.
#' @return a list with:
#' \describe{
#'   \item{PET}{physiologically equivalent temperature (degrees C)}
#'   \item{Tcore}{core temperature in the actual surroundings (degrees C)}
#'   \item{Tskin}{skin temperature in the actual surroundings (degrees C)}
#' }
#' @details The calculation follows Hoppe (1999) in two steps. First, metabolism
#' is held at \code{metmult} times its basal rate and the core and skin settle at
#' whatever temperatures balance the heat the animal makes against what it loses
#' in its actual surroundings. The body is that of [Endotherm()]: core, tissue,
#' coat and surroundings, with contact to the ground. Blood flow to the skin and
#' evaporation respond to core and skin temperature by the rules of Hoppe's model
#' of the human body, taken relative to the animal's normal temperatures: its body
#' temperature \code{Tb} for the core, and for the skin the temperature it has at
#' basal metabolism with blood flow midway between its limits. Blood flow rises as
#' the core warms and falls as the skin cools, within the limits set by
#' \code{kmin} and \code{kmax}. The water the animal seeks to evaporate rises with
#' body temperature. An animal that sweats evaporates it from its skin, and one
#' that does not pants, up to \code{pantmax} times its resting rate of breathing;
#' either way no faster than the coat, boundary layer and humidity allow. Skin that
#' is always wet (\code{wetfrac}) and breathing at rest also lose water. A warning
#' is given when core temperature is above \code{Tcmax} or below \code{Tcmin},
#' which the animal would not survive.
#'
#' Second, the animal is placed in the reference setting: shade, with sky and
#' ground at air temperature, wind \code{uref} and vapour pressure \code{eref}.
#' Core and skin keep their temperatures, so blood flow and the water the animal
#' seeks to evaporate are unchanged, and metabolism stays the same. PET is the air
#' temperature there at which the animal loses heat as fast as it makes it.
#' @references Hoppe P (1999) The physiological equivalent temperature: a
#' universal index for the biometeorological assessment of the thermal
#' environment. International Journal of Biometeorology 43: 71-75.
#' @seealso [Endotherm()] for what the animal must do to hold its body
#' temperature constant
#' @export
EndothermPET <- function(endo_params, tme, microcdata, lat, long, metmult = 2, uref = 0.1,
                         eref = 1.2) {
  p <- endo_params
  Tair <- microcdata$Tair
  pk <- microcdata$pk
  rcp <- rhoair(Tair, pk) * cpair(Tair)
  ea <- satvap(Tair) * microcdata$relhum / 100

  # Body size, the radiation the body absorbs and its size as the wind meets it
  solp <- sunposition(tme, lat, long)
  geo <- silhouette(solp$zen, solp$azi, p$height, p$width, p$len, p$adir, p$atilt, p$position)
  A <- geo$SurfArea
  Rabs <- animalRabs(p, solp, microcdata)
  d <- chardim(microcdata$winddir, p$height, p$width, p$len, p$adir, p$atilt, p$position,
               solp$zen, solp$azi)

  # Resistances of coat and of the layer between skin and ground, conductance of
  # the tissue from core to skin with vessels near the skin constricted and
  # dilated (W m^-2 K^-1), and basal and actual metabolism per unit body surface
  rc <- coatresistance(p$coat, p$kcoat, Tair, pk)
  rg <- coatresistance(p$contactdepth, p$kcontact, Tair, pk)
  R <- 3 * geo$Volume / A
  Kmin <- 2 * p$kmin / R
  Kmax <- 2 * p$kmax / R
  Mb <- p$basmet / A
  M <- metmult * Mb
  cf <- p$confrac

  # Normal temperatures: core at body temperature, and skin at the temperature
  # basal metabolism gives with blood flow midway between its limits
  Tcset <- p$Tb
  Tskset <- p$Tb - Mb * (1 / Kmin + 1 / Kmax) / 2

  # Responses from Hoppe's model of the human body. Blood flow to the skin
  # (L m^-2 h^-1) rises as the core warms and falls as the skin cools, each
  # litre carrying 3640 J per degree from core to skin. The sweat sought
  # (g m^-2 h^-1) rises with a mean body temperature weighted to the core, each
  # gram taking 2420 J to evaporate
  respond <- function(Tc, Tsk) {
    warm <- max(Tc - Tcset, 0)
    cool <- max(Tskset - Tsk, 0)
    blood <- min((6.3 + 75 * warm) / (1 + 0.5 * cool), 90)
    K <- min(Kmin + blood * 3640 / 3600, Kmax)
    a <- 0.0418 + 0.745 / (blood + 0.585) # share of the body at skin temperature
    Tbody <- a * Tsk + (1 - a) * Tc
    sweat <- min(304.94 * max(Tbody - (0.1 * Tskset + 0.9 * Tcset), 0), 500)
    list(K = K, Ed = 2420 * sweat / 3600)
  }

  # Step 1, actual surroundings. The core balances when the heat it makes, less
  # that lost by breathing, crosses the tissue to the skin; the skin balances
  # when that heat leaves it by dry exchange and evaporation. For each trial core
  # temperature the skin temperature that balances the skin is found, and the
  # core temperature is then adjusted until the core balances too. The boundary
  # layer depends on coat surface temperature, so all is repeated until that
  # settles
  Tsurf <- Tair
  for (i in 1:5) {
    rHa <- animalrHa(Tair, Tsurf - Tair, microcdata$windspeed, d)
    ot <- operativetemp(Rabs / (1 - cf), Tair, rHa, p$em, pk)
    ro <- 1 / ((1 - cf) / (rc + ot$re) + cf / rg)
    Tenv <- ro * ((1 - cf) * ot$Te / (rc + ot$re) + cf * microcdata$Tsurf / rg)
    skin <- function(Tc) {
      skinbalance <- function(Tsk) {
        r <- respond(Tc, Tsk)
        E <- .endoevap(p, Tc, Tsk, r$Ed, Tair, ea, rHa, rc, M, pk)
        r$K * (Tc - Tsk) - rcp * (Tsk - Tenv) / ro - E[["skin"]]
      }
      uniroot(skinbalance, c(-80, 80), tol = 1e-4)$root
    }
    corebalance <- function(Tc) {
      Tsk <- skin(Tc)
      r <- respond(Tc, Tsk)
      E <- .endoevap(p, Tc, Tsk, r$Ed, Tair, ea, rHa, rc, M, pk)
      M - E[["resp"]] - r$K * (Tc - Tsk)
    }
    Tc <- uniroot(corebalance, c(-80, 80), tol = 1e-4)$root
    Tsk <- skin(Tc)
    Tsurf <- ot$Te + (Tsk - ot$Te) * ot$re / (rc + ot$re)
  }
  if (Tc > p$Tcmax) {
    warning("core temperature of ", round(Tc, 1), " degrees C is above what the animal can survive")
  }
  if (Tc < p$Tcmin) {
    warning("core temperature of ", round(Tc, 1), " degrees C is below what the animal can survive")
  }

  # Step 2, reference setting. With core and skin at the same temperatures, PET
  # is the air temperature at which dry exchange and evaporation together remove
  # the heat the animal makes. Surroundings at air temperature make the
  # operative temperature equal to it
  r <- respond(Tc, Tsk)
  PET <- Tair
  Tsr <- Tsurf
  for (i in 1:20) {
    rHr <- animalrHa(PET, Tsr - PET, uref, d)
    rer <- operativetemp(p$em * 5.67e-8 * (PET + 273.15)^4, PET, rHr, p$em, pk)$re
    ror <- 1 / ((1 - cf) / (rc + rer) + cf / rg)
    E <- .endoevap(p, Tc, Tsk, r$Ed, PET, eref, rHr, rc, M, pk)
    PET <- Tsk - ror * (M - sum(E)) / rcp
    Tsr <- PET + (Tsk - PET) * rer / (rc + rer)
  }
  return(list(PET = PET, Tcore = Tc, Tskin = Tsk))
}
#' @title Endotherm parameters
#' @description Creates a list of parameters for [Endotherm()] and
#' [EndothermPET()] for a bird or
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
#' @param posture one of \code{"standing"} (the default), touching the ground
#' through the feet, or \code{"lying"}, touching it with the underside of the body
#' @return a list of parameters:
#' \describe{
#'   \item{height, width, len, coat, sweat, adir, atilt, position, refl}{as given}
#'   \item{confrac}{share of the body surface in contact with the ground: 0.02
#'   standing, the soles of the feet; 0.2 lying}
#'   \item{contactdepth, kcontact}{length (cm) and conductivity (W m^-1 K^-1) of
#'   the path heat takes between skin and ground. Standing, the legs: as long as
#'   the body is high, with 0.5, the conductivity of flesh with little blood
#'   flow. Lying, the coat pressed to half its depth, with conductivity \code{kcoat}.}
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
#'   \item{metmax}{highest metabolic rate the animal can sustain, as a multiple of
#'   basal: 6. The highest rates of heat production measured in the cold are
#'   mostly 4 to 8 times basal in mammals, and similar in birds (Minnaar et al.
#'   2014).}
#'   \item{Tcmax, Tcmin}{highest and lowest core temperatures the animal can
#'   survive (degrees C). \code{Tcmax} is 5 degrees above the default \code{Tb}:
#'   42 for mammals and 46 for birds, whose core temperatures rarely exceed 46
#'   (Freeman et al. 2020). \code{Tcmin} is 20, below which the hearts of mammals
#'   that do not hibernate fail. Change \code{Tcmax} with \code{Tb}.}
#' }
#' @references Kearney MR, Briscoe NJ, Mathewson PD, Porter WP (2021) NicheMapR:
#' an R package for biophysical modelling: the endotherm model. Ecography 44:
#' 1595-1605.
#'
#' Freeman MT, Czenze ZJ, Schoeman K, McKechnie AE (2020) Extreme hyperthermia
#' tolerance in the world's most abundant wild bird. Scientific Reports 10: 13098.
#'
#' Kleiber M (1947) Body size and metabolic rate. Physiological Reviews 27:
#' 511-541.
#'
#' Minnaar IA, Bennett NC, Chimimba CT, McKechnie AE (2014) Summit metabolism and
#' metabolic expansibility in Wahlberg's epauletted fruit bats (Epomophorus
#' wahlbergi): seasonal acclimatisation and effects of captivity. Journal of
#' Experimental Biology 217: 1363-1369.
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
                         adir = 0, atilt = 0, position = "randomdir", posture = "standing") {
  if (!type %in% c("mammal", "bird")) stop("type must be mammal or bird")
  if (!posture %in% c("standing", "lying")) stop("posture must be standing or lying")
  # Body mass from the volume of the ellipsoid, and basal metabolism from mass
  # (70 kcal per day is 3.39 W)
  volume <- (4 / 3) * pi * (height / 200) * (width / 200) * (len / 200)
  basmet <- 3.39 * (volume * 1000)^0.75
  kcoat <- 0.05
  # Contact with the ground: along the legs to the soles of the feet when
  # standing, through the compressed coat of the underside when lying
  if (posture == "standing") {
    confrac <- 0.02
    contactdepth <- height
    kcontact <- 0.5
  } else {
    confrac <- 0.2
    contactdepth <- 0.5 * coat
    kcontact <- kcoat
  }
  out <- list(height = height, width = width, len = len, coat = coat, sweat = sweat,
              adir = adir, atilt = atilt, position = position, refl = refl, em = 0.97,
              rho = 1000, Tb = if (type == "bird") 41 else 37, basmet = basmet,
              kmin = 0.9, kmax = 2.8, kcoat = kcoat, wetfrac = 0.005, pantmax = 5,
              metmax = 6, confrac = confrac, contactdepth = contactdepth,
              kcontact = kcontact)
  # Survivable core temperatures
  out$Tcmax <- out$Tb + 5
  out$Tcmin <- 20
  return(out)
}
