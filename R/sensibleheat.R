#' Calculate molar density of air
#'
#' @param tc air temperature (degrees C)
#' @param pk atmospheric pressure (kPa)
#' @returns molar density of air (mol/m3).
#' @export
rhoair <- function(tc, pk) {
  tk <- tc + 273.15
  ph <- 44.6 * (pk / 101.3) * (273.15 / tk)
  return(ph)
}

#' Calculate specific heat of air at constant pressure
#'
#' @param tc air temperature (degrees C)
#' @returns specific heat of air at constant pressure (J/mol/K)
#' @export
cpair <- function(tc) {
  cp <- 2e-05 * tc^2 + 0.0002 * tc + 29.119
  return(cp)
}

#' Calculate zero-plane displacement height
#' @param h vegetation height (m)
#' @param pai plant area index
#' @returns zero-plane displacement height (m)
#' @export
zeroplanedis <- function(h, pai) {
  pai[pai < 1e-9] <- 1e-9
  d <- (1 - (1 - exp(-sqrt(7.5 * pai))) / sqrt(7.5 * pai)) * h
  return(d)
}

# How closed a canopy is: 0 with no foliage, 1 when more foliage adds no more
# drag. The effects of the roughness sublayer are scaled by it, so that they are
# absent over bare ground and complete over a closed canopy
.canopyclosure <- function(pai) {
  Be <- pmin(sqrt(0.003 + 0.15 * pai), 0.3)
  (Be - sqrt(0.003)) / (0.3 - sqrt(0.003))
}
#' Calculate roughness length for momentum
#' @param h vegetation height (m)
#' @param pai plant area index
#' @param d zero-plane displacement height as returned by [zeroplanedis()] (m)
#' @param MOST Logical. If `TRUE`, the wind profile is taken to be logarithmic
#' down to canopy top. If `FALSE`, a roughness sublayer is assumed above the
#' canopy.
#' @param z0 roughness length of the surface beneath the vegetation (m). The
#' default is that of bare ground.
#' @returns roughness length (m)
#' @export
roughlength <- function(h, pai, d, MOST = TRUE, z0 = 0.004) {
  # Ratio of friction velocity to wind speed at canopy top. It rises with plant
  # area until the canopy is closed, when more foliage adds no more drag
  Be <- pmin(sqrt(0.003 + 0.15 * pai), 0.3)
  # Roughness length that gives this ratio at canopy top
  zm <- (h - d) * exp(-0.41 / Be)
  # Wind just above a canopy is stronger than a logarithmic profile gives. The
  # roughness length is raised to allow for this, by a factor of exp(0.193) over
  # a closed canopy and by less where foliage is sparse
  if (!MOST) {
    cw <- 1 + .canopyclosure(pai)
    zm <- zm * exp(log(cw) - 1 + 1 / cw)
  }
  # A canopy is no rougher than its height above the displacement height allows,
  # and no smoother than the surface beneath it, which is all that is left when
  # there is no vegetation
  zm <- pmax(pmin(zm, 0.9 * (h - d)), z0)
  return(zm)
}
#' Calculate integrated diabatic correction coefficient for momentum
#' @param ze atmospheric stability parameter (height divided by Obukhov length)
#' @returns diabatic correction coefficient for momentum
#' @export
dpsim <- function(ze) {
  psim <- rep(NA_real_, length(ze))
  un <- which(ze < 0)
  st <- which(ze >= 0)
  # Unstable air: rising warm air adds to mixing, up to a limit
  x <- (1 - 15 * ze[un])^0.25
  psim[un] <- pmin(log(((1 + x) / 2)^2 * (1 + x^2) / 2) - 2 * atan(x) + pi / 2, 3)
  # Stable air: mixing is suppressed, approaching a limit smoothly
  psim[st] <- -4 * tanh(4.7 * ze[st] / 4)
  return(psim)
}
#' Calculate integrated diabatic correction coefficient for heat
#' @param ze atmospheric stability parameter (height divided by Obukhov length)
#' @returns diabatic correction coefficient for heat
#' @export
dpsih <- function(ze) {
  psih <- rep(NA_real_, length(ze))
  un <- which(ze < 0)
  st <- which(ze >= 0)
  # Unstable air: rising warm air adds to mixing, approaching a limit smoothly
  y <- sqrt(1 - 9 * ze[un])
  ps <- 2 * log((1 + y) / 2)
  s <- ps > 1.5
  ps[s] <- 1.5 + 1.5 * tanh((ps[s] - 1.5) / 1.5)
  psih[un] <- ps
  # Stable air: mixing is suppressed, approaching a limit smoothly
  psih[st] <- -4 * tanh((4.7 / 0.74) * ze[st] / 4)
  return(psih)
}
#' Calculate diabatic influencing coefficient for heat
#'
#' The dimensionless temperature gradient that corresponds to [dpsih()]:
#' above 1 in stable air, where mixing is suppressed, and below 1 in unstable
#' air.
#' @param ze atmospheric stability parameter (height divided by Obukhov length)
#' @returns diabatic influencing coefficient for heat
#' @export
dphih <- function(ze) {
  phih <- rep(NA_real_, length(ze))
  un <- which(ze < 0)
  st <- which(ze >= 0)
  # Unstable air: returns towards 1 as the correction for heat nears its limit
  y <- sqrt(1 - 9 * ze[un])
  ph <- 1 / y
  ps <- 2 * log((1 + y) / 2)
  s <- ps > 1.5
  ph[s] <- 1 + (ph[s] - 1) / cosh((ps[s] - 1.5) / 1.5)^2
  phih[un] <- ph
  # Stable air: returns towards 1 as the correction for heat nears its limit
  a <- 4.7 / 0.74
  phih[st] <- 1 + a * ze[st] / cosh(a * ze[st] / 4)^2
  return(phih)
}
#' Calculate Obukhov length
#' @param tc air temperature (degrees C)
#' @param pk atmospheric pressure (kPa)
#' @param uf wind friction velocity (m/s)
#' @param H sensible heat flux (W/m^2)
#' @returns Obukhov length (m)
#' @export
Obukhov <- function(tc, pk, uf, H) {
  ph <- rhoair(tc, pk)
  cp <- cpair(tc)
  Tk <- tc + 273.15
  LL <- (ph * cp * uf^3 * Tk) / (-0.41 * 9.81 * H)
  return (LL)
}
#' Calculate effective wind speed when the surface heats the air
#'
#' When the surface heats the air above it, rising air stirs the surface layer
#' even in a calm. This function adds the velocity of that convection to the
#' measured wind speed (Beljaars 1994). The result is used in place of measured
#' wind speed when calculating friction velocity.
#' @param uref wind speed at reference height (m/s)
#' @param H sensible heat flux (W/m^2), positive when the surface heats the air
#' @param tc air temperature (degrees C)
#' @param pk atmospheric pressure (kPa)
#' @param zi depth of the layer of air stirred by convection (m)
#' @returns effective wind speed (m/s). Equal to `uref` when `H` is zero or
#' negative.
#' @export
effectivewind <- function(uref, H, tc, pk, zi = 1000) {
  ph <- rhoair(tc, pk)
  cp <- cpair(tc)
  Tk <- tc + 273.15
  # Velocity of the rising air, zero when heat flows from the air to the surface
  wstar <- ((9.81 / Tk) * zi * pmax(H, 0) / (ph * cp))^(1 / 3)
  # Measured wind and convection act together
  ueff <- sqrt(uref^2 + wstar^2)
  return(ueff)
}
# Amount by which wind in the roughness sublayer exceeds the logarithmic profile,
# as a multiple of friction velocity over 0.41. Over a closed canopy the sublayer
# reaches twice as far above the displacement height as canopy top does, and over
# sparse foliage less far. The excess is zero above it.
.sublayer <- function(z, h, d, pai) {
  zw <- (1 + .canopyclosure(pai)) * (h - d)
  zp <- pmin(z - d, zw)
  log(zw / zp) - 1 + zp / zw
}
#' Calculate wind friction velocity
#' @param uref wind speed at reference height (m/s)
#' @param zref height above ground of wind speed measurement (m)
#' @param h vegetation height (m)
#' @param pai plant area index
#' @param LL Obukhov length as returned by [Obukhov()] (m). The default, `Inf`,
#' is neutral air.
#' @param MOST Logical. If `TRUE`, profiles above the canopy are taken to be
#' logarithmic down to canopy top. If `FALSE`, a roughness sublayer is assumed
#' above the canopy.
#' @param z0 roughness length of the surface beneath the vegetation (m). The
#' default is that of bare ground.
#' @returns friction velocity (m/s)
#' @export
windfric <- function(uref, zref, h, pai, LL = Inf, MOST = TRUE, z0 = 0.004) {
  # Bare ground has no canopy to form a roughness sublayer
  if (!MOST && h == 0) {
    warning("No roughness sublayer over bare ground (h = 0): using MOST = TRUE")
    MOST <- TRUE
  }
  d <- zeroplanedis(h, pai)
  zm <- roughlength(h, pai, d, MOST, z0)
  # Diabatic correction for momentum between the surface and the reference height
  psim <- dpsim(zm / LL) - dpsim((zref - d) / LL)
  psis <- 0
  if (!MOST) psis <- .sublayer(zref, h, d, pai)
  uf <- (0.41 * uref) / (log((zref - d) / zm) + psim + psis)
  return(uf)
}
#' Calculate wind profile above canopy
#' @param z height above ground for which wind speed desired (m)
#' @param uf wind friction velocity as returned by [windfric()] (m/s)
#' @param h vegetation height (m)
#' @param pai plant area index
#' @param LL Obukhov length as returned by [Obukhov()] (m). The default, `Inf`,
#' is neutral air.
#' @param MOST Logical. If `TRUE`, profiles above the canopy are taken to be
#' logarithmic down to canopy top. If `FALSE`, a roughness sublayer is assumed
#' above the canopy.
#' @returns wind speed (m/s)
#' @export
windprofile_above <- function(z, uf, h, pai, LL = Inf, MOST = TRUE) {
  # Bare ground has no canopy to form a roughness sublayer
  if (!MOST && h == 0) {
    warning("No roughness sublayer over bare ground (h = 0): using MOST = TRUE")
    MOST <- TRUE
  }
  d <- zeroplanedis(h, pai)
  zm <- roughlength(h, pai, d, MOST)
  # Diabatic correction for momentum between the surface and each height
  psim <- dpsim(zm / LL) - dpsim((z - d) / LL)
  psis <- 0
  if (!MOST) psis <- .sublayer(z, h, d, pai)
  uz <- (uf / 0.41) * (log((z - d) / zm) + psim + psis)
  return(uz)
}
# Factor by which thermal diffusivity at canopy top exceeds what the logarithmic
# law gives there: three over a closed canopy, one over bare ground
.sublayerfactor <- function(pai) {
  1 + 2 * .canopyclosure(pai)
}
# Height above the displacement height at which the logarithmic law for heat is
# applied when there is a roughness sublayer. It exceeds the true height near
# canopy top, which gives the greater thermal diffusivity there, and the
# difference fades upward.
.stretchedheight <- function(z, h, d, pai) {
  zp <- z - d
  b <- (h - d) * sqrt(.sublayerfactor(pai)^2 - 1)
  0.5 * (zp + sqrt(zp^2 + b^2))
}
#' Calculate resistance to sensible heat loss of entire canopy
#' @param z height to which resistance is needed (m)
#' @param uf wind friction velocity as returned by [windfric()] (m/s)
#' @param h vegetation height (m)
#' @param pai plant area index
#' @param LL Obukhov length as returned by [Obukhov()] (m). The default, `Inf`,
#' is neutral air.
#' @param MOST Logical. If `TRUE`, profiles above the canopy are taken to be
#' logarithmic down to canopy top. If `FALSE`, a roughness sublayer is assumed
#' above the canopy.
#' @returns bulk surface aerodynamic resistance (s/m)
#' @export
canopyresistance <- function(z, uf, h, pai, LL = Inf, MOST = TRUE) {
  # Bare ground has no canopy to form a roughness sublayer
  if (!MOST && h == 0) {
    warning("No roughness sublayer over bare ground (h = 0): using MOST = TRUE")
    MOST <- TRUE
  }
  d <- zeroplanedis(h, pai)
  zh <- 0.2 * roughlength(h, pai, d) # roughness length for heat
  # Resistance is the integral of one over thermal diffusivity between two heights.
  # Where diffusivity rises in proportion to height above the displacement height,
  # that integral is the logarithm of the ratio of the two heights, with a
  # diabatic correction for heat
  if (MOST) {
    lnz <- log((z - d) / zh) + dpsih(zh / LL) - dpsih((z - d) / LL)
  } else {
    # In a roughness sublayer diffusivity above canopy top is enhanced, which
    # shortens the part of the integral that lies above the canopy
    s1 <- .stretchedheight(h, h, d, pai)
    s2 <- .stretchedheight(z, h, d, pai)
    lnz <- log((h - d) / zh) + dpsih(zh / LL) - dpsih((h - d) / LL) +
      log(s2 / s1) + dpsih(s1 / LL) - dpsih(s2 / LL)
  }
  rHa <- lnz / (0.41 * uf)
  return(rHa)
}
# Speed of the air movements that mix the air next to the ground (m/s). Gusts
# from above the canopy reach the ground weakened by the foliage. Where foliage
# shelters the ground, air warmed by the ground also rises of its own accord,
# which restores part of what the shelter removes and keeps the air next to the
# ground from becoming still in light wind
.floorgust <- function(uf, pai) {
  s <- exp(-4.18 * 0.25 * pai / (1 + 0.4 * pai)) # share of canopy-top gusts that reach the ground
  g <- 0.04                                      # speed of free convection at the ground (m/s)
  restored <- uf^2 * g^2 / (uf^2 + g^2)
  ug <- sqrt((s * uf)^2 + (1 - s^2) * restored)
  return(ug)
}
# Thermal diffusivity at any height between the ground and the reference height.
# `cf` is the factor by which thermal diffusivity at canopy top exceeds what the
# logarithmic law gives there: one for a standard profile
.diffusivityprofile <- function(z, uf, h, pai, LL, cf) {
  ka <- 0.41
  d <- zeroplanedis(h, pai)

  # Above the canopy, thermal diffusivity rises with height above the
  # zero-plane displacement. The eddies shed by a canopy add to it at canopy top,
  # an effect that fades with height
  zp <- pmax(z, h) - d
  b <- (h - d) * sqrt(cf^2 - 1)
  Ka <- ka * uf * sqrt(zp^2 + b^2) / dphih(0.5 * (zp + sqrt(zp^2 + b^2)) / LL)
  Ktop <- ka * uf * cf * (h - d) / dphih(0.5 * (1 + cf) * (h - d) / LL) # at canopy top

  # Within the canopy, heat is carried by gusts that penetrate from above. Their
  # strength falls from canopy top to the ground, by more the denser the foliage
  a1 <- 1.25                                           # gust strength at canopy top, relative to friction velocity
  a0 <- a1 * exp(-4.18 * 0.25 * pai / (1 + 0.4 * pai)) # gust strength at the ground
  gust <- 0.5 * (a1 + a0) - 0.5 * (a1 - a0) * cos(pi * pmin(z, h) / h)
  # Thermal diffusivity is gust strength squared times the time for which a parcel
  # of air keeps its vertical motion, set so that there is no jump at canopy top
  Kc <- Ktop * (gust / a1)^2
  Kc0 <- Ktop * (a0 / a1)^2 # at the ground

  # Close to the ground, eddies can be no larger than their height above it. The
  # air there behaves as it does over bare ground, driven by the gusts that reach
  # it and by free convection
  ug <- .floorgust(uf, pai)
  Kg <- ka * ug * z / dphih(z / LL)

  # Whichever of the two limits mixing sets thermal diffusivity. The ground adds
  # resistance only where it mixes the air less than the canopy gusts at the ground do
  Kcan <- ifelse(z > h, Ka, Kc)
  KH <- 1 / (1 / Kcan + pmax(1 / Kg - 1 / Kc0, 0))
  KH[z <= 0] <- 0
  return(KH)
}
# Resistance to heat transfer between the ground surface and a height above the
# canopy: the integral of one over thermal diffusivity between them
.groundtoheight <- function(z, uf, h, pai, LL, MOST, z0 = 0.004) {
  cf <- 1
  if (!MOST) cf <- .sublayerfactor(pai)
  # The integral starts at the height below which heat leaves the ground by
  # conduction through still air: a fifth of the surface's roughness length or,
  # where little turbulence reaches the ground, the height at which turbulent
  # mixing matches molecular diffusion
  ug <- .floorgust(uf, pai) # speed of the air movements that mix the air next to the ground
  zh <- max(0.2 * z0, 2.1e-5 / (0.41 * ug))
  invK <- function(zz) 1 / .diffusivityprofile(zz, uf, h, pai, LL, cf)
  # Through the canopy, where there is one, and then through the air above it
  rH <- 0
  if (h > zh) rH <- integrate(invK, zh, h, rel.tol = 1e-6, stop.on.error = FALSE)$value
  if (z > h) rH <- rH + integrate(invK, max(h, zh), z, rel.tol = 1e-6, stop.on.error = FALSE)$value
  return(rH)
}
#' Calculate thermal diffusivity at any height between the ground and the reference height
#'
#' Thermal diffusivity measures how readily turbulence carries heat vertically.
#' The resistance to heat transfer between two heights is the integral of one
#' over thermal diffusivity between them. A roughness sublayer is assumed above
#' the canopy, the effect of which fades as foliage becomes sparse.
#' @param z height above ground (m). Can be a vector.
#' @param uf friction velocity (m/s) as returned by [windfric()]
#' @param h vegetation height (m). Zero for bare ground.
#' @param pai plant area index
#' @param LL Obukhov length (m) as returned by [Obukhov()]. Infinite under neutral conditions.
#' @returns thermal diffusivity (m^2/s)
#' @export
thermaldiffusivity <- function(z, uf, h, pai, LL = Inf) {
  # Bare ground: eddies grow in size with height above the surface, so thermal
  # diffusivity rises in proportion to height
  if (h == 0) {
    KH <- 0.41 * uf * z / dphih(z / LL)
    KH[z <= 0] <- 0
    return(KH)
  }
  KH <- .diffusivityprofile(z, uf, h, pai, LL, .sublayerfactor(pai))
  return(KH)
}
#' Calculate resistance to sensible heat loss of individual canopy elements
#' @param tair air temperature (deg C)
#' @param dT difference between leaf and air temperature (deg C), as a positive value
#' @param uz wind speed at height of leaf (m/s)
#' @param len leaf length (m)
#' @param wid leaf width (m)
#' @param x ratio of horizontal to vertical projections of leaf foliage
#' @param rHmax maximum resistance (used to ensure stable convergence (s/m))
#' @returns leaf boundary layer resistance (s/m)
#' @export
leafresistance <- function(tair, dT, uz, len, wid, x, rHmax = 300.0) {
  Tk <- tair + 273.15
  # Properties of air that set how readily it carries heat and how readily it
  # flows: thermal diffusivity and kinematic viscosity (m^2/s)
  Kh <- 1.6667e-10 * Tk^2 + 2.9935e-8 * Tk - 1.7128e-6
  v <- 1.326e-5 * (Tk / 273.15)^2.5 * (393.15 / (Tk + 120))
  # Size of the leaf as a horizontal wind meets it (m): its dimensions scaled by
  # the fraction of its area that faces the wind, which is greatest for upright leaves
  Gx <- 1 / (x + 1.774 * (x + 1.182)^(-0.733))
  d <- sqrt(len * wid) * Gx

  # Heat carried away by the wind (forced convection). The Reynolds number
  # measures the strength of the flow over the leaf, and the Nusselt number how
  # much faster heat is lost than through still air. The leaf is treated as a
  # flat plate lying along the wind, for which the Nusselt number rises with the
  # square root of the Reynolds number
  Re <- (uz * d) / v
  Pr <- v / Kh # Prandtl number: how readily air flows relative to how readily it carries heat
  # Turbulence in outdoor wind raises heat loss to 1.4 times that in smooth flow
  Nuf <- 1.4 * 0.664 * sqrt(Re) * Pr^(1/3)

  # Heat carried away by air that rises or sinks because the leaf warms or cools
  # it (free convection), which matters when there is little wind. The Grashof
  # number measures the strength of that buoyancy
  Gr <- (9.81 * d^3 * dT) / (Tk * v^2)
  Nun <- (0.825 + (0.387 * (Gr * Pr)^(1/6)) /
            (1 + (0.492 / Pr)^(9/16))^(8/27))^2

  # The two act together, the stronger dominating
  Nu <- (Nuf^3 + Nun^3)^(1/3)
  # Resistance of the layer of air around the leaf to the loss of heat
  rHa <- d / (Kh * Nu)
  rHa[rHa > rHmax] <- rHmax
  return (rHa)
}
#' Calculate wind profile below canopy
#' @param hgt height of canopy(m)
#' @param paii vector of plant area index values per unit ground area ordered
#' sequentially from canopy bottom to top.
#' @param uh wind speed at top of canopy as derived using [windprofile_above()] (m/s)
#' @param z0 roughness length of the ground surface (m)
#' @returns wind speed at the middle of each canopy layer (m/s)
#' @export
windprofile_below <- function(hgt, paii, uh, z0 = 0.004) {
  n <- length(paii)
  dz <- hgt / n # depth of each layer (m)
  pai <- sum(paii)
  # Ratio of friction velocity to wind speed at canopy top. It rises with plant
  # area until the canopy is closed, when more foliage adds no more drag
  Be <- min(sqrt(0.003 + 0.15 * pai), 0.3)

  # Wind falls with the foliage passed on the way down from canopy top. What
  # counts is the cube root of foliage density (plant area per unit volume),
  # summed from canopy top to the middle of each layer
  a <- paii / dz
  a13 <- a^(1 / 3)
  cumA <- rev(cumsum(rev(a13 * dz))) - 0.5 * a13 * dz

  # How far the eddies that carry momentum downward reach (the mixing length)
  # depends on the density of the foliage where the wind's momentum is absorbed.
  # That depends on the wind profile itself, so the two are found by iteration.
  # In a canopy of uniform foliage the result is an exponential decline in wind
  cd <- 0.25 # drag coefficient of a leaf
  lmix <- 2 * Be^3 / (cd * pai / hgt)
  for (iter in 1:60) {
    ui <- exp(-(cd / (2 * lmix^2))^(1 / 3) * cumA) # wind as a fraction of that at canopy top
    wgt <- a * ui^2                                # momentum absorbed by each layer
    lnew <- 2 * Be^3 / (cd * sum(a * wgt) / sum(wgt))
    if (abs(lnew - lmix) < 1e-10 * lmix) break
    lmix <- lmix + 0.6 * (lnew - lmix)
  }

  # In the lowest tenth of the canopy the drag of the ground takes over, and wind
  # follows a logarithmic profile down to the ground surface
  n2 <- n %/% 10
  if (n2 >= 1 && 0.5 * dz > z0) {
    zl <- ((1:n2) - 0.5) * dz # height of the middle of each of these layers
    za <- (n2 + 0.5) * dz     # height of the layer above them
    ui[1:n2] <- ui[n2 + 1] * log(zl / z0) / log(za / z0)
  }
  uz <- uh * ui
  return(uz)
}


# Strength of the gusts that carry heat and vapour within a canopy (m/s): the
# typical vertical speed of the air. It is greatest at canopy top and falls
# towards the ground, by more the denser the foliage
.gustspeed <- function(z, uf, h, pai) {
  a1 <- 1.25                                           # at canopy top, relative to friction velocity
  a0 <- a1 * exp(-4.18 * 0.25 * pai / (1 + 0.4 * pai)) # at the ground
  uf * (0.5 * (a1 + a0) - 0.5 * (a1 - a0) * cos(pi * pmin(z, h) / h))
}
# Resistance to heat transfer between the ground surface and each of a set of
# heights within the canopy, given in ascending order, and between the ground
# surface and canopy top: the integral of one over thermal diffusivity
.groundtolayers <- function(z, uf, h, pai, LL) {
  cf <- .sublayerfactor(pai)
  # The integral starts at the height below which heat leaves the ground by
  # conduction through still air
  ug <- .floorgust(uf, pai) # speed of the air movements that mix the air next to the ground
  zh <- max(0.2 * 0.004, 2.1e-5 / (0.41 * ug))
  invK <- function(zz) 1 / .diffusivityprofile(zz, uf, h, pai, LL, cf)
  # Resistance of the air between each height and the next, added up from the ground
  edges <- c(zh, pmax(z, zh), h)
  steps <- numeric(length(edges) - 1)
  for (i in seq_along(steps)) {
    if (edges[i + 1] > edges[i]) {
      steps[i] <- integrate(invK, edges[i], edges[i + 1], rel.tol = 1e-6, stop.on.error = FALSE)$value
    }
  }
  rH <- cumsum(steps)
  return(list(rHz = rH[seq_along(z)], rHh = rH[length(rH)]))
}
# Dilogarithm on (0,1), by its defining series with the reflection identity
# applied on the upper half where that series converges too slowly. Needed to
# integrate the near-field kernel in closed form.
.dilog01 <- function(x) {
  x <- pmin(pmax(x, 0), 1)
  out <- numeric(length(x))
  k <- 1:60
  ser <- function(y) sum(y^k / k^2)
  lo <- x > 0 & x <= 0.5
  hi <- x > 0.5 & x < 1
  if (any(lo)) out[lo] <- vapply(x[lo], ser, numeric(1))
  if (any(hi)) {
    out[hi] <- pi^2 / 6 - log(x[hi]) * log1p(-x[hi]) -
      vapply(1 - x[hi], ser, numeric(1))
  }
  out[x >= 1] <- pi^2 / 6
  out
}
# Integral of the near-field kernel from zero to x. The kernel is the same
# either side of a source, so the integral changes sign with x
.nearfield_integral <- function(x) {
  e <- exp(-abs(x))
  sign(x) * (-0.39894 * (.dilog01(e) - pi^2 / 6) - 0.15623 * (1 - e))
}
# Factor that converts a concentration of heat (J/m^3) into a temperature
# difference (deg C), or a concentration of latent heat into a difference in
# vapour pressure (kPa)
.fromconcentration <- function(tair, pk, vapour) {
  ph <- rhoair(tair, pk) # molar density of air (mol/m^3)
  if (vapour) {
    # Latent heat to moles of vapour, and moles of vapour per mole of air to vapour pressure
    cv <- pk / (latvap(tair) * ph)
  } else {
    cv <- 1 / (ph * cpair(tair))
  }
  return(cv)
}
#' Calculate the far-field part of air temperature or vapour pressure within a canopy
#'
#' The air within a canopy receives heat and water vapour from the foliage of
#' every layer and from the ground. The far field is the part of its temperature
#' or vapour pressure that results from these being mixed gradually through the
#' canopy by turbulence, in the way that heat diffuses: each layer of air
#' exchanges with the layers above and below it across the resistance between
#' them. Adding the near field, as returned by [nearfield()], gives air
#' temperature or vapour pressure.
#' @param Fz vector of the flux of sensible heat from the foliage of each canopy
#' layer, per unit ground area (W/m^2), ordered from canopy bottom to top. Flux
#' of latent heat if `vapour = TRUE`.
#' @param hgt height of canopy (m)
#' @param pai plant area index of the canopy as a whole
#' @param uf friction velocity (m/s) as returned by [windfric()] or [solve_wholecanopy()]
#' @param top air temperature at canopy top (deg C). Vapour pressure of the air
#' at canopy top (kPa) if `vapour = TRUE`.
#' @param ground temperature of the ground surface (deg C). Vapour pressure at
#' the ground surface (kPa) if `vapour = TRUE`: the saturated vapour pressure at
#' ground surface temperature, as returned by [satvap()], times the fractional
#' relative humidity of the soil.
#' @param tair air temperature within the canopy (deg C), from which the density
#' and heat capacity of air are calculated. A single value or one for each
#' height at which values are returned.
#' @param pk atmospheric pressure (kPa)
#' @param LL Obukhov length (m) as returned by [Obukhov()]. Infinite under neutral conditions.
#' @param vapour optional logical indicating whether to calculate vapour
#' pressure (TRUE) or temperature (FALSE)
#' @param z optional vector of heights above ground (m) at which values are
#' returned, between the ground surface (0) and canopy top (`hgt`). The default
#' is the middle of each canopy layer.
#' @returns a vector of the far-field air temperature (deg C) at each height
#' `z`. Vapour pressure (kPa) if `vapour = TRUE`. It equals `ground` at the
#' ground surface and `top` at canopy top.
#' @seealso [nearfield()], [thermaldiffusivity()]
#' @export
farfield <- function(Fz, hgt, pai, uf, top, ground, tair, pk = 101.3, LL = Inf, vapour = FALSE,
                     z = NULL) {
  n <- length(Fz)
  zs <- ((1:n) - 0.5) * hgt / n # height of the middle of each layer, where its heat is released (m)
  if (is.null(z)) z <- zs
  if (any(z < 0 | z > hgt)) stop("z must lie between the ground and canopy top\n")
  # Resistance between the ground surface and each height at which heat is
  # released or a value is wanted, and between the ground surface and canopy top
  zall <- sort(unique(c(zs, z)))
  res <- .groundtolayers(zall, uf, hgt, pai, LL)
  rHs <- res$rHz[match(zs, zall)] # to each layer
  rHz <- res$rHz[match(z, zall)]  # to each height wanted
  rHh <- res$rHh

  # With no foliage, the air at each height lies between the ground and
  # canopy-top values, in proportion to the resistance either side of it
  base <- (rHz * top + (rHh - rHz) * ground) / rHh

  # Heat released by a layer raises the temperature of the air there until it is
  # carried away as fast as it is released, part upward to canopy top and part
  # downward to the ground. The rise is greatest at the source and falls away
  # to nothing at canopy top and at the ground. Element [i, j] is the rise at
  # height i for each W/m^2 released by layer j
  lower <- outer(rHz, rHs, pmin)
  upper <- outer(rHz, rHs, pmax)
  W <- lower * (rHh - upper) / rHh
  Cf <- as.vector(W %*% Fz) # concentration of heat added by the foliage (J/m^3)

  far <- base + Cf * .fromconcentration(tair, pk, vapour)
  return(far)
}
#' Calculate the near-field part of air temperature or vapour pressure within a canopy
#'
#' Heat or vapour released by foliage stays concentrated close to where it is
#' released for a short while, because air moving away from a leaf at first keeps
#' the direction in which it set off. The air next to dense foliage is therefore
#' warmer (or cooler, or more humid) than gradual mixing alone would give. This
#' excess is the near field, calculated following Raupach (1989). Adding it to
#' the far field, as returned by [farfield()], gives air temperature or vapour
#' pressure.
#' @param Fz vector of the flux of sensible heat from the foliage of each canopy
#' layer, per unit ground area (W/m^2), ordered from canopy bottom to top. Flux
#' of latent heat if `vapour = TRUE`.
#' @param hgt height of canopy (m)
#' @param pai plant area index of the canopy as a whole
#' @param uf friction velocity (m/s) as returned by [windfric()] or [solve_wholecanopy()]
#' @param tair air temperature within the canopy (deg C), from which the density
#' and heat capacity of air are calculated. A single value or one for each
#' height at which values are returned.
#' @param pk atmospheric pressure (kPa)
#' @param LL Obukhov length (m) as returned by [Obukhov()]. Infinite under neutral conditions.
#' @param vapour optional logical indicating whether to calculate vapour
#' pressure (TRUE) or temperature (FALSE)
#' @param z optional vector of heights above ground (m) at which values are
#' returned, between the ground surface (0) and canopy top (`hgt`). The default
#' is the middle of each canopy layer.
#' @returns a vector of the near-field part of air temperature (deg C) at each
#' height `z`, relative to the near field at canopy top. Vapour pressure (kPa) if
#' `vapour = TRUE`.
#' @details The near field is given relative to that at canopy top because the
#' air temperature or vapour pressure at canopy top passed to [farfield()] is
#' that of the air itself, which includes any near field there.
#' @references Raupach MR (1989) A practical Lagrangian method for relating
#' scalar concentrations to source distributions in vegetation canopies.
#' Quarterly Journal of the Royal Meteorological Society 115: 609-632.
#' @seealso [farfield()], [thermaldiffusivity()]
#' @export
nearfield <- function(Fz, hgt, pai, uf, tair, pk = 101.3, LL = Inf, vapour = FALSE, z = NULL) {
  n <- length(Fz)
  dz <- hgt / n
  zs <- ((1:n) - 0.5) * dz # height of the middle of each layer, where its heat is released (m)
  if (is.null(z)) z <- zs
  if (any(z < 0 | z > hgt)) stop("z must lie between the ground and canopy top\n")

  # How far the near field of a source reaches (m): the distance air travels
  # while it keeps its vertical motion, which is thermal diffusivity over gust speed
  ow <- .gustspeed(zs, uf, hgt, pai)
  reach <- thermaldiffusivity(zs, uf, hgt, pai, LL) / ow

  # Concentration of heat at a set of heights for each W/m^2 released by each
  # layer. It is high close to a source and falls away over the distance its
  # near field reaches. The source is spread evenly through its layer, and the
  # ground returns what reaches it, as if from a mirror image of the source
  # beneath the ground. Element [i, j] is for height i and layer j
  perunit <- function(zobs) {
    m <- length(zobs)
    zo <- matrix(zobs, m, n)
    lo <- matrix(zs - 0.5 * dz, m, n, byrow = TRUE) # bottom of each source layer
    hi <- matrix(zs + 0.5 * dz, m, n, byrow = TRUE) # top of each source layer
    rc <- matrix(reach, m, n, byrow = TRUE)
    direct <- .nearfield_integral((zo - lo) / rc) - .nearfield_integral((zo - hi) / rc)
    mirror <- .nearfield_integral((zo + hi) / rc) - .nearfield_integral((zo + lo) / rc)
    (direct + mirror) * rc / dz / matrix(ow, m, n, byrow = TRUE)
  }
  Cn <- as.vector(perunit(z) %*% Fz)    # at each height wanted (J/m^3)
  Ch <- as.vector(perunit(hgt) %*% Fz)  # at canopy top

  near <- (Cn - Ch) * .fromconcentration(tair, pk, vapour)
  return(near)
}
#' Calculate air temperature within a canopy
#'
#' Air temperature at any height within a canopy, from the heat released by the
#' foliage of each layer and the temperatures of the air at canopy top and of
#' the ground surface. It is the sum of the far field, as returned by
#' [farfield()], and the near field, as returned by [nearfield()].
#' @param Hz vector of the flux of sensible heat from the foliage of each canopy
#' layer, per unit ground area (W/m^2), ordered from canopy bottom to top
#' @param hgt height of canopy (m)
#' @param pai plant area index of the canopy as a whole
#' @param uf friction velocity (m/s) as returned by [windfric()] or [solve_wholecanopy()]
#' @param ttop air temperature at canopy top (deg C)
#' @param tground temperature of the ground surface (deg C)
#' @param tair estimate of air temperature within the canopy (deg C), from which
#' the density and heat capacity of air are calculated. A single value or one
#' for each height at which values are returned.
#' @param pk atmospheric pressure (kPa)
#' @param LL Obukhov length (m) as returned by [Obukhov()]. Infinite under neutral conditions.
#' @param z optional vector of heights above ground (m) at which values are
#' returned, between the ground surface (0) and canopy top (`hgt`). The default
#' is the middle of each canopy layer.
#' @returns a vector of air temperatures (deg C) at each height `z`
#' @seealso [relhum_below()] for humidity
#' @export
temp_below <- function(Hz, hgt, pai, uf, ttop, tground, tair, pk = 101.3, LL = Inf, z = NULL) {
  # Heat mixed gradually through the canopy, between the ground and canopy top
  far <- farfield(Hz, hgt, pai, uf, ttop, tground, tair, pk, LL, z = z)
  # Heat still concentrated close to the foliage that released it
  near <- nearfield(Hz, hgt, pai, uf, tair, pk, LL, z = z)
  return(far + near)
}
#' Calculate resistance from ground surface to some height zref above canopy
#' @param hgt height of canopy (m)
#' @param pai plant area index
#' @param uf wind friction velocity as returned by [windfric()] (m/S)
#' @param LL Obukhov length as returned by [Obukhov()]
#' @param psih diabatic influencing factor for heat as returned by [dpsih()]
#' @param zref height to which resistance is calculated (m)
#' @returns resistance to heat loss (s/m)
#' @export
groundresistance <- function(hgt, pai, uf, LL, psih, zref) {
  if (zref < hgt) stop("zref must be greater than or equal to hgt\n")
  # Calculate a2
  d <- zeroplanedis(hgt, pai)
  phih <- dphih((zref - d) / LL)
  a2 <- (phih * 0.41 * (1 - d / hgt)) / 1.5625
  # Resistance from ground to canopy top
  rHc <- 4.293251 / (a2 * uf)
  if (zref > hgt) {
    # resistance from zref to hes
    zm <- roughlength(hgt, pai, d)
    rHz <- (log((zref - d) / zm) + psih) / (0.41 * uf)
    # resistance from hgt to hes
    psihh <- dpsih(zm / LL) - dpsih((hgt - d) / LL)
    rHh <- (log((hgt - d) / zm) + psihh) / (0.41 * uf)
    rHa <- rHz - rHh
  } else {
    rHa <- 0
  }
  rHg <- rHc + rHa
  return(rHg)
}


