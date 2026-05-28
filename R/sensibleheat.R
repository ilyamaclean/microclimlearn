#' Calculate molar density of air
#'
#' @param tc air temperature (degrees C)
#' @param pk atmospheric pressure (kPa)
#' @returns molar density of air (kPa).
#' @export
rhohair <- function(tc, pk) {
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
  pai[pai < 0.001] <- 0.001
  d <- (1 - (1 - exp(-sqrt(7.5 * pai))) / sqrt(7.5 * pai)) * h
  return(d)
}

#' Calculate roughness length for momentum
#' @param h vegetation height (m)
#' @param pai plant area index
#' @param psi_h diabatic correction for heat
#' @returns roughness length (m)
#' @export
roughlength <- function(h, pai, d, psih) {
  Be   <- sqrt(0.003 + 0.1 * pai)
  zm   <- (h - d) * exp(0.41*psih - 0.4/Be)
  zmlim<- 0.9 * (h - d)
  zm   <- pmin(pmax(zm, 5e-4), zmlim)
  return (zm)
}

#' Calculate integrated diabatic correction coefficient for momentum
#' @param ze atmospheric stability parameter
#' @returns diabatic correction coefficient for momentum
#' @export
dpsim <- function(ze) {
  out <- rep(NA_real_, length(ze))
  zneg <- ze < 0L
  # unstable branch
  if (any(zneg)) {
    x  <- pmax(0, 1 - 15 * ze[zneg])^(0.25)     # guard against round-off
    out[zneg] <- log(((1 + x) / 2)^2 * (1 + x^2) / 2) -
      2 * atan(x) + pi/2
  }
  # stable branch
  if (any(!zneg)) {
    out[!zneg] <- -4.7 * ze[!zneg]
  }
  # clamp to [-4, 3]
  pmax(pmin(out, 3), -4)
}
#' Calculate integrated diabatic correction coefficient for heat
#' @param ze atmospheric stability parameter
#' @returns diabatic correction coefficient for heat
#' @export
dpsih <- function(ze) {
  out <- rep(NA_real_, length(ze))
  zneg <- ze < 0
  # unstable branch
  y <- sqrt(1 - 9 * ze[zneg])          # domain OK because ze[zneg] < 0
  out[zneg] <- log(((1 + y) / 2)^2)
  # stable branch
  out[!zneg] <- -(4.7 * ze[!zneg]) / 0.74
  # clamp to [-4, 3]
  pmax(pmin(out, 3), -4)
}
#' Calculate diabatic influencing coefficient for heat
#' @param ze atmospheric stability parameter
#' @returns diabatic influencing coefficient for heat
dphih <- function(ze) {
  phih <- rep(NA_real_, length(ze))
  zneg <- ze < 0
  # unstable branch
  phim <- 1 / ((1.0 - 16.0 * ze[zneg])^0.25)
  phih[zneg] <- phim^2
  # stable branch
  phih[!zneg] <- 1 + ((6 * ze[!zneg]) / (1 + ze[!zneg]))
  phih[phih > 1.5] <- 1.5
  phih[phih < 0.5] <- 0.5
  return (phih)
}
#' Calculate Obukhov length
#' @param tc air temperature (degrees C)
#' @param pk atmospheric pressure (kPa)
#' @param uf wind friction velocity (m/s)
#' @param H sensible heat flux (W/m^2)
#' @returns Obukhov length (m)
#' @export
Obukhov <- function(tc, pk, uf, H) {
  ph <- rhohair(tc, pk)
  cp <- cpair(tc)
  Tk <- tc + 273.15
  LL <- (ph * cp * uf^3 * Tk) / (-0.41 * 9.81 * H)
  return (LL)
}
#' Calculate wind friction velocity
#' @param uref wind speed at reference height (m/s)
#' @param zref height above ground of wind speed measurement (m)
#' @param d zero-plane displacement height as returned by [zeroplanedis()] (m)
#' @param zm roughness length for momentum as returned by [roughlength()] (m)
#' @param psim diabatic correction for momentum as returned by [dpsim()]
#' @returns friction velocity (m/s)
#' @export
windfric <- function(uref, zref, d, zm, psim) {
   uf <- (0.41 * uref)  / (log((zref - d) / zm) + psim)
   return(uf)
}
#' Calculate wind profile above canopy
#' @param z height above ground for which wind speed desired (m)
#' @param uf wind friction velocity as returned by [windfric()] (m/s)
#' @param d zero-plane displacement height as returned by [zeroplanedis()] (m)
#' @param zm roughness length for momentum as returned by [roughlength()] (m)
#' @param psim diabatic correction for momentum as returned by [dpsim()]
#' @returns wind speed (m/s)
#' @export
windprofile_above <- function(z, uf, d, zm, psim) {
  uz <- (uf / 0.41) * (log((z - d) / zm) + psim)
}
#' Calculate resistance to sensible heat loss of entire canopy
#' @param z height to which resistance is needed (m)
#' @param d zero-plane displacement height as returned by [zeroplanedis()] (m)
#' @param zm roughness length for momentum as returned by [roughlength()] (m)
#' @param uf wind friction velocity as returned by [windfric()] (m/s)
#' @param psih diabatic correction for heat as returned by [dpsih()]
#' @returns bulk surface aerodynamic resistance (m/s)
#' @export
canopyresistance <- function(z, d, zm, uf, psih) {
  zh <- 0.2 * zm # roughness length for heat
  rHa <- (log((z - d) / zh) + psih) / (0.41 * uf)
  return(rHa)
}
#' Calculate resistance to sensible heat loss of individual canopy elements
#' @param tair air temperature (deg C)
#' @param dT difference between leaf and air temperature (deg C)
#' @param uz wind speed at height of leaf (m/S)
#' @param len leaf length (m)
#' @param wid leaf width (m)
#' @param x ratio of vertical to horizontal projections of leaf foliage
#' @param rHmax maximum resistance (used to ensure stable convergence (s/m))
#' @returns leaf boundary layer resistance
#' @export
leafresistance <- function(tair, dT, uz, len, wid, x, rHmax = 300.0) {
  Tk <- tair + 273.15
  # Compute thermal diffusivity
  Kh <- 1.6667e-10 * Tk^2 + 2.9935e-8 * Tk - 1.7128e-6
  # Compute kinematic viscosity
  v <- 1.326e-5 * (Tk / 273.15)^1.5 * (393.55 / (Tk + 120))
  # Compute projected area in direction of wind flow (horizontal)
  y <- 1 / x
  Gy <- sqrt(y^2 + (1 - y^2)) / (y + 1.774 * (y + 1.182)^(-0.733))
  # Compute characteristic dimension
  d <- sqrt(len * wid) * Gy
  # Compute Reynolds number
  Re <- (uz * d) / v
  # Compute Prandlt number
  Pr <- v / Kh;
  # Compute Nusselt number for forced convection
  Nuf <- ifelse(Re > 2e5, 0.37 * Re^0.6 * Pr^(1/3) + 9.08,
    ifelse(Re > 1000, 2 + (0.48 * Re^0.6 - 11.31) * Pr^(1/3),
      2 + 0.6 * sqrt(Re) * Pr^(1/3)))
  # Compute Grashof number
  Gr <- (9.81 * d^3 * dT) / (Tk * v^2)
  # Compute Nusselt number for free convection
  Nun <- (0.825 + (0.387 * (Gr * Pr)^(1/6)) /
            (1 + (0.492 / Pr)^(9/16))^(8/27))^2
  # Compute Nusselt number (mixed forced and free)
  Nu <- (Nuf^3 + Nun^3)^(1/3)
  # Compute boundary layer resistance
  rHa <- d / (Kh * Nu)
  rHa[rHa > rHmax] <- rHmax
  return (rHa)
}
#' Calculate resistance from ground surface to some height zref above canopy
#' @param hgt height of canopy (m)
#' @param pai plant area index
#' @param uf wind friction velocity as returned by [windfric()] (m/S)
#' @param LL Obukhov length as returned by [Obukhov()]
#' @param psih diabatic influencing factor for heat as returned by [dpsih()]
#' @param zref height to which resistance is calculated (m)
#' @returns resistance to heat loss (s/m)
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
    zm <- roughlength(hgt, pai, d, psih)
    rHz <- (log((zref - d) / zm) + psih) / (0.41 * uf)
    # resistance from hgt to hes
    psihh <- dpsih(zm / LL) - dpsim((hgt - d) / LL)
    rHh <- (log((hgt - d) / zm) + psihh) / (0.41 * uf)
    rHa <- rHz - rHh
  } else {
    rHa <- 0
  }
  rHg <- rHc + rHa
  return(rHg)
}
#' Calculate wind profile above canopy
#' @param hgt height of canopy(m)
#' @param paii vector of plant area index values per unit ground area ordered
#' sequentially from canopy bottom to top.
#' @param uh wind speed at top of canopy as derived using [windprofile_above()] (m/s)
#' @returns wind speed (m/s)
#' @export

windprofile_below <- function(hgt, paii, uh) {
  # Calculate whole canopy attenuation coefficient
  pai <- sum(paii)
  Be <- sqrt(0.003 + 0.1 * pai)
  a <- pai / hgt
  Lc <- (0.25 * a)^(-1)
  Lm <-  2 * Be^3 * Lc
  att <- Be * hgt / Lm
  # Calculate attenuation coefficient for canopy elements
  Bei <- sqrt(0.003 + 0.1 * paii)
  ai <- paii / hgt
  Lci <- (0.25 * ai)^(-1)
  Lmi <-  2 * Bei^3 * Lci
  ati <- Bei * hgt / Lmi
  ati <- (ati / sum(ati)) * att
  # Calculate scaled profile
  n <- length(paii)
  n2 <- trunc(n / 10)
  ui <- rep(1, n)
  # From canopy top down to the lower 10%
  ui[n2:n] <- c(rev(cumprod(rev(1 - ati[n2:(n - 1)]))), 1)
  # Bottom 10%
  zm <- hgt / (20 * n2)
  z2 <- (1:n2) * hgt / (10 * n2)
  uf <- (0.41 * ui[n2]) / log(hgt / (10 * zm))
  ui[1:n2] <- (uf / 0.41) * log(z2 / zm)
  uz <- uh * ui
  return(uz)
}
#' Calculate air temperature below canopy
#' @param hgt height of canopy(m)
#' @param paii vector of plant area index values per unit ground area ordered
#' sequentially from canopy bottom to top.
#' @param TL Langrangian timescale as returned by [solve_wholecanopy()]
#' @param uf wind friction velocity as returned by [solve_wholecanopy()]
#' @param pk atmospheric pressure (kPa)
#' @param Hz vector of sensible heat fluxes from each canopy layer (W/m^2)
#' @param tair vector of initial estimates of air temperatures (deg C).
#' @param tleaf vector of leaf temperatures (deg C)
#' @param tground ground surface temperature (deg C)
#' @returns a vector of air temperatures (deg C)
#' @export
temp_below <- function(hgt, paii, TL, uf, pk, Hz, tair, tleaf, tground) {
  # Generate z vector
  n <- length(paii)
  z <- (c(1:n) / n) * hgt
  dz <- hgt / n
  # extract canopy top temperature
  th <- tair[n]
  tmn <- min(c(tground, th, tleaf)) - 2
  tmx <- max(c(tground, th, tleaf)) + 2
  # Calculate position variance
  ow <- uf * (0.75 + 0.5 * cos(pi * (1 - z / hgt)))
  # Calculate thermal diffusivity
  KH <- TL * ow^2
  SS <- paii * Hz
  # Compute near-field correction factor for small sample size
  n <- length(paii)
  mu <- 1 + 0.894 * exp(-0.01386 * n) + 9.82 * exp(-0.15 * n)
  # Compute near-field concentrations at the top of the canopy
  dz1 <- (hgt - z[-n]) / (ow[-n] * TL)
  kn  <- -0.39894 * log1p(-exp(-dz1)) - 0.15623 * exp(-dz1)
  Cnh <- sum((SS[-n] / ow[-n]) * kn * (2 * hgt / (ow[-n] * TL))) * mu
  # Compute far field concentration at top of canopy
  Cfh <- rhohair(th, pk) * cpair(th) * th
  # Compute near and far-field concentrations for each canopy element
  sumRH <- 0
  for (i in 1:n) {
    # Compute resistance from ground to z
    RH <- 1 / KH[i]
    sumRH <- sumRH + RH
    rHa = sumRH * dz;
    if (rHa < 2) rHa = 2.0;
    ph <- rhohair(tair[i], pk)
    cp <- cpair(tair[i])
    GT <- (ph * cp / rHa) * (tground - tair[i]) * dz;
    H <- sum(SS[1:i]) + GT
    # Compute far-field concentration
    Cf <- sum((H / KH[i:n]) * dz)
    # Compute near-field concentration
    dz1 <- (z[i] - z) / (ow * TL)
    dz2 <- (z[i] + z) / (ow * TL)
    zeta <- abs(dz1)
    e <- exp(-zeta)
    kn <- -0.39894 * log1p(-e) - 0.15623 * e
    Cn <- sum((SS[-i] / ow[-i]) * kn[-i] * (dz1[-i] + dz2[-i]))
    # Compute total source concentration
    CT <- Cfh - Cnh + Cf + Cn * mu
    tair[i] <- CT / (ph * cp)
  }
  tair <- pmax(tair, tmn)
  tair <- pmin(tair, tmx)
  tair[n] <- th
  return(tair)
}
#' Weighted Aitken relaxation for profile convergence
#'
#' @details Applies a weighted Aitken acceleration step to update a profile (`newv`)
#' towards convergence with respect to a previous estimate (`oldv`). The method
#' adaptively updates a relaxation parameter (`omega`) based on successive
#' residuals, with vertical weighting that emphasises lower canopy layers. The
#' function operates on a single list (`aitk`) containing both the profile
#' variables and the persistent Aitken state, and returns the updated list. This
#' design supports iterative use without external state management.
#'
#' @param aitk A list containing:
#' \describe{
#'   \item{oldv}{Numeric vector of previous values.}
#'   \item{newv}{Numeric vector of current (unrelaxed) values.}
#'   \item{z}{Numeric vector of heights corresponding to each element.}
#'   \item{hgt}{Scalar canopy height used to normalise \code{z}.}
#'   \item{state}{(Optional) List containing Aitken state:
#'     \describe{
#'       \item{r_prev}{Residuals from the previous iteration.}
#'       \item{omega}{Current relaxation parameter.}
#'       \item{have_prev}{Logical indicating whether a previous iteration exists.}
#'     }
#'   }
#' }
#' @param omega_min Minimum allowed relaxation parameter.
#' @param omega_max Maximum allowed relaxation parameter.
#' @param w_bot Weight applied at the canopy base.
#' @param w_top Weight applied at the canopy top.
#'
#' @return The updated \code{aitk} list with:
#' \describe{
#'   \item{newv}{Relaxed profile values.}
#'   \item{state}{Updated Aitken state for reuse in subsequent iterations.}
#' }
#'
#' @details
#' On the first iteration, a fixed relaxation parameter is applied. On subsequent
#' iterations, \code{omega} is updated using a weighted Aitken scheme based on
#' changes in residuals between iterations. Vertical weights vary smoothly with
#' height, increasing from \code{w_bot} to \code{w_top}. The method is intended for iterative solvers where \code{newv} is repeatedly
#' updated and relaxed until convergence.
#' @export
aitken_weightdif <- function(aitk, omega_min = 0.02, omega_max = 0.90,
    w_bot = 0.05, w_top = 0.80) {
    # Required fields
    req <- c("oldv", "newv", "z", "hgt")
    miss <- setdiff(req, names(aitk))
  if (length(miss) > 0) {
    stop("Input list is missing: ", paste(miss, collapse = ", "))
  }
  oldv <- aitk$oldv
  newv <- aitk$newv
  z    <- aitk$z
  hgt  <- aitk$hgt

  if (length(oldv) != length(newv) || length(oldv) != length(z)) {
    stop("oldv, newv, and z must have the same length.")
  }
  if (!is.numeric(hgt) || length(hgt) != 1 || hgt == 0) {
    stop("hgt must be a single non-zero numeric value.")
  }
  # Initialise nested state if needed
  if (is.null(aitk$state)) {
    aitk$state <- list(
      r_prev = NULL,
      omega = 0.3,
      have_prev = FALSE
    )
  }
  n <- length(oldv)
  if (is.null(aitk$state$r_prev) || length(aitk$state$r_prev) != n) {
    aitk$state$r_prev <- numeric(n)
    aitk$state$have_prev <- FALSE
  }
  if (is.null(aitk$state$omega)) {
    aitk$state$omega <- 0.3
  }
  if (is.null(aitk$state$have_prev)) {
    aitk$state$have_prev <- FALSE
  }
  aitk$state$omega <- max(omega_min, min(aitk$state$omega, omega_max))
  dw <- w_top - w_bot
  # ---- First iteration ----
  if (!isTRUE(aitk$state$have_prev)) {
    r <- newv - oldv
    aitk$state$r_prev <- r
    s <- z / hgt
    wz <- w_bot + dw * s^2
    aitk$newv <- oldv + (aitk$state$omega * wz) * r
    aitk$state$have_prev <- TRUE
    return(aitk)
  }
  # ---- Learn omega ----
  beta <- 10.0
  r  <- newv - oldv
  dr <- r - aitk$state$r_prev
  s <- z / hgt
  t <- 1.0 - s
  g <- 1.0 + beta * t^2
  num <- sum(g * aitk$state$r_prev * dr)
  den <- sum(g * dr^2)
  omega <- aitk$state$omega
  if (den > 0) {
    omega <- -aitk$state$omega * (num / den)
  }
  omega <- max(omega_min, min(omega, omega_max))
  # ---- Apply update ----
  aitk$state$r_prev <- r
  wz <- w_bot + dw * s^2
  aitk$newv <- oldv + (omega * wz) * r
  aitk$state$omega <- omega
  aitk
}


