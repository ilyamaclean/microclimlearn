#' Calculate saturated vapour pressure
#'
#' @param tc temperature (degrees C)
#' @returns saturated vapour pressure (kPa).
#' @export
satvap <- function(tc) {
  neg <- tc < 0L
  es <- rep(NA_real_, length(tc))
  es[!neg] <- 0.61078 * exp(17.27 * tc[!neg] / (tc[!neg] + 237.3))
  es[neg] <- 0.61078 * exp(21.875 * tc[neg] / (tc[neg] + 265.5))
  es
}
#' Calculate latent heat of vapourisation (J / mol)
#'
#' @param Ts surface temperature (degrees C)
#' @returns latent heat of vapourisation (J / mol)
#' @export
latvap <- function(Ts) {
  la <- 45068.7 - 42.8428 * Ts
  la[Ts < 0] <- 51078.69 - 4.338 * Ts[Ts < 0] - 0.06367 * Ts[Ts < 0]^2
  return(la)
}
#' Creates a list of vegetation parameter inputs for a given plant function type
#' @param PFT a plant function type as listed in `PFTtable`
#' @returns a list of vegetation parameters.
#' @export
createplant_inputs <- function(PFT) {

  s  <- match(PFT, names(PFTtable))
  if (is.na(s)) stop("PFT not found in PFTtable")

  v  <- PFTtable[, s]
  mu <- PFTtable$multiplier

  plant_inputs <- c(
    h       = v[1]  * mu[1],   # plant height (m)
    pai     = v[2]  * mu[2],   # Plant area index
    LAIfrac = v[3]  * mu[3],   # Fraction of PAI that is LAI
    x       = v[4]  * mu[4],   # Campbell x parameter
    lref    = v[6]  * mu[6],   # leaf reflectance (PAR)
    ltra    = v[7]  * mu[7],   # leaf transmittance (PAR)
    lrefp   = v[8]  * mu[8],   # leaf reflectance (PAR)
    ltrap   = v[9]  * mu[9],   # leaf transmittance (PAR)
    Vcmax25 = v[10] * mu[10],  # max Rubisco carboxylation @ 25 °C
    Tup     = v[11] * mu[11],  # high-T photosynthesis limit (°C)
    Tlw     = v[12] * mu[12],  # low-T photosynthesis limit (°C)
    Dcrit   = v[13] * mu[13],  # onset of strong stomatal limitation (kPa)
    alpha   = v[14] * mu[14],  # quantum efficiency (mol CO2 / mol IPAR)
    Kxmx    = v[15] * mu[15],  # max xylem hydraulic conductance
    hv      = v[16] * mu[16],  # Huber value
    f0      = v[17] * mu[17],  # ci / (ci − Gma) under low D
    fd      = v[18] * mu[18],  # fraction of Vcmax to Rd
    psi50   = v[19] * mu[19],  # water potential at 50% K loss (MPa)
    apsi    = v[20] * mu[20],  # xylem conductance param
    len     = v[21] * mu[21],  # leaf length (m)
    wid     = v[22] * mu[22],  # leaf width (m)
    pTAW    = v[25] * mu[25],  # water stress coefficient
    rootskew = v[26] * mu[26], # skew towards top of soil profile of roots
    rpmin = v[27] * mu[27], # Minimum resistance
    em = 0.97
  )
  plant_inputs <- as.list(plant_inputs)
  plant_inputs
}
#' Internal plant photosynthesis model
.photosynth_model <-function(Ca, Rswabs, tair, tleaf, rh, pk, psi_r, plant_inputs, C3 = TRUE) {
  # Extract inputs: plant data
  h <- plant_inputs$h
  x <- plant_inputs$x
  lrefp <- plant_inputs$lrefp
  ltrap <- plant_inputs$ltrap
  Vcmax25 <- plant_inputs$Vcmax25
  Tup <- plant_inputs$Tup
  Tlw <- plant_inputs$Tlw
  Dcrit <- plant_inputs$Dcrit
  alpha <- plant_inputs$alpha
  Kxmx <- plant_inputs$Kxmx
  hv <- plant_inputs$hv
  f0 <- plant_inputs$f0
  fd <- plant_inputs$fd
  om <- lrefp + ltrap # leaf scattering coefficient for PAR
  IPAR <- Rswabs * (0.48  / 0.219) * 1e-6 # PAR converstion to (mol photons / m^2 / s)
  # Vcmax temperature response
  Q10leaf = 2
  Vcmax = Vcmax25 * (Q10leaf ^ (0.1 * (tleaf - 25))) / ((1 + exp(0.3 * (tleaf - Tup))) * (1 + exp(0.3 * (Tlw - tleaf))))
  Rd <- fd * Vcmax
  # Photocompensation
  Q10rs = 0.57
  Oa <- 0.2095 * pk * 1000  # Pa
  photocomp <- Oa / (2 * 2600 * (Q10rs ^ (0.1 * (tleaf - 25))))  # Photo compensation point (Pa)
  # Carbon conversion
  ea <- satvap(tair) * (rh / 100)
  es <- satvap(tleaf)
  DD <- es - ea
  ca <- Ca * pk / 1000  # convert to Pa
  ci <- f0 * (1 - DD / Dcrit) * (ca - photocomp) # in Pa
  # C3 photosynthetic pathway
  if (C3) {
    # Kc & Ko
    Q10Kc <- 2.1
    Kc <- 30 * (Q10Kc ^ (0.1 * (tair - 25)))
    Q10Ko = 1.2
    Ko <- 30000 * (Q10Ko ^ (0.1 * (tair - 25)))
    # Assimilation
    om <- lrefp + ltrap
    Wc = Vcmax * ((ci - photocomp) / (ci + Kc * (1 + Oa / Ko)))
    #Wl = alpha * (1 - om) * IPAR * ((ci - photocomp) / (ci + 2 * photocomp))
    Wl = alpha * IPAR * ((ci - photocomp) / (ci + 2.0 * photocomp))
    if (length(Vcmax) == 1) {
      We <- rep(0.5 * Vcmax, length(Wc))
    } else We <- 0.5 * Vcmax
    Wc[Wc < 0] <- 0
    Wl[Wl < 0] <- 0
    We[We < 0] <- 0
    # Calculate rate limited assimilation
    W <- pmin(Wc, Wl, We)
    A <- W - Rd
    # Calculate co-limiting assimilation
    Wcol <- ((We + Wl) - sqrt((We + Wl)^2 - 4 * 0.93 * (We * Wl))) / (2 * 0.93)  # Point at which no longer limited by ci
    Acol <- Wcol - Rd
    cicol <- (-Vcmax * photocomp - Kc * (1 + Oa / Ko) * Wcol) / (Wcol - Vcmax)
  } else {
    k <- 2e-4
    Wc <- rep(Vcmax, length(Ca))
    Wl <- rep(alpha * IPAR, length(Ca))
    We <- k * Vcmax * (ci / (pk * 1000))
    W <- pmin(Wc, Wl, We)
    A <- W - Rd
    Wcol <- ((Wc + Wl) + sqrt((Wc + Wl)^2 - 4 * Be * (Wc * Wl))) / (2 * Be)
    Acol <- Wcol - Rd
    cicol <- (-Vcmax * photocomp - Kc * (1 + Oa / Ko) * Wcol) / (Wcol - Vcmax)
  }
  return(list(A = A, ci = ci, Acol = Acol, cicol = cicol))
}
#' internal function for calculating minimum resistance
.rpmin_calc <- function(h, hv = 0.002, Kxmx) {
  n_ext <- 2 # number of daughter branches produced by each parent branch
  Lpet <- 0.04 # petiole length
  r_intpet <- 10 # petiole conduit radius (micro m)
  r_intref <- 22 # radius of conduits of the terminal branches
  Nh <- (3 * log(1 - (h/Lpet) * (1 - n_ext^(1/3)))) / log(n_ext)  # Number of branching levels
  Nh[Nh < 1] <- 1
  Chi_1 <- 2.888503
  Chi_h <- (6.6e-13 * Nh^1.85) / (7.2e-13 * Nh^1.32)
  Chi_tap <- Chi_h / Chi_1 # tapering factor
  Kpmx <- Kxmx * (r_intpet / r_intref)^2  # maximum petiole level hydraulic conductivity (mol / m / s / MPa)
  rpmin <- h / (Kpmx * hv  * Chi_tap) # minimum m^2 s MPa mol^-1 H20
  return(rpmin)
}
#' internal function for calculating psi50
.psi50toa <- function(psi50){
  stem_Slope <- 65.15*(-psi50)^(-1.25)
  a <- -4*stem_Slope/100*psi50
  return(a)
}
#' Calculate leaf stomatal conductance
#' @param Ca atmospheric carbon dioxide concentration (ppm)
#' @param Rswabs flux density of PAR absorbed by leaf (W/m^2)
#' @param tair temperature of air (deg C)
#' @param tleaf temperature of leaf (deg C)
#' @param rh relative humidity of air (percentage)
#' @param pk atmospheric pressure (kPa)
#' @param psi_r # mean water potential in root zone (MPa)
#' @param plant_inputs a vector of vegetation parameters as returned by [createplant_inputs()]
#' @param z height above ground
#' @param C3 optional logical indicating whether vegetation has C3 or C4 photosynthetic pathway
#' @returns stomatal conductance (mol / m^2 / s)
#' @export
stomatalcond_calc <- function(Ca, Rswabs, tair, tleaf, rh, pk, psi_r, plant_inputs, z, C3 = TRUE) {
  # Extract inputs:
  h <- plant_inputs$h
  Kxmx <- plant_inputs$Kxmx
  hv <- plant_inputs$hv
  psi_50 <- plant_inputs$psi50
  apsi <- plant_inputs$apsi
  # Run photosynthesis model
  photomodel <- .photosynth_model(Ca, Rswabs, tair, tleaf, rh, pk, psi_r, plant_inputs)
  # Compute change in assimilation per ci gradient
  dadc <- with(photomodel, (A - Acol) / (ci - cicol))
  dadc[dadc < 1e-99] <- 1e-99
  # Compute change in conductivity per psi gradient
  if (apsi < 0) apsi <- .psi50toa(psi_50)
  rhow <- 1000 * (1 - ((tleaf + 288.9414) * (tleaf - 3.9863)^2) / (508929.2 * (tleaf + 68.12963)))
  psi_pd <- psi_r - z * 9.81 * rhow * 1e-6
  K_psi_pd <- 1 / (1 + (psi_pd / psi_50)^apsi)
  K_50f <- 0.5 / (1 + ((psi_pd + psi_50) / psi_50)^apsi)
  psi_50f <- (psi_pd + psi_50) / 2
  dKdpKi <- ((K_psi_pd - K_50f) / (psi_pd - psi_50f)) * (1 / K_psi_pd)
  # Compute rp
  rpmin <- .rpmin_calc(h, hv, Kxmx) # mol^-1 m^2 s MPa
  rp <- rpmin / K_psi_pd   # This term appears in equation S1.6
  # Compute zeta
  ea <- satvap(tair) * (rh / 100)
  es <- satvap(tleaf)
  DD <- (es - ea) / pk # divide by pk to convert to mol / mol
  #DD <- pmax(DD, 1e-5)
  DD[DD < 0.05] <- 0.05
  zeta <- 2 / (dKdpKi * rp * 1.6 * DD)
  mu <- 1 + (4 * zeta) / dadc
  mu[mu < 1] <- 1
  # Compute stomatal conductance (reproduce stomatal closure at high ca and low light)
  gs <- 0.5 * dadc * (sqrt(mu) - 1)
  gs <- as.numeric(gs)
  gs[Rswabs == 0] <- 0
  return(gs)
}
#' Calculate bulk surface stomatal resistance
#' @param solp a list of solar positions as returned by [sunposition()]
#' @param PAI total one-sided plant area per unit ground area (m^2 / m^2)
#' @param Ca atmospheric carbon dioxide concentration (ppm)
#' @param Rsw flux density of shortwave radiation on the horizontal (W/m^2)
#' @param Rdif flux density of diffuse radiation on the horizontal (W/m^2)
#' @param tair temperature of air (deg C)
#' @param tleaf temperature of leaf (deg C)
#' @param rh relative humidity of air (percentage)
#' @param pk atmospheric pressure (kPa)
#' @param psi_r # mean water potential in root zone (MPa)
#' @param plant_inputs a vector of vegetation parameters as returned by [createplant_inputs()]
#' @param C3 optional logical indicating whether vegetation has C3 or C4 photosynthetic pathway
#' @returns bul surface stomatal resistance (s/m)
#' @export
bulkstomatalresist_calc <- function(solp, Ca, Rsw, Rdif, tair, tcanopy, rh, pk, psi_r, plant_inputs, C3 = TRUE) {
  # Extract parameters form plant inputs
  z <- plant_inputs$h / 2
  x <- plant_inputs$x
  lrefp <- plant_inputs$lrefp
  ltrap <- plant_inputs$ltrap
  LAIfrac <- plant_inputs$LAIfrac
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
  L_sun <- (1 - exp(-k* PAI * LAIfrac)) / k
  L_shade <- (PAI * LAIfrac) - L_sun
  # Absorbed - sun and shaded parts
  Rshade_abs <- Rdif * ((1 - exp(-PAI)) / PAI) * (1 - om)
  Rsun_abs <- (Rsw - Rdif) * k  * (1 - om) + Rshade_abs
  # Calculate stomatal conductance
  gs_sun <- stomatalcond_calc(Ca, Rsun_abs, tair, tcanopy, rh, pk, psi_r, plant_inputs, z, C3)
  gs_shade <- stomatalcond_calc(Ca, Rshade_abs, tair, tcanopy, rh, pk, psi_r, plant_inputs, z, C3)
  Gs <- gs_sun * L_sun + gs_shade * L_shade
  # Convert to bulk surface resistance
  ph <- rhohair(tair, pk)
  rS <- ph / Gs
  return(rS)
}
#' Dew-point temperature
#'
#' Calculate dew-point temperature (°C) from air temperature (°C) and
#' relative humidity (%). Uses Magnus equations with coefficients for
#' water and ice, applying the ice formulation when dew point is below 0 °C.
#'
#' @param tc Air temperature (°C).
#' @param rh Relative humidity (%).
#'
#' @return Dew-point temperature (°C).
#'
#' @seealso [satvap()]
#' @export
dewpoint <- function(tc, rh) {
  # actual vapour pressure (kPa)
  ea <- satvap(tc) * rh / 100
  # Magnus constants
  a_w <- 17.27
  b_w <- 237.7    # water
  a_i <- 21.875
  b_i <- 265.5    # ice
  # first guess assuming water
  gamma <- log(ea / 0.6108)
  td <- (b_w * gamma) / (a_w - gamma)
  # correct for ice if dewpoint < 0
  ice <- td < 0
  if (any(ice)) {
    gamma_i <- log(ea[ice] / 0.6108)
    td[ice] <- (b_i * gamma_i) / (a_i - gamma_i)
  }
  return(td)
}
#' Calculate relative humidity below canopy
#' @param hgt height of canopy(m)
#' @param paii vector of plant area index values per unit ground area ordered
#' sequentially from canopy bottom to top.
#' @param TL Langrangian timescale as returned by [solve_wholecanopy()]
#' @param uf wind friction velocity as returned by [windfric()] (m/s)
#' @param pk atmospheric pressure (kPa)
#' @param Lz vector of latent heat fluxes from each canopy layer (W/m^2)
#' @param rh vector of initial estimates of relative humidity (percentage)
#' @param tair vector of air temperatures (deg C)
#' @param tleaf vector of leaf temperatures (deg C)
#' @param tground ground surface temperature (deg C)
#' @param soilrh soil effective relative humidity (percentage)
#' @returns a vector of air temperatures (deg C)
#' @export
relhum_below <- function(hgt, paii, TL, uf, pk, Lz, rh, tair, tleaf, tground, soilrh) {
  # Generate z vector
  n <- length(paii)
  z <- (c(1:n) / n) * hgt
  dz <- hgt / n
  # Calculate Langrangian timescale
  pai <- sum(paii)
  d <- zeroplanedis(hgt, pai)
  a2 <- (0.41 * (1 - d / hgt)) / 1.5625
  TL <- a2 * hgt / uf
  # Calculate humidity limits
  ez <- satvap(tair) * (rh / 100)
  eh <- ez[n]
  th <- tair[n]
  eg <- satvap(tground) * (soilrh / 100)
  emn <- min(c(eg, eh, ez))
  emx <- max(c(eg, eh, ez))
  # Calculate position variance
  ow <- uf * (0.75 + 0.5 * cos(pi * (1 - z / hgt)))
  # Calculate thermal diffusivity
  KH <- TL * ow^2
  SS <- paii * Lz
  # Compute near-field correction factor for small sample size
  n <- length(paii)
  mu <- 1 + 0.894 * exp(-0.01386 * n) + 9.82 * exp(-0.15 * n)
  # Compute near-field concentrations at the top of the canopy
  dz1 <- (hgt - z[-n]) / (ow[-n] * TL)
  kn  <- -0.39894 * log1p(-exp(-dz1)) - 0.15623 * exp(-dz1)
  Cnh <- sum((SS[-n] / ow[-n]) * kn * (2 * hgt / (ow[-n] * TL))) * mu
  # Compute far field concentration at top of canopy
  lah <- ifelse(th < 0, 51078.69 - 4.338 * th - 0.06367 * th^2,
                45068.7 - 42.8428 * th)
  Cfh <- eh * rhohair(th, pk) * lah / pk
  # Compute near and far-field concentrations for each canopy element
  sumRH <- 0
  for (i in 1:n) {
    # Compute resistance from ground to z
    RH <- 1 / KH[i]
    sumRH <- sumRH + RH
    rHa = sumRH * dz;
    if (rHa < 2) rHa = 2.0;
    ph <- rhohair(tair[i], pk)
    la <- ifelse(th < 0, 51078.69 - 4.338 * tleaf[i] - 0.06367 * tleaf[i]^2,
                 45068.7 - 42.8428 * tleaf[i])
    GL <- (la / (rHa * pk)) * (eg - ez[i]) * dz
    L <- sum(SS[1:i]) + GL
    # Compute far-field concentration
    Cf <- sum((L / KH[i:n]) * dz)
    # Compute near-field concentration
    dz1 <- (z[i] - z) / (ow * TL)
    dz2 <- (z[i] + z) / (ow * TL)
    zeta <- abs(dz1)
    e <- exp(-zeta)
    kn <- -0.39894 * log1p(-e) - 0.15623 * e
    Cn <- sum((SS[-i] / ow[-i]) * kn[-i] * (dz1[-i] + dz2[-i]))
    # Compute total source Concetration
    CL <- Cfh - Cnh + Cf + Cn * mu
    ean <- (CL * pk) / (la * ph)
    if (ean > emx) ean <- emx
    if (ean < emn) ean <- emn
    rh[i] = (ean / satvap(tair[i])) * 100.0
    if (rh[i] > 100) rh[i] <- 100
    if (rh[i] <  20) rh[i] <-  20
  }
  return(rh)
}
