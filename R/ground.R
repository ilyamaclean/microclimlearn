#' Internal function for calculating geometric spacings
.geometric <- function(n, totalDepth) {
  i <- seq_len(n)
  weightSum <- sum(i^2)
  dz_unit <- totalDepth / weightSum
  z <- numeric(n + 2)
  z[1] <- 0
  z[2] <- dz_unit
  if (n > 1) z[3:(n + 1)] <- dz_unit + cumsum(dz_unit * i[-1]^2)
  z
}
#' Internal function for used my soil water model
.dTheta_dPsi <- function(psiw, psie, b, Smax) {
  psie <- -abs(psie)
  Se <- (psiw / psie)^(-1 / b)
  theta <- Smax * Se
  out <- -theta / (b * psiw)
  out[psiw >= psie] <- 0
  return(out)
}
#' Internal function for used my soil water model
.dvapor_dPsi <- function(psiw, psie, theta, b, Smax, Tk) {
  humidity <- exp(0.018015 * psiw / (8.314 * Tk))
  dTdP <- .dTheta_dPsi(psiw, psie, b, Smax)
  capacity_vapor <- (Smax - theta) * 0.017 * humidity *
    (0.018015 / (8.314 * Tk)) - dTdP * 0.017 * humidity
  return(capacity_vapor)
}
#' Internal function used to distribute roots zone
.root_distribute <- function(dz, totalDepth, skew) {
  n <- length(dz)
  if (skew == 0) skew <- 1e-12
  i <- seq_len(n) - 1
  z <- (i / n) * totalDepth
  v <- exp(-z * skew) * dz
  # fix last element (matches C++ behaviour)
  if (n > 1) v[n] <- exp(-z[n] * skew) * dz[n - 1]
  v / sum(v)
}
#' Internal function used to distribute transpiration in root zone
.transpiration_distribute <- function(soilp, rootfrac, totalTransp_mm, dT, psiw, p = 0.5) {
  .alpha_wet <- function(psiw, psie) {
    ifelse(psiw >= 0, 0, ifelse(psiw <= psie, 1, abs(psiw) / abs(psie)))
  }
  .alpha_dry <- function(psiw, psi_dry, psi_wilt) {
    ifelse(psiw <= psi_wilt, 0, ifelse(psiw >= psi_dry, 1,
                                       (psiw - psi_wilt) / (psi_dry - psi_wilt)))
  }
  n <- length(psiw)
  S <- numeric(n)
  if (totalTransp_mm > 0) {
    Trate <- totalTransp_mm / dT
    theta_dry <- soilp$Smin[seq_len(n)] +
      p * (soilp$Smax[seq_len(n)] - soilp$Smin[seq_len(n)])
    psi_wilt <- PsiFromtheta(
      theta = soilp$Smin[seq_len(n)],
      psie = soilp$psie[seq_len(n)],
      b = soilp$b[seq_len(n)],
      Smax = soilp$Smax[seq_len(n)],
      MPa = FALSE
    )
    psi_dry <- PsiFromtheta(
      theta = theta_dry,
      psie = soilp$psie[seq_len(n)],
      b = soilp$b[seq_len(n)],
      Smax = soilp$Smax[seq_len(n)],
      MPa = FALSE
    )
    aw <- .alpha_wet(psiw, soilp$psie[seq_len(n)])
    ad <- .alpha_dry(psiw, psi_dry, psi_wilt)
    alpha <- aw * ad
    w <- rootfrac * alpha
    sumw <- sum(w)
    if (sumw > 0) {
      w <- w / sumw
    } else {
      w[] <- 0
    }
    S <- w * Trate
  }
  return(S)
}
#' Create soil parameter list for running the soil heat and water models
#'
#' Creates a list of soil physical parameters for a selected soil type and
#' replicates them across soil layers, with optional adjustment for slope,
#' aspect, surface organic matter, and deep boundary saturation.
#'
#' @param soiltype Character string giving the soil type. See details. Default is \code{"Clay loam"}.
#' @param nlayers Integer. Number of soil layers.
#' @param totalDepth Numeric. Total soil depth (m). Internally this is multiplied
#'   by 1.5 to allow for a boundary layer.
#' @param slope Numeric. Surface slope (degrees).
#' @param aspect Numeric. Surface aspect (degrees).
#' @param surface_organicmu Numeric. Multiplication factor controlling the
#'   vertical profile of soil organic matter near the surface.
#' @param FreeDrain Logical. If `TRUE`, a free-drainage lower boundary is
#' assumed, allowing water to drain out of the bottom of the soil profile under
#' gravity. If `FALSE`, the lower boundary is treated as saturated, meaning the
#' soil is assumed to be sitting above a water table, which limits drainage from
#' the deepest layer.
#' @details
#' The selected soil parameter values are replicated across \code{nlayers + 1}
#' soil nodes. Organic matter content is adjusted vertically with corresponding
#' rescaling of mineral and related soil constituents to preserve total volume.
#' Bulk density is also adjusted slightly with depth.
#'
#' Available soil types (argument \code{soiltype}) are:
#' \itemize{
#'   \item Sand
#'   \item Loamy sand
#'   \item Sandy loam
#'   \item Loam
#'   \item Silt loam
#'   \item Sandy clay loam
#'   \item Clay loam
#'   \item Silty clay loam
#'   \item Sandy clay
#'   \item Silty clay
#'   \item Clay
#'   \item Silt
#' }
#' @return A list of soil parameters:
#' \describe{
#'   \item{Smax}{Volumetric soil water content at saturation (m^3 m^-3), repeated for each soil layer}
#'   \item{Smin}{Residual volumetric soil water content (m^3 m^-3), repeated for each soil layer}
#'   \item{n}{Pore size distribution parameter (unitless), repeated for each soil layer}
#'   \item{Ksat}{Saturated hydraulic conductivity (kg s / m^3), repeated for each soil layer}
#'   \item{Vq}{Quartz fraction by volume (m^3 m^-3), repeated for each soil layer}
#'   \item{Vm}{Mineral fraction by volume (m^3 m^-3), repeated for each soil layer}
#'   \item{Vo}{Organic fraction by volume (m^3 m^-3), repeated for each soil layer and adjusted with depth}
#'   \item{Mc}{Coarse fragment or gravel fraction (m^3 m^-3), repeated for each soil layer}
#'   \item{rho}{Bulk density (kg m^-3), repeated for each soil layer and adjusted with depth}
#'   \item{b}{Campbell soil water retention parameter (unitless), repeated for each soil layer}
#'   \item{psi_e}{Air-entry water potential (m), repeated for each soil layer}
#'   \item{gref}{Ground shortwave reflectance (unitless)}
#'   \item{groundem}{Ground longwave emissivity (unitless)}
#'   \item{grefPAR}{Ground reflectance for photosynthetically active radiation (unitless)}
#'   \item{slope}{Surface slope (degrees)}
#'   \item{aspect}{Surface aspect (degrees)}
#'   \item{nLayers}{Number of soil layers}
#'   \item{totalDepth}{Total soil depth used internally (m), after multiplying input depth by 1.5}
#'   \item{FreeDrain}{Logical indicating whether the deep boundary condition is saturated}
#' }
#' @export
createsoilplist <- function(soiltype = "Clay loam", nlayers = 15, totalDepth = 2, slope = 0, aspect = 180,
                            surface_organicmu = 3, FreeDrain = TRUE) {
  soilparams <- micropoint::newsoilparamstable
  totalDepth <- 1.5 * totalDepth
  s <- which(soilparams$Soil.type == soiltype)
  v <- as.numeric(soilparams[s, 2:13])
  nms <- names(soilparams)[2:13]
  soilc <- as.list(v)
  for (i in 1:12) {
    soilc[[i]] <- rep(v[i], nlayers + 1)
  }
  names(soilc) <- nms
  mu <- rev(.geometric(nlayers, surface_organicmu))[2:(nlayers + 2)]
  soilc$Vo <- soilc$Vo * mu
  dif <- soilc$Vo - v[8]
  sm <- v[6] + v[7]
  mu <- (sm - dif) / sm
  soilc$Vm <- soilc$Vm * mu
  soilc$Vq <- soilc$Vq * mu
  soilc$Mc <- soilc$Mc * mu
  mu <- seq(0.9, 1.1, length.out = nlayers + 1)
  soilc$rho <- soilc$rho * mu
  soilc$gref <- 0.15
  soilc$groundem <- 0.97
  soilc$grefPAR <- 0.75 * soilc$gref
  soilc$slope <- slope
  soilc$aspect <- aspect
  soilc$nLayers <- nlayers
  soilc$totalDepth <- totalDepth
  soilc$FreeDrain <- FreeDrain
  soilc$alpha <- NULL
  soilc$psi_min <- PsiFromtheta(soilc$Smin, soilc$psi_e, soilc$b, soilc$Smax, MPa = FALSE)
  soilc$psie <- -abs(soilc$psi_e)
  soilc$Psie <- NULL
  return(soilc)
}
#' @title Initalizes input for the soil heat model
#' @description Creates the soil state object used by the soil heat model
#' @param soilp list of soil parameters as returned by [createsoilplist()]
#' @param Te vector of soil layer temperatures (degrees C)
#' @param wc vector of volumetric soil water contents (m3 m-3)
#' @return a list comprising
#' \describe{
#'   \item{n}{number of soil layers}
#'   \item{z}{depth of soil nodes (m)}
#'   \item{dz}{soil layer thicknesses (m)}
#'   \item{zCenter}{depth of layer centres (m)}
#'   \item{vol}{soil layer volumes per unit ground area (m)}
#'   \item{wc}{volumetric soil water content of each layer (m3 m-3)}
#'   \item{Te}{soil layer temperature (degrees C)}
#'   \item{oldTe}{soil layer temperature at previous time step (degrees C)}
#'   \item{Gflux}{soil heat flux}
#'   \item{iters}{number of solver iterations}
#' }
#' @export
InitSoilheatmod <- function(soilp, Te, wc) {
  n <- soilp$nLayers
  totalDepth <- soilp$totalDepth
  z <- .geometric(n, totalDepth)
  dz <- z[2:(n + 2)] - z[1:(n + 1)]
  zCenter <- z[1:(n + 1)] + dz * 0.5
  state <- list(n = n, z = z, dz = dz, zCenter = zCenter, vol = dz,
                wc = wc, oldTe = Te, Te = Te, Gflux = 0, iters = 0)
}

#' @title Initializes input for the soil water model
#' @description Creates the soil state object used by the soil water model from
#' an existing soil heat model state
#' @param soilp list of soil parameters as returned by [createsoilplist()]
#' @param soilheatmod soil heat model state as returned by [InitSoilheatmod()]
#' @param rootskew parameter controlling vertical distribution of roots
#' @return a list comprising
#' \describe{
#'   \item{n}{number of soil layers}
#'   \item{z}{depth of soil nodes (m)}
#'   \item{dz}{soil layer thicknesses (m)}
#'   \item{zCenter}{depth of layer centres (m)}
#'   \item{Tc}{soil temperature profile (degrees C)}
#'   \item{theta}{volumetric soil water content (m3 m-3)}
#'   \item{vol}{soil layer volumes per unit ground area (m)}
#'   \item{psiw}{soil water potential (J kg-1)}
#'   \item{k}{hydraulic conductivity (m s-1)}
#'   \item{vapor}{soil vapour content}
#'   \item{oldvapor}{soil vapour content at previous time step}
#'   \item{oldtheta}{soil water content at previous time step (m3 m-3)}
#'   \item{rootfrac}{fraction of roots in each soil layer}
#'   \item{success}{logical indicating convergence success - set to TRUE on initalization}
#'   \item{iterations}{number of solver iterations - set to zero initially}
#' }
#' @export
InitSoilwatermod <- function(soilp, soilheatmod, rootskew) {
  out <- list()
  out$n = soilheatmod$n # number of layers
  out$z = soilheatmod$z # node depths
  out$dz = soilheatmod$dz # layer thickness
  out$zCenter = soilheatmod$zCenter # centre between nodes
  out$Tc = soilheatmod$Te # temperature of soil profile
  out$oldTc = soilheatmod$oldTe # old temperature of soil profile
  out$theta = soilheatmod$wc # volumetric water content of soil profile
  idx1 <- seq_len(out$n)
  idx2 <- seq_len(out$n) + 2
  out$vol <- c(0, (out$z[idx2] - out$z[idx1]) / 2)
  out$psiw <- PsiFromtheta(out$theta, soilp$psie, soilp$b, soilp$Smax, MPa = FALSE)
  out$k <- hydraulicConductivityFromTheta(out$theta, soilp$b, soilp$Smax, soilp$Ksat);
  out$vapor <- vaporFromPsi(out$psiw, out$theta, out$Tc + 273.15, soilp$Smax)
  out$oldvapor <- out$vapor
  out$oldtheta <- soilheatmod$wc
  totalDepth <- max(out$z)
  out$rootfrac <- .root_distribute(out$dz, totalDepth, rootskew)
  out$success <- TRUE
  out$iterations <- 0
  out$iterations <- 0
  return(out)
}
#' @title Calculates soil surface relative humidity
#' @param theta volumetric soil water content (m^3 / m^3)
#' @param tsoil soil temperature (degrees C)
#' @param psie air-entry water potential
#' @param b Campbell soil water retention parameter
#' @param Smax saturated volumetric water content (m^3 / m^3)
#' @return relative humidity (Percentage)
#' @export
soilrelhum <- function(theta, tsoil, psie, b, Smax) {
  psie <- -abs(psie)
  psiw <- psie * (theta / Smax)^-b
  Tk <- tsoil + 273.15
  hr <- exp(0.018015 * psiw / (8.314 * Tk)) * 100
  hr
}
#' @title Calculates soil thermal conductivity
#' @description Calculates soil thermal conductivity (W m-1 K-1)
#' @param Vq volumetric quartz fraction of soil (m^3 / m^3)
#' @param Vm volumetric mineral fraction of soil (m^3 / m^3)
#' @param Vo volumetric organic fraction of soil (m^3 / m^3)
#' @param Vw volumetric water content of soil (m^3 / m^3)
#' @param Mc mass fraction of clay (kg / kg)
#' @param Tc soil temperature (degrees C)
#' @param pk air pressure (kPa)
#' @return soil thermal conductivity (W m-1 K-1)
#' @export
thermalConductivity <- function(Vq, Vm, Vo, Vw, Mc, Tc, pk) {
  # solids
  Vsolid <- Vq + Vm + Vo
  kSolid <- 2.5^(Vq / Vsolid) * 8.8^(Vm / Vsolid) * 0.25^(Vo / Vsolid)
  # shape factor
  Ga <- 0.088 * ((Vq + Vm) / Vsolid) + 0.5 * (Vo / Vsolid)
  Q <- 7.25 * Mc + 2.52
  xwo <- 0.33 * Mc + 0.078
  porosity <- 1 - Vsolid
  gasPorosity <- pmax(porosity - Vw, 0)
  Tk <- Tc + 273.15
  Lv <- 45144 - 48 * Tc
  svp <- 0.611 * exp(17.502 * Tc / (Tc + 240.97))
  slope <- 17.502 * 240.97 * svp / (240.97 + Tc)^2
  Dv <- 0.0000212 * (101.3 / pk) * (Tk / 273.15)^1.75
  rhoAir <- 44.65 * (pk / 101.3) * (273.15 / Tk)
  stcor <- 1 - svp / pk
  stcor <- pmax(stcor, 0.3)
  kWater <- 0.56 + 0.0018 * Tc
  wf <- ifelse(Vw < 0.01 * xwo, 0, 1 / (1 + (Vw / xwo)^(-Q)))
  kGas <- 0.0242 + 0.00007 * Tc + wf * Lv * rhoAir * Dv * slope / (pk * stcor)
  Gc <- 1 - 2 * Ga
  kFluid <- kGas + (kWater - kGas) * (Vw / porosity)^2
  ka <- (2 / (1 + (kGas / kFluid - 1) * Ga) +
           1 / (1 + (kGas / kFluid - 1) * Gc)) / 3
  kw <- (2 / (1 + (kWater / kFluid - 1) * Ga) +
           1 / (1 + (kWater / kFluid - 1) * Gc)) / 3
  ks <- (2 / (1 + (kSolid / kFluid - 1) * Ga) +
           1 / (1 + (kSolid / kFluid - 1) * Gc)) / 3

  kout <- (kw * kWater * Vw + ka * kGas * gasPorosity + ks * kSolid * Vsolid) /
    (kw * Vw + ka * gasPorosity + ks * Vsolid)
  return(kout)
}
#' @title Calculates soil volumetric heat capacity
#' @description Calculates soil volumetric heat capacity (J m-3 K-1)
#' @param Vq volumetric quartz fraction of soil (m^3 / m^3)
#' @param Vm volumetric mineral fraction of soil (m^3 / m^3)
#' @param Vo volumetric organic fraction of soil (m^3 / m^3)
#' @param Vw volumetric water content of soil (m^3 / m^3)
#' @param Tc soil temperature (degrees C)
#' @param pk air pressure (kPa)
#' @return soil volumetric heat capacity (J m-3 K-1)
#' @export
heatCapacity <- function(Vq, Vm, Vo, Vw, Tc, pk) {
  Vsum <- Vq + Vm + Vo + Vw
  Va <- ifelse(Vsum > 1, 0, 1 - Vsum)
  Vq <- ifelse(Vsum > 1, Vq / Vsum, Vq)
  Vm <- ifelse(Vsum > 1, Vm / Vsum, Vm)
  Vo <- ifelse(Vsum > 1, Vo / Vsum, Vo)
  Vw <- ifelse(Vsum > 1, Vw / Vsum, Vw)
  CHa <- cpair(Tc) * rhohair(Tc, pk) / 1e6
  if (length(Tc) == 1L) Tc <- rep(Tc, length(Vw))
  ChW <- 1.93 + 0.0067 * Tc
  ChW[Tc >= 0] <- 4.18
  CH <- Vq * 2.13 + Vm * 2.31 + Vo * 2.50 + Va * CHa + Vw * ChW
  CH * 1e6
}
#' @title Calculates soil surface energy balance
#' @description Calculates the soil surface energy balance (W m-2)
#' @param soilp list of soil parameters as returned by [createsoilplist()]
#' @param Rabs absorbed radiation at the surface (W m-2)
#' @param Tref reference air temperature (degrees C)
#' @param Tsurface soil surface temperature (degrees C)
#' @param pk air pressure (kPa)
#' @param relhum air relative humidity (Percentage)
#' @param rHa aerodynamic resistance to heat from ground to reference height
#' as returned by [groundresistance()] (s m-1)
#' @param theta volumetric soil water content (m3 m-3)
#' @return soil surface energy balance (W m-2)
#' @export
soilsurfaceEB <- function(soilp, Rabs, Tref, Tsurface, pk, relhum, rHa, theta) {
  # net radiation
  sb <- 5.67e-8
  Rnet <- Rabs - soilp$groundem * sb * (Tsurface + 273.15)^4
  # sensible heat
  cp <- cpair(Tref)
  ph <- rhohair(Tref, pk)
  H <- (ph * cp / rHa) * (Tsurface - Tref)
  # latent heat
  hr <- soilrelhum(theta, Tsurface, soilp$psie[1], soilp$b[1], soilp$Smax[1])
  es <- satvap(Tsurface) * (hr / 100)
  ea <- satvap(Tref) * (relhum / 100)
  la <- latvap(Tsurface)
  L <- (la * ph / (rHa * pk)) * (es - ea)
  Rnet - H - L
}
#' @title Calculates soil evaporation flux
#' @param theta volumetric soil surface water content (m3 m-3)
#' @param Tsurface soil surface temperature (degrees C)
#' @param Tair air temperature (degrees C)
#' @param relhum air relative humidity (Percentage)
#' @param rHa aerodynamic resistance to heat as returned by [groundresistance()] (s m-1)
#' @param psie air-entry water potential (J kg-1)
#' @param b Campbell soil water retention parameter
#' @param Smax volumetric soil water content at saturation (m3 m-3)
#' @param dT time step (s)
#' @return soil evaporation (mm over timestep)
#' @export
evaporation_flux <- function(theta, Tsurface, Tair, relhum, rHa, psie, b, Smax, dT = 3600) {
  psiw <- PsiFromtheta(theta, Smax, psie, b, MPa = FALSE)
  Tk <- Tair + 273.15 # Temperature deg C -> K
  hs <- exp(0.018015 * psiw / (8.314 * Tk)) #  dimensionless soil effective relative humidity (0 to 1 range)
  es <- 1000 * satvap(Tsurface) * hs # kPa -> Pa
  ea <- 1000 * satvap(Tair) * (relhum / 100.0)  # kPa -> Pa
  Ev <- (0.018015 / (rHa * 8.314 * Tk)) * (es - ea) * dT # Bare soil evaporation (kg m^-2 s^-1 -> mm)
  return (Ev)
}
#' @title Solves a tridiagonal system of equations
#' @description Solves a tridiagonal system of linear equations using the Thomas algorithm
#' @param aa vector of sub-diagonal coefficients
#' @param bb vector of diagonal coefficients
#' @param cc vector of super-diagonal coefficients
#' @param dd vector of right-hand-side values
#' @param x vector of solution values
#' @param first index of first equation to solve
#' @param last index of last equation to solve
#' @return a list comprising
#' \describe{
#'   \item{bb}{updated diagonal coefficients}
#'   \item{cc}{updated super-diagonal coefficients}
#'   \item{dd}{updated right-hand-side values}
#'   \item{x}{solution vector}
#' }
#' @export
ThomasBoundaryCondition <- function(aa, bb, cc, dd, x, first = 1, last = length(x)) {
  for (i in first:(last - 1)) {
    cc[i] <- cc[i] / bb[i]
    dd[i] <- dd[i] / bb[i]
    bb[i + 1] <- bb[i + 1] - aa[i + 1] * cc[i]
    dd[i + 1] <- dd[i + 1] - aa[i + 1] * dd[i]
  }
  x[last] <- dd[last] / bb[last]
  if (last - 1 >= first) {
    for (i in seq(last - 1, first, by = -1)) {
      x[i] <- dd[i] - cc[i] * x[i + 1]
    }
  }
  list(bb = bb, cc = cc, dd = dd, x = x)
}
#' @title Updates soil temperatures using the soil heat model
#' @description Solves soil heat conduction through the soil profile over one time step,
#' updating layer temperatures by combining surface energy balance from [soilsurfaceEB()],
#' thermal conductivity from [thermalConductivity()], volumetric heat capacity from
#' [heatCapacity()], and tridiagonal solution with [ThomasBoundaryCondition()].
#' @param soilheatmod list of variables provided as input and updated by soil
#' heat model as returned by [InitSoilheatmod()] on first step or by previous model run
#' @param soilp list of soil parameters as returned by [createsoilplist()]
#' @param Rabs absorbed radiation at the ground surface (W m-2)
#' @param Tref reference air temperature (degrees C)
#' @param relhum air relative humidity (Percentage)
#' @param atmPressure air pressure (kPa)
#' @param rHa aerodynamic resistance to heat from ground to reference height
#' as returned by [groundresistance()] (s m-1)
#' @param dT time step (s)
#' @param Fact weighting factor for the implicit solution scheme
#' @param maxNrIterations maximum number of nonlinear iterations
#' @param tolerance convergence criterion for maximum change in layer temperature
#' between successive iterations (degrees C)
#' @return a list comprising
#' \describe{
#'   \item{n}{number of soil layers}
#'   \item{z}{depth of soil nodes (m)}
#'   \item{dz}{soil layer thicknesses (m)}
#'   \item{zCenter}{depth of layer centres (m)}
#'   \item{vol}{soil layer volumes per unit ground area (m)}
#'   \item{wc}{volumetric soil water content of each layer (m3 m-3)}
#'   \item{Te}{updated soil layer temperature (degrees C)}
#'   \item{oldTe}{soil layer temperature at previous time step (degrees C)}
#'   \item{Gflux}{integrated soil heat storage flux over the profile (W m-2)}
#'   \item{iters}{number of iterations used to reach convergence}
#' }
#' @export
SoilHeatModel <- function(soilheatmod, soilp, Rabs, Tref, relhum, atmPressure, rHa, dT = 3600, Fact = 0.5,
                          maxNrIterations = 100, tolerance = 1e-2) {

  n <- soilheatmod$n
  gg <- 1 - Fact
  boundaryT <- soilheatmod$oldTe[n + 1]
  oldTe_fixed <- soilheatmod$oldTe
  Te_new <- soilheatmod$Te
  wc <- soilheatmod$wc
  dz <- soilheatmod$dz
  Vq <- soilp$Vq
  Vm <- soilp$Vm
  Vo <- soilp$Vo
  Mc <- soilp$Mc
  max_qsurface <- 0
  nrIterations <- 0
  maxdT <- Inf
  qsurface <- 0
  qsurf1 <- NA_real_
  qsurf2 <- NA_real_
  while (maxdT > tolerance && nrIterations < maxNrIterations) {
    if (nrIterations < 10) {
      Tav <- 0.5 * (oldTe_fixed[1] + Te_new[1])
      qsurface <- soilsurfaceEB(soilp, Rabs, Tref, Tav, atmPressure, relhum, rHa, wc[1])
      if (nrIterations < 2 && abs(qsurface) > abs(max_qsurface)) {
        max_qsurface <- abs(qsurface)
      }
      if (nrIterations > 1 && abs(qsurface) > abs(max_qsurface)) {
        qsurface <- sign(qsurface) * abs(max_qsurface)
      }
      if (nrIterations == 8) qsurf1 <- qsurface
      if (nrIterations == 9) qsurf2 <- qsurface
    } else if (nrIterations == 10) {
      qsurface <- 0.5 * (qsurf1 + qsurf2)
    }
    lambda <- thermalConductivity(soilp$Vq, soilp$Vm, soilp$Vo, wc, Mc, Te_new, atmPressure)
    CT <- heatCapacity(Vq, Vm, Vo, wc, Te_new, atmPressure) * soilheatmod$vol
    ff <- numeric(n + 1)
    ff[1] <- if (lambda[1] == lambda[2]) {
      lambda[1] / dz[1]
    } else {
      ((lambda[1] - lambda[2]) / log(lambda[1] / lambda[2])) / dz[1]
    }
    if (n > 1) {
      ff[2:n] <- lambda[2:n] / dz[2:n]
    }
    aa <- numeric(n + 1)
    bb <- numeric(n + 1)
    cc <- numeric(n + 1)
    dd <- numeric(n + 1)
    if (n > 2) {
      idx <- 2:(n - 1)
      aa[idx] <- -ff[idx - 1] * Fact
      bb[idx] <- CT[idx] / dT + (ff[idx - 1] + ff[idx]) * Fact
      cc[idx] <- -ff[idx] * Fact
      dd[idx] <- CT[idx] / dT * oldTe_fixed[idx] + gg * (ff[idx - 1] *
                                                           oldTe_fixed[idx - 1] + ff[idx] * oldTe_fixed[idx + 1] -
                                                           (ff[idx - 1] + ff[idx]) * oldTe_fixed[idx])
    }
    aa[1] <- 0
    bb[1] <- CT[1] / dT + ff[1]
    cc[1] <- -ff[1]
    dd[1] <- CT[1] / dT * oldTe_fixed[1] + qsurface
    aa[n] <- 0
    bb[n] <- 1
    cc[n] <- 0
    dd[n] <- boundaryT
    Te_prev <- Te_new
    TBC <- ThomasBoundaryCondition(aa, bb, cc, dd, Te_new, 1, n - 1)
    Te_new <- TBC$x
    if (n > 2) {
      idx <- 2:(n - 1)
      lo <- pmin(Te_new[idx - 1], Te_new[idx + 1])
      hi <- pmax(Te_new[idx - 1], Te_new[idx + 1])
      Te_new[idx] <- pmax(pmin(Te_new[idx], hi), lo)
    }
    maxdT <- max(abs(Te_new - Te_prev))
    nrIterations <- nrIterations + 1
  }
  G <- if (n > 1) {
    sum(CT[2:n] * (Te_new[2:n] - oldTe_fixed[2:n]) / dT)
  } else {
    0
  }
  soilheatmod$Te <- Te_new
  soilheatmod$Gflux <- G
  soilheatmod$iters <- nrIterations
  return(soilheatmod)
}
#' Create list of Campbell soil hydraulic parameters
#' @param soiltype Character string giving the soil type. See details
#' @return A list of soil parameters:
#' \describe{
#'   \item{psie}{Air-entry water potential (m)}
#'   \item{n}{Pore size distribution parameter (unitless)}
#'   \item{b}{Campbell soil water retention parameter (unitless)}
#'   \item{Smax}{Volumetric soil water content at saturation (m^3 m^-3)}
#'   \item{Smin}{Residual volumetric soil water content (m^3 m^-3)}
#'   \item{Ksat}{Saturated hydraulic conductivity (kg s / m^3)}
#' }
#' @details:
#' Available soil types (argument \code{soiltype}) are:
#' \itemize{
#'   \item Sand
#'   \item Loamy sand
#'   \item Sandy loam
#'   \item Loam
#'   \item Silt loam
#'   \item Sandy clay loam
#'   \item Clay loam
#'   \item Silty clay loam
#'   \item Sandy clay
#'   \item Silty clay
#'   \item Clay
#'   \item Silt
#' }
#' @export
Campbellparams <- function(soiltype) {
  soilparams <- micropoint::newsoilparamstable
  s <- which(soilparams$Soil.type == soiltype)
  out <- list()
  out$psie <- soilparams$psi_e[s]
  out$n <- soilparams$n[s]
  out$b <- soilparams$b[s]
  out$Smin <- soilparams$Smin[s]
  out$Smax <- soilparams$Smax[s]
  out$Ksat <- soilparams$Ksat[s]
  return(out)
}
#' Create list of Van Genuchten soil hydraulic parameters
#' @param soiltype Character string giving the soil type. See details
#' @return A list of soil parameters:
#' \describe{
#'   \item{psie}{Air-entry water potential (m)}
#'   \item{alpha}{Van Genuchten shapte parameter (unitless)}
#'   \item{n}{Pore size distribution parameter (unitless)}
#'   \item{Smax}{Volumetric soil water content at saturation (m^3 m^-3)}
#'   \item{Smin}{Residual volumetric soil water content (m^3 m^-3)}
#'   \item{Ksat}{Saturated hydraulic conductivity (kg s / m^3)}
#' }
#' @details:
#' Available soil types (argument \code{soiltype}) are:
#' \itemize{
#'   \item Sand
#'   \item Loamy sand
#'   \item Sandy loam
#'   \item Loam
#'   \item Silt loam
#'   \item Sandy clay loam
#'   \item Clay loam
#'   \item Silty clay loam
#'   \item Sandy clay
#'   \item Silty clay
#'   \item Clay
#'   \item Silt
#' }
#' @export
VGparams <- function(soiltype) {
  soilparams <- micropoint::newsoilparamstable
  s <- which(soilparams$Soil.type == soiltype)
  out <- list()
  out$psie <- soilparams$VGpsie[s]
  out$alpha <- soilparams$alpha[s]
  out$n <- soilparams$VGn[s]
  out$Smin <- soilparams$Smin[s]
  out$Smax <- soilparams$Smax[s]
  out$Ksat <- soilparams$Ksat[s]
  return(out)
}
#' @title Derives soil water potential from volumetric water content using the Campbell method
#' @param theta volumetric soil water fraction (m^3 / m^3)
#' @param Smax volumetric soil water fraction at saturation (m^3 / m^3)
#' @param psi_e air entry water potential (J/m^3)
#' @param b Campbell soil water retention parameter
#' @param MPa optional logical indicating whether to return outputs in Mpa (TRUE) or kPa (FALSE)
#' @returns soil water potential (MPa or kPa ~ J/Kg)
#' @export
PsiFromtheta <- function(theta, psie, b, Smax, MPa = TRUE) {
  psie <- -abs(psie)
  Se <- theta / Smax
  Se[Se > 1] <- 1
  psiw <- psie * Se^(-b) # J/kg
  if (MPa) psiw <- psiw * 0.001 # Convert to MPa
  return(psiw)
}
#' @title Derives soil water potential from volumetric water content using the Van Genuchten method
#' @param theta volumetric soil water fraction (m^3 / m^3)
#' @param Smin Residual volumetric soil water fraction (m^3 / m^3)
#' @param Smax volumetric soil water fraction at saturation (m^3 / m^3)
#' @param alpha Van Genuchten shape parameter
#' @param n Van Genuchten pore size distribution parameter
#' @returns soil water potential (kPa ~ J/Kg)
#' @export
PsiFromthetaVG <- function(theta, Smin, Smax, alpha, n) {
  # Calculate water saturation at the air-entry potential
  m <- 1 - (1 / n)
  # Calculate degree of saturation
  Se <- (theta - Smin) / (Smax - Smin)
  Se <- pmin(1 - 1e-12, pmax(1e-12, Se))
  psiw <- -(1 / alpha) * (Se^(-1 / m) - 1)^(1 / n)
  return(psiw)
}
#' @title Calculates soil vapour content from water potential
#' @description Calculates soil vapour content from soil water potential using the Kelvin equation
#' @param psiw soil water potential (J kg-1)
#' @param theta volumetric soil water content (m3 m-3)
#' @param Tk soil temperature (K)
#' @param Smax volumetric soil water content at saturation (m3 m-3)
#' @return soil vapour content
#' @export
vaporFromPsi <- function(psiw, theta, Tk, Smax) {
  humidity <- exp(0.018015 * psiw / (8.314 * Tk))
  vapor <- (Smax - theta) * 0.017 * humidity
  return(vapor)
}
#' @title Calculates vapour conductivity from soil water potential
#' @description Calculates soil vapour conductivity as a function of water potential and moisture content
#' @param psiw soil water potential (J kg-1)
#' @param theta volumetric soil water content (m3 m-3)
#' @param Tk soil temperature (K)
#' @param Smax volumetric soil water content at saturation (m3 m-3)
#' @return soil vapour conductivity (kg^2 m-1 s-1 J-1)
#' @export
vaporConductivityFromPsiTheta <- function(psiw, theta, Tk, Smax) {
  humidity <- exp(0.018015 * psiw / (8.314 * Tk))
  k <- 0.66 * (Smax - theta) * 0.000024 * 0.017 * humidity * 0.018015 / (8.314 * Tk)
  return(k)
}
#' @title Calculates volumetric soil water fraction from water potential using the Campbell method
#' @param psiw soil water potential (J kg-1)
#' @param psie air-entry water potential (J kg-1)
#' @param b Campbell soil water retention parameter
#' @param Smax volumetric soil water content at saturation (m3 m-3)
#' @return volumetric soil water content (m3 m-3)
#' @export
thetaFromPsi <- function(psiw, psie, b, Smax) {
  psie <- -abs(psie)
  Se <- pmin((psiw / psie) ^ (-1 / b), 1)
  theta <- Se * Smax
  return(theta)
}
#' @title Calculates volumetric soil water fraction from water potential using the Van Genuchten method
#' @param psiw soil water potential (kPa ~ J/Kg)
#' @param Smin Residual volumetric soil water fraction (m^3 / m^3)
#' @param Smax volumetric soil water fraction at saturation (m^3 / m^3)
#' @param psie air entry water potential (J/m^3)
#' @param alpha Van Genuchten shape parameter
#' @param n Van Genuchten pore size distribution parameter
#' @returns volumetric soil water content (m3 m-3)
#' @export
thetaFromPsiVG <- function(psiw, Smin, Smax, psie, alpha, n) {
  psie <- -abs(psie)
  # Calculate water saturation at the air-entry potential
  m <- 1 - (1 / n)
  Sc <- (1 + (alpha * abs(psie))^n)^(-m)
  # Calculate degree of saturation
  Se <- (1 / Sc) * (1 + (alpha * abs(psiw))^n)^(-m)
  Se <- pmin(1, pmax(0, Se))
  # Calculate theta
  theta <- Se * (Smax - Smin) + Smin
  return(theta)
}
#' @title Calculates hydraulic conductivity from soil water content using Campbell method
#' @param theta volumetric soil water content (m3 m-3)
#' @param psie air-entry water potential (J kg-1)
#' @param b Campbell soil water retention parameter
#' @param Smax volumetric soil water content at saturation (m3 m-3)
#' @param Ksat saturated hydraulic conductivity (kg s / m^3)
#' @param n Campbell pore size distribution parameter
#' @return hydraulic conductivity (kg s / m^3)
#' @export
hydraulicConductivityFromTheta <- function(theta, b, Smax, Ksat) {
  n <- 2*b+3
  k <- Ksat * (theta / Smax)^n
  return(k)
}
#' @title Calculates hydraulic conductivity from soil water content using Van Genuchten method
#' @param theta volumetric soil water content (m3 m-3)
#' @param Smin Residual volumetric soil water fraction (m^3 / m^3)
#' @param Smax volumetric soil water fraction at saturation (m^3 / m^3)
#' @param psie air-entry water potential (J kg-1)
#' @param alpha van Genuchten shape parameter
#' @param n van Genuchten pore size distribution parameter
#' @param Ksat saturated hydraulic conductivity (kg s / m^3)
#' @return hydraulic conductivity (kg s / m^3)
#' @export
hydraulicConductivityFromThetaVG <- function(theta, Smin, Smax, psie, alpha, n, Ksat) {
  psie <- -abs(psie)
  # Calculate degree of saturation
  m <- 1 - (1 / n)
  Se <- (theta - Smin) / (Smax - Smin)
  Se <- pmin(1, pmax(0, Se))
  # Calculate water saturation at the air-entry potential
  Sc <- (1 + (alpha * abs(psie))^n)^(-m)
  # Calculate hydrualic conductivity
  top <- 1 - (1 - (Se * Sc)^(1 / m))^m
  btm <- 1 - (1 - Sc^(1/m))^m
  k <- (Ksat * Se^0.5 * (top / btm)^2)
  return(k)
}
#' @title Updates soil water status using the soil water model
#' @description Solves vertical soil water redistribution through the soil profile
#' over one time step, updating soil water potential, volumetric water content,
#' vapour content, and hydraulic conductivity.
#' @param soilwatermod list of variables provided as input and updated by the soil
#' water model, as returned by [InitSoilwatermod()] on the first step or by a
#' previous model run
#' @param soilp list of soil parameters as returned by [createsoilplist()]
#' @param climdata list of driving variables for the current time step. Must contain:
#' \describe{
#'   \item{Tair}{air temperature (degrees C)}
#'   \item{relhum}{air relative humidity (Percentage)}
#'   \item{rHa}{aerodynamic resistance to heat from the ground surface to the
#'   reference height (s m^-1)}
#'   \item{precip}{precipitation over the time step, expressed as an equivalent
#'   water flux (mm over timestep)}
#'   \item{Et}{plant transpiration demand (mm over timestep)}
#' }
#' @param dT time step (s)
#' @param pTAW fraction of available water below saturation at which soil begins
#' to limit transpiration in the root zone
#' @param maxNrIterations maximum number of iterations
#' @param tolerance convergence criterion
#' @return a list comprising
#' \describe{
#'   \item{n}{number of soil layers}
#'   \item{z}{depth of soil nodes (m)}
#'   \item{dz}{soil layer thicknesses (m)}
#'   \item{zCenter}{depth of layer centres (m)}
#'   \item{Tc}{soil temperature profile (degrees C)}
#'   \item{theta}{updated volumetric soil water content (m^3 m^-3)}
#'   \item{vol}{soil layer volumes per unit ground area (m)}
#'   \item{psiw}{updated soil water potential (J kg^-1)}
#'   \item{k}{total soil water conductivity including liquid and vapour components}
#'   \item{vapor}{updated soil vapour content}
#'   \item{oldvapor}{soil vapour content at previous time step}
#'   \item{oldtheta}{soil water content at previous time step (m^3 m^-3)}
#'   \item{rootfrac}{fraction of roots in each soil layer}
#'   \item{success}{logical indicating whether the solver converged within the
#'   specified tolerance}
#'   \item{iterations}{number of iterations used}
#'   \item{Evapmmhr}{soil evaporation over the timestep (mm)}
#' }
#' @export
SoilWaterModel <- function(soilwatermod, soilp, climdata, dT = 3600, pTAW = 0.5,
                           maxNrIterations = 20, tolerance = 1e-4,
                           useDamping = TRUE) {
  g <- 9.80665
  rho <- 1000
  n <- soilp$nLayers
  oldtheta <- soilwatermod$oldtheta
  oldvapor <- soilwatermod$oldvapor
  psiw <- soilwatermod$psiw
  theta <- soilwatermod$theta
  vapor <- soilwatermod$vapor
  k <- numeric(n + 1)
  aa <- bb <- cc <- dd <- ff <- u <- du <- Ca <- dpsi <- numeric(n)
  iter <- 1
  massBalance <- 1
  # keep surface flux external/fixed for this solve
  surfaceFlux <- (climdata$Evapmmhr - climdata$precip) / dT
  while (massBalance > tolerance && iter < maxNrIterations) {
    # transpiration from current iterate
    STr <- .transpiration_distribute(
      soilp, soilwatermod$rootfrac, climdata$Et, dT, psiw, pTAW
    )
    # suppress root uptake from frozen layers
    STr[soilwatermod$Tc <= 0] <- 0
    # bottom boundary for psi/theta
    if (soilp$FreeDrain == FALSE) {
      psiw[n + 1] <- soilp$psie[n]
      theta[n + 1] <- soilp$Smax[n]
      k[n + 1] <- soilp$Ksat[n]
    } else {
      psiw[n + 1] <- psiw[n]
      theta[n + 1] <- theta[n]
    }
    idx <- seq_len(n)
    kh <- hydraulicConductivityFromTheta(
      theta[idx], soilp$b[idx], soilp$Smax[idx], soilp$Ksat[idx]
    )
    kv <- vaporConductivityFromPsiTheta(
      psiw[idx], theta[idx], soilwatermod$Tc[idx] + 273.15, soilp$Smax[idx]
    )
    # suppress liquid flow in frozen soil
    frozen <- soilwatermod$Tc[idx] <= 0
    kh[frozen] <- 0
    k[idx] <- kh + kv
    # free drainage ghost conductivity must use CURRENT bottom-layer k
    if (!soilp$FreeDrain == FALSE) {
      k[n + 1] <- k[n]
    }
    u[idx] <- g * k[idx]
    du[idx] <- -u[idx] * soilp$n[idx] / psiw[idx]
    Cw <- .dTheta_dPsi(
      psiw[idx], soilp$psie[idx], soilp$b[idx], soilp$Smax[idx]
    )
    Cv <- .dvapor_dPsi(
      psiw[idx], soilp$psie[idx], theta[idx], soilp$b[idx],
      soilp$Smax[idx], soilwatermod$Tc[idx] + 273.15
    )
    Ca[idx] <- soilwatermod$vol[idx] * (rho * Cw + Cv) / dT
    ff[idx] <- ((psiw[idx + 1] * k[idx + 1] - psiw[idx] * k[idx]) /
                  (soilwatermod$dz[idx] * (1 - soilp$n[idx]))) - u[idx]
    massBalance <- 0
    # top layer
    aa[1] <- 0
    cc[1] <- -k[2] / soilwatermod$dz[1]
    bb[1] <- k[1] / soilwatermod$dz[1] + Ca[1] + du[1]
    dd[1] <- surfaceFlux + STr[1] - ff[1] +
      soilwatermod$vol[1] * (rho * (theta[1] - oldtheta[1]) +
                               (vapor[1] - oldvapor[1])) / dT
    massBalance <- massBalance + abs(dd[1])
    # remaining layers
    if (n > 1) {
      idx2 <- 2:n
      aa[idx2] <- -k[idx2 - 1] / soilwatermod$dz[idx2 - 1] - du[idx2 - 1]
      cc[idx2] <- -k[idx2 + 1] / soilwatermod$dz[idx2]
      bb[idx2] <- k[idx2] / soilwatermod$dz[idx2 - 1] +
        k[idx2] / soilwatermod$dz[idx2] + Ca[idx2] + du[idx2]
      dd[idx2] <- ff[idx2 - 1] + STr[idx2] - ff[idx2] +
        soilwatermod$vol[idx2] *
        (rho * (theta[idx2] - oldtheta[idx2]) +
           (vapor[idx2] - oldvapor[idx2])) / dT
      massBalance <- massBalance + sum(abs(dd[idx2]))
    }
    TBC <- ThomasBoundaryCondition(aa, bb, cc, dd, dpsi, 1, n)
    dpsi <- TBC$x
    if (useDamping) {
      dryness <- (psiw[idx] - soilp$psie[idx]) /
        (soilp$psi_min[idx] - soilp$psie[idx])
      dryness <- pmax(0, pmin(1, dryness))
      lambda <- 0.2 + 0.8 * dryness
    } else {
      lambda <- rep(1, n)
    }
    psiw[idx] <- psiw[idx] - lambda * dpsi[idx]
    psiw[idx] <- pmin(psiw[idx], soilp$psie[idx] - 1e-8)
    psiw[idx] <- pmax(psiw[idx], soilp$psi_min[idx])
    theta[idx] <- thetaFromPsi(
      psiw[idx], soilp$psie[idx], soilp$b[idx], soilp$Smax[idx]
    )
    vapor[idx] <- vaporFromPsi(
      psiw[idx], theta[idx], soilwatermod$Tc[idx] + 273.15, soilp$Smax[idx]
    )
    iter <- iter + 1
  }
  soilwatermod$psiw <- psiw
  soilwatermod$theta <- theta
  soilwatermod$vapor <- vapor
  soilwatermod$k <- k
  soilwatermod$Evapmmhr <- climdata$Evapmmhr
  soilwatermod$success <- (massBalance < tolerance)
  soilwatermod$iterations <- iter
  soilwatermod$error <- massBalance
  return(soilwatermod)
}
