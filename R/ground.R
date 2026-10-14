# Depths of the soil nodes (m): thin layers near the surface, thickening with
# depth. The deepest node, the lower boundary, is at the soil depth given, and one
# further depth is returned below it
.geometric <- function(n, totalDepth) {
  i <- seq_len(n + 1)
  dz_unit <- totalDepth / sum(seq_len(n)^2)
  # The top layer is one unit thick and each layer below it its number squared
  thick <- dz_unit * c(1, i[-1]^2)
  z <- c(0, cumsum(thick))
  return(z)
}
# Slope of the water retention curve: change in volumetric water content per
# unit change in water potential. Zero in saturated soil
.dTheta_dPsi <- function(psiw, psie, b, Smax) {
  psie <- -abs(psie)
  Se <- (psiw / psie)^(-1 / b)
  theta <- Smax * Se
  out <- -theta / (b * psiw)
  out[psiw >= psie] <- 0
  return(out)
}
# Density of water vapour in saturated air (kg/m^3) at temperature Tk (K), from
# saturated vapour pressure by the ideal gas law
.satVapDensity <- function(Tk) {
  tc <- Tk - 273.15
  es_Pa <- 1000 * satvap(tc)  # kPa -> Pa
  0.018015 * es_Pa / (8.314 * Tk)
}
# Change in the water vapour held in the air-filled pores per unit change in
# water potential: drier soil has more pore air, but that air is less humid
.dvapor_dPsi <- function(psiw, psie, theta, b, Smax, Tk) {
  humidity <- exp(0.018015 * psiw / (8.314 * Tk))
  dTdP <- .dTheta_dPsi(psiw, psie, b, Smax)
  vp <- .satVapDensity(Tk)
  capacity_vapor <- (Smax - theta) * vp * humidity *
    (0.018015 / (8.314 * Tk)) - dTdP * vp * humidity
  return(capacity_vapor)
}
# Fraction of roots in the soil held by each node, from the depths above which
# half and 95 percent of roots lie
.root_distribute <- function(z, root50, root95) {
  n <- length(z) - 2
  # Fraction of roots above a given depth rises along an S-shaped curve, which
  # passes through a half at root50 and 0.95 at root95
  p <- log(19) / log(root95 / root50)
  zmid <- 0.5 * (z[1:n] + z[2:(n + 1)]) # depth of the base of the soil held by each node
  x <- (zmid / root50)^p
  above <- x / (1 + x)
  above[n] <- 1 # the deepest rooted node holds all roots below it
  # Fraction of roots in the soil held by each node, with none at the lower boundary
  v <- c(diff(c(0, above)), 0)
  return(v)
}
# Shares transpiration among soil nodes by root fraction, reduced where the soil
# is too wet or too dry for roots to take up water. Returns the rate at which
# water is taken from each node (kg m^-2 s^-1)
.transpiration_distribute <- function(soilp, rootfrac, totalTransp_mm, dT, psiw, p = 0.5) {
  # Uptake in wet soil: none at saturation, rising to full once air enters the pores
  .alpha_wet <- function(psiw, psie) {
    ifelse(psiw >= 0, 0, ifelse(psiw <= psie, 1, abs(psiw) / abs(psie)))
  }
  # Uptake in dry soil: full down to psi_dry, falling to none at the wilting point
  .alpha_dry <- function(psiw, psi_dry, psi_wilt) {
    ifelse(psiw <= psi_wilt, 0, ifelse(psiw >= psi_dry, 1,
                                       (psiw - psi_wilt) / (psi_dry - psi_wilt)))
  }
  n <- length(psiw)
  S <- numeric(n)
  if (totalTransp_mm > 0) {
    Trate <- totalTransp_mm / dT
    idx <- seq_len(n)
    # Roots can take up the water held between field capacity (-33 J/kg) and
    # the wilting point (-1500 J/kg)
    psi_fc <- -33
    psi_wilt <- -1500
    theta_fc <- thetaFromPsi(psi_fc, soilp$psie[idx], soilp$b[idx], soilp$Smax[idx])
    theta_wilt <- thetaFromPsi(psi_wilt, soilp$psie[idx], soilp$b[idx], soilp$Smax[idx])
    # Uptake starts to fall once the fraction p of that water has been used
    theta_dry <- theta_fc - p * (theta_fc - theta_wilt)
    psi_dry <- PsiFromtheta(
      theta = theta_dry,
      psie = soilp$psie[idx],
      b = soilp$b[idx],
      Smax = soilp$Smax[idx],
      MPa = FALSE
    )
    aw <- .alpha_wet(psiw, soilp$psie[idx])
    ad <- .alpha_dry(psiw, psi_dry, psi_wilt)
    w <- rootfrac * aw * ad
    # If no rooted node can supply water, uptake is shared by root fraction alone
    if (sum(w) <= 0) w <- rootfrac
    if (sum(w) > 0) w <- w / sum(w)
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
#' @param totalDepth Numeric. Depth of the lower boundary of the soil profile (m),
#'   where the deepest node lies.
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
#'   \item{Ksat}{Saturated hydraulic conductivity (kg s / m^3), repeated for each soil layer}
#'   \item{Vq}{Quartz fraction by volume (m^3 m^-3), repeated for each soil layer}
#'   \item{Vm}{Mineral fraction by volume (m^3 m^-3), repeated for each soil layer}
#'   \item{Vo}{Organic fraction by volume (m^3 m^-3), repeated for each soil layer and adjusted with depth}
#'   \item{Mc}{Mass fraction of clay (kg kg^-1), repeated for each soil layer}
#'   \item{b}{Campbell soil water retention parameter (unitless), repeated for each soil layer.
#'     The Campbell hydraulic conductivity exponent used internally by the soil water solver
#'     is derived from this as \code{2 * b + 3}; it is not a separate element of this list.}
#'   \item{psi_e}{Air-entry water potential (J/kg, numerically equal to kPa), repeated for each soil layer}
#'   \item{gref}{Ground shortwave reflectance (unitless)}
#'   \item{groundem}{Ground longwave emissivity (unitless)}
#'   \item{grefPAR}{Ground reflectance for photosynthetically active radiation (unitless)}
#'   \item{slope}{Surface slope (degrees)}
#'   \item{aspect}{Surface aspect (degrees)}
#'   \item{nLayers}{Number of soil layers}
#'   \item{totalDepth}{Depth of the lower boundary of the soil profile (m)}
#'   \item{FreeDrain}{Logical indicating whether the lower boundary is free draining}
#'   \item{psi_min}{Water potential of oven-dry soil (J/kg), the driest the soil is allowed to get, repeated for each soil layer}
#'   \item{psie}{Air-entry water potential as a negative value (J/kg), repeated for each soil layer}
#' }
#' @export
createsoilplist <- function(soiltype = "Clay loam", nlayers = 15, totalDepth = 3, slope = 0, aspect = 180,
                            surface_organicmu = 3, FreeDrain = TRUE) {
  # Typical physical and hydraulic properties of the soil type, the same at
  # every node to begin with
  soilparams <- micropoint::newsoilparamstable
  s <- which(soilparams$Soil.type == soiltype)
  campbellcols <- c("Smax", "Smin", "Ksat", "Vq", "Vm", "Vo", "Mc", "b", "psi_e")
  v <- as.numeric(soilparams[s, campbellcols])
  names(v) <- campbellcols
  soilc <- as.list(v)
  for (nm in campbellcols) {
    soilc[[nm]] <- rep(v[[nm]], nlayers + 1)
  }
  # Organic matter is concentrated near the surface: its fraction is multiplied
  # by surface_organicmu at the surface, falling to nothing at the lower boundary
  mu <- rev(.geometric(nlayers, surface_organicmu))[2:(nlayers + 2)]
  soilc$Vo <- soilc$Vo * mu
  # Quartz, other minerals and clay give way to the organic matter gained, or
  # take the place of that lost
  dif <- soilc$Vo - v[["Vo"]]
  sm <- v[["Vq"]] + v[["Vm"]]
  mu <- (sm - dif) / sm
  soilc$Vm <- soilc$Vm * mu
  soilc$Vq <- soilc$Vq * mu
  soilc$Mc <- soilc$Mc * mu
  # Radiative properties of the soil surface, and its slope and aspect
  soilc$gref <- 0.15
  soilc$groundem <- 0.97
  soilc$grefPAR <- 0.75 * soilc$gref
  soilc$slope <- slope
  soilc$aspect <- aspect
  soilc$nLayers <- nlayers
  soilc$totalDepth <- totalDepth
  soilc$FreeDrain <- FreeDrain
  # Driest the soil is allowed to get: the water potential of oven-dry soil (J/kg)
  soilc$psi_min <- rep(-1e6, nlayers + 1)
  # Air-entry water potential as a negative value
  soilc$psie <- -abs(soilc$psi_e)
  return(soilc)
}
#' @title Initializes input for the soil heat model
#' @description Creates the soil state object used by the soil heat model. The
#' soil is represented by a column of nodes, closely spaced near the surface
#' and further apart with depth. The first node is at the soil surface and the
#' last, the lower boundary, is at the depth given to [createsoilplist()].
#' @param soilp list of soil parameters as returned by [createsoilplist()]
#' @param Te vector of the temperatures of the soil nodes (degrees C), from the
#' surface to the lower boundary: one more than the number of soil layers. The
#' temperature of the lower boundary stays fixed.
#' @param wc vector of the volumetric water contents of the soil nodes (m3 m-3),
#' of the same length as `Te`
#' @return a list comprising
#' \describe{
#'   \item{n}{number of soil layers}
#'   \item{z}{depth of soil nodes (m), with one further depth below the lower boundary}
#'   \item{dz}{distance between each node and the node beneath it (m)}
#'   \item{vol}{thickness of the soil for which each node stores heat (m)}
#'   \item{wc}{volumetric water content of the soil at each node (m3 m-3)}
#'   \item{Te}{temperature of the soil at each node (degrees C)}
#'   \item{oldTe}{temperature of the soil at each node at the start of the time step (degrees C)}
#'   \item{Gflux}{heat flux into the soil at its surface (W m-2)}
#'   \item{iters}{number of solver iterations}
#' }
#' @export
InitSoilheatmod <- function(soilp, Te, wc) {
  n <- soilp$nLayers
  # Depths of the nodes, and the distance from each node to the one beneath it
  z <- .geometric(n, soilp$totalDepth)
  dz <- z[2:(n + 2)] - z[1:(n + 1)]
  # Each node stores heat for the soil around it: from halfway to the node above
  # to halfway to the node below. The surface node has only the soil beneath it
  vol <- c(0.5 * dz[1], 0.5 * (z[3:(n + 2)] - z[1:n]))
  state <- list(n = n, z = z, dz = dz, vol = vol,
                wc = wc, oldTe = Te, Te = Te, Gflux = 0, iters = 0)
  return(state)
}

#' @title Initializes input for the soil water model
#' @description Creates the soil state object used by the soil water model from
#' an existing soil heat model state
#' @param soilp list of soil parameters as returned by [createsoilplist()]
#' @param soilheatmod soil heat model state as returned by [InitSoilheatmod()]
#' @param root50 depth above which half of all roots lie (m)
#' @param root95 depth above which 95 percent of all roots lie (m)
#' @return a list comprising
#' \describe{
#'   \item{n}{number of soil layers}
#'   \item{z}{depth of soil nodes (m), with one further depth below the lower boundary}
#'   \item{dz}{distance between each node and the node beneath it (m)}
#'   \item{Tc}{temperature of the soil at each node (degrees C)}
#'   \item{oldTc}{temperature of the soil at each node at the start of the time step (degrees C)}
#'   \item{theta}{volumetric water content of the soil at each node (m3 m-3)}
#'   \item{vol}{thickness of the soil for which each node stores water (m)}
#'   \item{psiw}{soil water potential (J kg-1)}
#'   \item{k}{hydraulic conductivity (kg s m-3)}
#'   \item{vapor}{water vapour held in the air-filled pores (kg m-3)}
#'   \item{oldvapor}{water vapour held in the air-filled pores at the previous time step (kg m-3)}
#'   \item{oldtheta}{soil water content at previous time step (m3 m-3)}
#'   \item{rootfrac}{fraction of roots in the soil held by each node}
#'   \item{success}{logical indicating convergence success - set to TRUE on initalization}
#'   \item{iterations}{number of solver iterations - set to zero initially}
#' }
#' @export
InitSoilwatermod <- function(soilp, soilheatmod, root50, root95) {
  out <- list()
  out$n = soilheatmod$n # number of layers
  out$z = soilheatmod$z # node depths
  out$dz = soilheatmod$dz # distance between each node and the one beneath
  out$Tc = soilheatmod$Te # temperature of soil profile
  out$oldTc = soilheatmod$oldTe # old temperature of soil profile
  out$theta = soilheatmod$wc # volumetric water content of soil profile
  # Each node stores water for the same soil as it stores heat for: half the
  # distance to each neighbour, and for the surface node half the top layer
  out$vol <- soilheatmod$vol
  out$psiw <- PsiFromtheta(out$theta, soilp$psie, soilp$b, soilp$Smax, MPa = FALSE)
  out$k <- hydraulicConductivityFromTheta(out$theta, soilp$b, soilp$Smax, soilp$Ksat)
  out$vapor <- vaporFromPsi(out$psiw, out$theta, out$Tc + 273.15, soilp$Smax)
  out$oldvapor <- out$vapor
  out$oldtheta <- soilheatmod$wc
  out$rootfrac <- .root_distribute(out$z, root50, root95)
  out$success <- TRUE
  out$iterations <- 0
  return(out)
}
#' @title Calculates soil surface relative humidity
#' @description Calculates the relative humidity of the air in the pores of the
#' soil, which falls as the soil dries and holds its water more tightly.
#' @param theta volumetric soil water content (m^3 / m^3)
#' @param tsoil soil temperature (degrees C)
#' @param psie air-entry water potential (J / kg)
#' @param b Campbell soil water retention parameter
#' @param Smax saturated volumetric water content (m^3 / m^3)
#' @return relative humidity (Percentage)
#' @export
soilrelhum <- function(theta, tsoil, psie, b, Smax) {
  # Water potential of the soil (J/kg), from the water retention curve
  psie <- -abs(psie)
  psiw <- psie * (theta / Smax)^-b
  # Humidity of air in equilibrium with water held at that potential
  Tk <- tsoil + 273.15
  hr <- exp(0.018015 * psiw / (8.314 * Tk)) * 100
  return(hr)
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
  # Conductivity of the solids (W/m/K): quartz conducts heat far better than
  # other minerals (8.8 against 2.5), and organic matter poorly (0.25)
  Vsolid <- Vq + Vm + Vo
  kSolid <- 8.8^(Vq / Vsolid) * 2.5^(Vm / Vsolid) * 0.25^(Vo / Vsolid)
  # Shape of the soil particles, which sets how readily heat passes between them
  Ga <- 0.088 * ((Vq + Vm) / Vsolid) + 0.5 * (Vo / Vsolid)
  Gc <- 1 - 2 * Ga
  # How wet the soil must be for water vapour to move freely through its pores,
  # which depends on its clay content
  Q <- 7.25 * Mc + 2.52
  xwo <- 0.33 * Mc + 0.078
  # Space between the solids, and the part of it that is filled with air
  porosity <- 1 - Vsolid
  gasPorosity <- pmax(porosity - Vw, 0)
  # Conductivity of water, and of the air in the pores. In moist soil, water
  # evaporates on the warm side of a pore and condenses on the cool side, which
  # carries heat across it and adds to the conductivity of the air
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
  # Conductivity of the soil as a whole: that of water, air and solids, each
  # weighted by its volume and by how it is arranged around the particles
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
  # Air fills whatever volume the solids and water leave
  Vsum <- Vq + Vm + Vo + Vw
  Va <- ifelse(Vsum > 1, 0, 1 - Vsum)
  Vq <- ifelse(Vsum > 1, Vq / Vsum, Vq)
  Vm <- ifelse(Vsum > 1, Vm / Vsum, Vm)
  Vo <- ifelse(Vsum > 1, Vo / Vsum, Vo)
  Vw <- ifelse(Vsum > 1, Vw / Vsum, Vw)
  # Heat capacity of each constituent (MJ m-3 K-1). Water below freezing is
  # taken to be ice
  CHa <- cpair(Tc) * rhoair(Tc, pk) / 1e6
  if (length(Tc) == 1L) Tc <- rep(Tc, length(Vw))
  ChW <- 1.93 + 0.0067 * Tc
  ChW[Tc >= 0] <- 4.18
  # That of the soil is the sum of these, each weighted by its volume
  CH <- Vq * 2.13 + Vm * 2.31 + Vo * 2.50 + Va * CHa + Vw * ChW
  return(CH * 1e6)
}
#' @title Calculates soil surface energy balance
#' @description Calculates the energy available at the soil surface (W m-2):
#' radiation absorbed, less that emitted and the sensible and latent heat given
#' to the air. It is the heat flux into the soil.
#' @param soilp list of soil parameters as returned by [createsoilplist()]
#' @param Rabs absorbed radiation at the surface (W m-2)
#' @param Tref reference air temperature (degrees C)
#' @param Tsurface soil surface temperature (degrees C)
#' @param pk air pressure (kPa)
#' @param relhum air relative humidity (Percentage)
#' @param rHa aerodynamic resistance to heat from ground to reference height
#' as returned by [groundresistance()] (s m-1)
#' @param theta volumetric water content of the soil at its surface (m3 m-3)
#' @return heat flux into the soil (W m-2)
#' @export
soilsurfaceEB <- function(soilp, Rabs, Tref, Tsurface, pk, relhum, rHa, theta) {
  # net radiation
  sb <- 5.67e-8
  Rnet <- Rabs - soilp$groundem * sb * (Tsurface + 273.15)^4
  # sensible heat
  cp <- cpair(Tref)
  ph <- rhoair(Tref, pk)
  H <- (ph * cp / rHa) * (Tsurface - Tref)
  # latent heat
  hr <- soilrelhum(theta, Tsurface, soilp$psie[1], soilp$b[1], soilp$Smax[1])
  es <- satvap(Tsurface) * (hr / 100)
  ea <- satvap(Tref) * (relhum / 100)
  la <- latvap(Tsurface)
  L <- (la * ph / (rHa * pk)) * (es - ea)
  return(Rnet - H - L)
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
  psiw <- PsiFromtheta(theta, psie, b, Smax, MPa = FALSE) # water potential of the soil (J/kg)
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
  # Working down the rows, remove from each the unknown it shares with the row above
  for (i in seq(first, length.out = last - first)) {
    cc[i] <- cc[i] / bb[i]
    dd[i] <- dd[i] / bb[i]
    bb[i + 1] <- bb[i + 1] - aa[i + 1] * cc[i]
    dd[i + 1] <- dd[i + 1] - aa[i + 1] * dd[i]
  }
  # The last row then holds one unknown. Working back up, each row gives the next
  x[last] <- dd[last] / bb[last]
  for (i in rev(seq(first, length.out = last - first))) {
    x[i] <- dd[i] - cc[i] * x[i + 1]
  }
  return(list(bb = bb, cc = cc, dd = dd, x = x))
}
#' @title Updates soil temperatures using the soil heat model
#' @description Solves heat conduction through the soil profile over one time step.
#' The energy balance of the soil surface, from [soilsurfaceEB()], is the upper
#' boundary, and the deepest node is held at a fixed temperature. Heat is
#' conducted between nodes according to [thermalConductivity()] and stored
#' according to [heatCapacity()], with temperatures at the end of the time step
#' found using [ThomasBoundaryCondition()].
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
#' @param maxNrIterations maximum number of iterations
#' @param tolerance convergence criterion for maximum change in node temperature
#' between successive iterations (degrees C)
#' @details The model steps forward from the temperatures in `oldTe`, which are
#' those at the start of the time step, and returns those at its end in `Te`.
#' `oldTe` is left as it is, so that the same time step can be solved again with
#' different inputs. To move on to the next time step, set `oldTe` to `Te`.
#' @return a list comprising
#' \describe{
#'   \item{n}{number of soil layers}
#'   \item{z}{depth of soil nodes (m), with one further depth below the lower boundary}
#'   \item{dz}{distance between each node and the node beneath it (m)}
#'   \item{vol}{thickness of the soil for which each node stores heat (m)}
#'   \item{wc}{volumetric water content of the soil at each node (m3 m-3)}
#'   \item{Te}{temperature of the soil at each node at the end of the time step (degrees C)}
#'   \item{oldTe}{temperature of the soil at each node at the start of the time step (degrees C)}
#'   \item{Gflux}{heat flux into the soil at its surface (W m-2)}
#'   \item{iters}{number of iterations used to reach convergence}
#' }
#' @export
SoilHeatModel <- function(soilheatmod, soilp, Rabs, Tref, relhum, atmPressure, rHa, dT = 3600,
                          maxNrIterations = 100, tolerance = 1e-3) {
  n <- soilheatmod$n
  wc <- soilheatmod$wc
  dz <- soilheatmod$dz
  oldTe <- soilheatmod$oldTe
  # Temperatures at the end of the time step are found by iteration, starting
  # from the latest estimate. The deepest node keeps its temperature
  Te <- soilheatmod$Te
  Te[n + 1] <- oldTe[n + 1]
  idx <- seq_len(n) # the nodes solved

  nrIterations <- 0
  maxdT <- Inf
  while (maxdT > tolerance && nrIterations < maxNrIterations) {
    # Thermal properties of the soil at its current temperature. C is the heat
    # needed to warm the soil held by each node by one degree (J m-2 K-1)
    lambda <- thermalConductivity(soilp$Vq, soilp$Vm, soilp$Vo, wc, soilp$Mc, Te, atmPressure)
    C <- heatCapacity(soilp$Vq, soilp$Vm, soilp$Vo, wc, Te, atmPressure) * soilheatmod$vol

    # Conductance between each node and the one beneath it (W m-2 K-1), from the
    # conductivity of the upper node. Close to the surface conductivity changes
    # quickly with depth, so the logarithmic mean of the top two nodes is used
    ff <- lambda[idx] / dz[idx]
    if (lambda[1] != lambda[2]) {
      ff[1] <- ((lambda[1] - lambda[2]) / log(lambda[1] / lambda[2])) / dz[1]
    }
    above <- c(0, ff[-n]) # conductance between each node and the one above it

    # Heat flux into the soil at its surface, and how it changes as the surface
    # warms: a warmer surface emits more radiation and gives more heat to the air
    qsurface <- soilsurfaceEB(soilp, Rabs, Tref, Te[1], atmPressure, relhum, rHa, wc[1])
    dqsurface <- (soilsurfaceEB(soilp, Rabs, Tref, Te[1] + 0.01, atmPressure, relhum, rHa, wc[1]) - qsurface) / 0.01

    # Energy balance of each node over the time step: the heat it stores equals
    # the heat conducted to it from the nodes above and below, at the
    # temperatures of the end of the time step
    aa <- -above
    cc <- -ff
    bb <- C[idx] / dT + above + ff
    dd <- C[idx] / dT * oldTe[idx]
    # The surface node also gains the heat flux into the soil
    bb[1] <- bb[1] - dqsurface
    dd[1] <- dd[1] + qsurface - dqsurface * Te[1]
    # The node above the lower boundary exchanges heat with a fixed temperature
    cc[n] <- 0
    dd[n] <- dd[n] + ff[n] * Te[n + 1]

    Tnew <- ThomasBoundaryCondition(aa, bb, cc, dd, numeric(n))$x
    maxdT <- max(abs(Tnew - Te[idx]))
    Te[idx] <- Tnew
    nrIterations <- nrIterations + 1
  }
  soilheatmod$Te <- Te
  soilheatmod$Gflux <- soilsurfaceEB(soilp, Rabs, Tref, Te[1], atmPressure, relhum, rHa, wc[1])
  soilheatmod$iters <- nrIterations
  return(soilheatmod)
}
#' Create list of Campbell soil hydraulic parameters
#' @param soiltype Character string giving the soil type. See details
#' @return A list of soil parameters:
#' \describe{
#'   \item{psie}{Air-entry water potential (J/kg, numerically equal to kPa), given as a positive value}
#'   \item{b}{Campbell soil water retention parameter (unitless)}
#'   \item{Smin}{Residual volumetric soil water content (m^3 m^-3)}
#'   \item{Smax}{Volumetric soil water content at saturation (m^3 m^-3)}
#'   \item{Ksat}{Saturated hydraulic conductivity (kg s / m^3)}
#' }
#' @details
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
#'   \item{alpha}{Van Genuchten shape parameter (kg/J): the inverse of a water potential in J/kg}
#'   \item{n}{Pore size distribution parameter (unitless)}
#'   \item{Smin}{Residual volumetric soil water content (m^3 m^-3)}
#'   \item{Smax}{Volumetric soil water content at saturation (m^3 m^-3)}
#'   \item{Ko}{Matching point hydraulic conductivity of the unsaturated conductivity curve (kg s / m^3)}
#'   \item{L}{Pore connectivity parameter of the unsaturated conductivity curve (unitless)}
#' }
#' @references Schaap MG, Leij FJ (2000) Improved prediction of unsaturated
#' hydraulic conductivity with the Mualem-van Genuchten model. Soil Science
#' Society of America Journal 64: 843-851.
#'
#' Schaap MG, Leij FJ, van Genuchten MTh (2001) Rosetta: a computer program for
#' estimating soil hydraulic parameters with hierarchical pedotransfer
#' functions. Journal of Hydrology 251: 163-176.
#' @details
#' The values are the class averages of the Rosetta pedotransfer functions.
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
  # Parameters of the water retention and unsaturated conductivity curves,
  # averaged over the soils of each texture class
  vg <- data.frame(
    Soil.type = c("Sand", "Loamy sand", "Sandy loam", "Loam", "Silt loam", "Sandy clay loam",
                  "Clay loam", "Silty clay loam", "Sandy clay", "Silty clay", "Clay", "Silt"),
    alpha = c(0.359, 0.354, 0.272, 0.113, 0.052, 0.215, 0.161, 0.086, 0.341, 0.165, 0.153, 0.067),
    n     = c(3.177, 1.746, 1.449, 1.472, 1.663, 1.330, 1.416, 1.521, 1.208, 1.321, 1.253, 1.679),
    Smin  = c(0.053, 0.049, 0.039, 0.061, 0.065, 0.063, 0.079, 0.090, 0.117, 0.111, 0.098, 0.050),
    Smax  = c(0.375, 0.390, 0.387, 0.399, 0.439, 0.384, 0.442, 0.482, 0.385, 0.481, 0.459, 0.489),
    Ko    = c(2.89e-4, 2.87e-4, 1.83e-4, 4.36e-5, 2.06e-5, 8.18e-5, 5.9e-5, 2.64e-5, 5.11e-5,
              3.74e-5, 3.5e-5, 3.94e-5),
    L     = c(-0.930, -0.874, -0.861, -0.371, 0.365, -1.280, -0.763, -0.156, -3.665, -1.287,
              -1.561, 0.624)
  )
  s <- which(vg$Soil.type == soiltype)
  out <- as.list(vg[s, c("alpha", "n", "Smin", "Smax", "Ko", "L")])
  return(out)
}
#' @title Derives soil water potential from volumetric water content using the Campbell method
#' @param theta volumetric soil water fraction (m^3 / m^3)
#' @param psie air entry water potential (J/kg)
#' @param b Campbell soil water retention parameter
#' @param Smax volumetric soil water fraction at saturation (m^3 / m^3)
#' @param MPa optional logical indicating whether to return outputs in MPa (TRUE) or kPa (FALSE)
#' @returns soil water potential (MPa or kPa ~ J/Kg)
#' @export
PsiFromtheta <- function(theta, psie, b, Smax, MPa = TRUE) {
  psie <- -abs(psie)
  # Degree of saturation: the fraction of the pore space filled with water
  Se <- theta / Smax
  Se[Se > 1] <- 1
  # Water is held more tightly as the soil dries, following a power law
  psiw <- psie * Se^(-b) # J/kg
  if (MPa) psiw <- psiw * 0.001 # Convert to MPa
  return(psiw)
}
#' @title Derives soil water potential from volumetric water content using the Van Genuchten method
#' @param theta volumetric soil water fraction (m^3 / m^3)
#' @param Smin Residual volumetric soil water fraction (m^3 / m^3)
#' @param Smax volumetric soil water fraction at saturation (m^3 / m^3)
#' @param alpha Van Genuchten shape parameter (kg/J)
#' @param n Van Genuchten pore size distribution parameter
#' @returns soil water potential (kPa ~ J/Kg)
#' @export
PsiFromthetaVG <- function(theta, Smin, Smax, alpha, n) {
  m <- 1 - (1 / n)
  # Degree of saturation: the fraction of the water between residual and
  # saturation that the soil holds
  Se <- (theta - Smin) / (Smax - Smin)
  Se <- pmin(1 - 1e-12, pmax(1e-12, Se))
  # Water potential, which falls smoothly from zero as the soil dries
  psiw <- -(1 / alpha) * (Se^(-1 / m) - 1)^(1 / n)
  return(psiw)
}
#' @title Calculates soil vapour content from water potential
#' @description Calculates the water vapour held in the air-filled pores of
#' the soil from soil water potential using the Kelvin equation
#' @param psiw soil water potential (J kg-1)
#' @param theta volumetric soil water content (m3 m-3)
#' @param Tk soil temperature (K)
#' @param Smax volumetric soil water content at saturation (m3 m-3)
#' @return soil vapour content (kg of vapour per m3 of soil)
#' @export
vaporFromPsi <- function(psiw, theta, Tk, Smax) {
  # Relative humidity of air in equilibrium with water held at that potential
  humidity <- exp(0.018015 * psiw / (8.314 * Tk))
  # Vapour in the part of the pore space that is filled with air
  vapor <- (Smax - theta) * .satVapDensity(Tk) * humidity
  return(vapor)
}
#' @title Calculates vapour conductivity from soil water potential
#' @description Calculates the conductivity of the soil for water vapour,
#' which diffuses through its air-filled pores, as a function of water
#' potential and water content
#' @param psiw soil water potential (J kg-1)
#' @param theta volumetric soil water content (m3 m-3)
#' @param Tk soil temperature (K)
#' @param Smax volumetric soil water content at saturation (m3 m-3)
#' @return soil vapour conductivity (kg s / m^3)
#' @export
vaporConductivityFromPsiTheta <- function(psiw, theta, Tk, Smax) {
  humidity <- exp(0.018015 * psiw / (8.314 * Tk))
  # Diffusivity of vapour in air (0.000024 m2/s), reduced to 0.66 of that by the
  # winding path through the air-filled pores
  k <- 0.66 * (Smax - theta) * 0.000024 * .satVapDensity(Tk) * humidity * 0.018015 / (8.314 * Tk)
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
  # Degree of saturation. The soil is saturated at potentials above air entry
  Se <- (pmin(psiw, psie) / psie)^(-1 / b)
  theta <- Se * Smax
  return(theta)
}
#' @title Calculates volumetric soil water fraction from water potential using the Van Genuchten method
#' @param psiw soil water potential (kPa ~ J/Kg)
#' @param Smin Residual volumetric soil water fraction (m^3 / m^3)
#' @param Smax volumetric soil water fraction at saturation (m^3 / m^3)
#' @param alpha Van Genuchten shape parameter (kg/J)
#' @param n Van Genuchten pore size distribution parameter
#' @returns volumetric soil water content (m3 m-3)
#' @export
thetaFromPsiVG <- function(psiw, Smin, Smax, alpha, n) {
  m <- 1 - (1 / n)
  # Degree of saturation: the fraction of the water between residual and
  # saturation that the soil holds
  Se <- (1 + (alpha * abs(psiw))^n)^(-m)
  theta <- Se * (Smax - Smin) + Smin
  return(theta)
}
#' @title Calculates hydraulic conductivity from soil water content using Campbell method
#' @param theta volumetric soil water content (m3 m-3)
#' @param b Campbell soil water retention parameter
#' @param Smax volumetric soil water content at saturation (m3 m-3)
#' @param Ksat saturated hydraulic conductivity (kg s / m^3)
#' @return hydraulic conductivity (kg s / m^3)
#' @export
hydraulicConductivityFromTheta <- function(theta, b, Smax, Ksat) {
  # Conductivity falls steeply as the soil dries and water is left only in the
  # smallest pores
  n <- 2 * b + 3
  k <- Ksat * (theta / Smax)^n
  return(k)
}
#' @title Calculates hydraulic conductivity from soil water content using Van Genuchten method
#' @param theta volumetric soil water content (m3 m-3)
#' @param Smin residual volumetric soil water fraction (m^3 / m^3)
#' @param Smax volumetric soil water fraction at saturation (m^3 / m^3)
#' @param n van Genuchten pore size distribution parameter
#' @param Ko matching point hydraulic conductivity (kg s / m^3): the conductivity
#' the curve gives at saturation
#' @param L pore connectivity parameter (unitless)
#' @return hydraulic conductivity (kg s / m^3)
#' @export
hydraulicConductivityFromThetaVG <- function(theta, Smin, Smax, n, Ko, L = 0.5) {
  m <- 1 - (1 / n)
  # Degree of saturation
  Se <- (theta - Smin) / (Smax - Smin)
  Se <- pmin(1, pmax(1e-12, Se))
  # Conductivity of the water-filled pores, which falls as the soil dries and
  # water is left only in the smallest and least connected of them
  k <- Ko * Se^L * (1 - (1 - Se^(1/m))^m)^2
  return(k)
}
#' @title Updates soil water status using the soil water model
#' @description Solves the movement of water through the soil profile over one
#' time step, updating soil water potential, volumetric water content, vapour
#' content and conductivity. Liquid water moves between nodes along gradients in
#' water potential and drains under gravity, and water vapour diffuses through
#' the air-filled pores. Rain enters and water evaporates at the surface, roots
#' take up the water transpired, and water drains from the base of the profile.
#' @param soilwatermod list of variables provided as input and updated by the soil
#' water model, as returned by [InitSoilwatermod()] on the first step or by a
#' previous model run
#' @param soilp list of soil parameters as returned by [createsoilplist()]
#' @param climdata list of driving variables for the current time step. Must contain:
#' \describe{
#'   \item{Tair}{air temperature (degrees C)}
#'   \item{relhum}{air relative humidity (Percentage)}
#'   \item{rHa}{aerodynamic resistance to heat from the ground surface to the
#'   reference height (s m^-1). Evaporation from the soil surface is calculated
#'   from this using [evaporation_flux()]; an infinite value gives none.}
#'   \item{precip}{precipitation over the time step (mm)}
#'   \item{Et}{plant transpiration over the time step (mm)}
#' }
#' @param dT time step (s)
#' @param pTAW fraction of the water held between field capacity and the
#' wilting point that roots can take up before uptake starts to fall
#' @param maxNrIterations maximum number of iterations
#' @param tolerance convergence criterion: the water unaccounted for, summed
#' over nodes (kg m^-2 s^-1)
#' @details The model steps forward from the water and vapour contents in
#' `oldtheta` and `oldvapor`, which are those at the start of the time step.
#' These are left as they are, so that the same time step can be solved again
#' with different inputs. To move on to the next time step, set them to `theta`
#' and `vapor`. Soil temperature (`Tc`) is taken from the soil heat model. Water
#' at or below 0 degrees C is taken to be frozen: it does not flow as liquid and
#' roots do not take it up.
#' @return a list comprising
#' \describe{
#'   \item{n}{number of soil layers}
#'   \item{z}{depth of soil nodes (m), with one further depth below the lower boundary}
#'   \item{dz}{distance between each node and the node beneath it (m)}
#'   \item{Tc}{temperature of the soil at each node (degrees C)}
#'   \item{theta}{updated volumetric soil water content (m^3 m^-3)}
#'   \item{vol}{thickness of the soil for which each node stores water (m)}
#'   \item{psiw}{updated soil water potential (J kg^-1)}
#'   \item{k}{conductivity of the soil for liquid water and vapour together (kg s m^-3)}
#'   \item{vapor}{updated water vapour held in the air-filled pores (kg m^-3)}
#'   \item{oldvapor}{vapour held in the air-filled pores at the start of the time step (kg m^-3)}
#'   \item{oldtheta}{soil water content at the start of the time step (m^3 m^-3)}
#'   \item{rootfrac}{fraction of roots in the soil held by each node}
#'   \item{success}{logical indicating whether the solver converged within the
#'   specified tolerance}
#'   \item{iterations}{number of iterations used}
#'   \item{Evapmmhr}{evaporation from the soil surface over the time step (mm)}
#'   \item{surplus}{water that a saturated soil surface could not take in over
#'   the time step, which runs off (mm)}
#'   \item{drainage}{water leaving through the base of the profile over the time
#'   step (mm). Negative where water enters.}
#'   \item{error}{water unaccounted for when the solver stopped (kg m^-2 s^-1)}
#' }
#' @export
SoilWaterModel <- function(soilwatermod, soilp, climdata, dT = 3600, pTAW = 0.5,
                           maxNrIterations = 100, tolerance = 1e-6) {
  g <- 9.80665     # gravity (m/s^2)
  rho <- 1000      # density of water (kg/m^3)
  MwR <- 0.018015 / 8.314 # molar mass of water over the gas constant
  n <- soilp$nLayers
  idx <- seq_len(n) # the nodes solved. Node n + 1 is the lower boundary
  dz <- soilwatermod$dz[idx]
  vol <- soilwatermod$vol[idx]
  Tk <- soilwatermod$Tc[idx] + 273.15
  frozen <- soilwatermod$Tc[idx] <= 0
  psie <- soilp$psie[idx]
  b <- soilp$b[idx]
  Smax <- soilp$Smax[idx]
  Ksat <- soilp$Ksat[idx]
  psi_min <- soilp$psi_min[idx]
  rootfrac <- soilwatermod$rootfrac[idx]
  # Water and vapour held at the start of the time step
  oldtheta <- soilwatermod$oldtheta[idx]
  oldvapor <- soilwatermod$oldvapor[idx]
  # Water potential at the end of the time step is found by iteration, starting
  # from the latest estimate. Water content and vapour follow from it
  psiw <- soilwatermod$psiw[idx]
  theta <- soilwatermod$theta[idx]
  vapor <- soilwatermod$vapor[idx]
  # Exponent with which liquid conductivity falls as water potential falls
  nPsi <- 2 + 3 / b
  sh <- dz * (1 - nPsi)

  iter <- 1
  haveStep <- FALSE
  stepFraction <- 1
  while (iter < maxNrIterations) {
    # ---- Water leaving the soil: evaporation less rain, and root uptake ----
    Evapmmhr <- evaporation_flux(theta[1], soilwatermod$Tc[1], climdata$Tair, climdata$relhum,
                                 climdata$rHa, psie[1], b[1], Smax[1], dT)
    surfaceFlux <- (Evapmmhr - climdata$precip) / dT
    # A wetter surface evaporates more. How much more, per unit of water potential
    hp <- max(1e-4 * abs(psiw[1]), 1)
    Evprobe <- evaporation_flux(thetaFromPsi(psiw[1] + hp, psie[1], b[1], Smax[1]), soilwatermod$Tc[1],
                                climdata$Tair, climdata$relhum, climdata$rHa, psie[1], b[1], Smax[1], dT)
    dEvap <- max((Evprobe - Evapmmhr) / (hp * dT), 0)
    STr <- .transpiration_distribute(soilp, rootfrac, climdata$Et, dT, psiw, pTAW)
    STr[frozen] <- 0

    # ---- Hydraulic properties of each node ----
    # Conductivity for liquid water, which does not flow where frozen, and for vapour
    kh <- hydraulicConductivityFromTheta(theta, b, Smax, Ksat)
    kh[frozen] <- 0
    kv <- vaporConductivityFromPsiTheta(psiw, theta, Tk, Smax)
    # Water and vapour that each node gains per unit rise in its water potential
    Cw <- .dTheta_dPsi(psiw, psie, b, Smax)
    Cv <- .dvapor_dPsi(psiw, psie, theta, b, Smax, Tk)
    Ca <- vol * (rho * Cw + Cv) / dT
    # Change in each conductivity per unit rise in water potential
    dkh <- -nPsi * kh / psiw
    airFilled <- Smax - theta
    dkv <- kv * MwR / Tk - ifelse(airFilled > 0, kv / airFilled * Cw, 0)

    # ---- Lower boundary ----
    # Free drainage: the soil beneath is as wet as the deepest node, so water
    # leaves under gravity alone. Otherwise the soil beneath is saturated
    if (soilp$FreeDrain) {
      psiB <- psiw[n]
      khB <- kh[n]
      kvB <- kv[n]
    } else {
      psiB <- psie[n]
      khB <- Ksat[n]
      kvB <- 0
    }
    # Values at the node beneath each node
    psib <- c(psiw[-1], psiB)
    khb <- c(kh[-1], khB)
    kvb <- c(kv[-1], kvB)
    dkhb <- c(dkh[-1], 0)
    dkvb <- c(dkv[-1], 0)

    # ---- Upward flow of water between each node and the one beneath (kg m^-2 s^-1) ----
    # Liquid water moves towards lower water potential. Because conductivity
    # changes steeply between two nodes, the flow is found by integrating it
    # between them. Gravity adds a downward flow
    liquid <- (psib * khb - psiw * kh) / sh - g * kh
    # Vapour diffuses towards lower water potential
    dpsiFace <- psib - psiw
    kvFace <- 0.5 * (kv + kvb)
    ff <- liquid + kvFace * dpsiFace / dz
    # Change in that flow per unit rise in the water potential of the upper node
    # (A) and of the lower node (B)
    dFaceA <- -(kh + psiw * dkh) / sh - g * dkh - kvFace / dz + 0.5 * dkv * dpsiFace / dz
    dFaceB <- (khb + psib * dkhb) / sh + kvFace / dz + 0.5 * dkvb * dpsiFace / dz
    dFaceB[n] <- 0
    if (soilp$FreeDrain) dFaceA[n] <- -g * dkh[n]

    # ---- Water balance of each node ----
    # Water arriving from below, less water leaving upwards and taken up by
    # roots, should equal the change in water and vapour stored. What is left
    # over is the imbalance
    ffabove <- c(0, ff[-n]) # flow across the top of each node
    dd <- ffabove + STr - ff + vol * (rho * (theta - oldtheta) + (vapor - oldvapor)) / dT
    dd[1] <- dd[1] + surfaceFlux
    # How the imbalance of each node changes with its own water potential (bb)
    # and with that of the nodes above (aa) and below (cc)
    aa <- c(0, dFaceA[-n])
    bb <- c(0, dFaceB[-n]) - dFaceA + Ca
    bb[1] <- bb[1] + dEvap
    cc <- -dFaceB
    massBalance <- sum(abs(dd))
    # A saturated surface can take in no more water. Its water potential is
    # held and the excess runs off
    surplusRate <- 0
    if (psiw[1] >= psie[1] - 2e-8 && dd[1] < 0) {
      surplusRate <- -dd[1]
      massBalance <- massBalance - abs(dd[1])
      aa[1] <- 0
      bb[1] <- 1
      cc[1] <- 0
      dd[1] <- 0
    }
    if (massBalance <= tolerance) break

    # ---- Step towards balance ----
    if (haveStep && massBalance > originResidual && stepFraction > 1 / 128) {
      # The last step made the imbalance worse, so it is retaken at half its length
      stepFraction <- 0.5 * stepFraction
      psiw <- stepOrigin - stepFraction * stepTaken
    } else {
      # Change in water potential that would bring every node into balance if
      # the imbalance changed in proportion (Newton's method)
      stepOrigin <- psiw
      stepTaken <- ThomasBoundaryCondition(aa, bb, cc, dd, numeric(n))$x
      psiw <- psiw - stepTaken
      originResidual <- massBalance
      stepFraction <- 1
      haveStep <- TRUE
    }
    # Water potential stays between oven dryness and air entry
    psiw <- pmax(pmin(psiw, psie - 1e-8), psi_min)
    theta <- thetaFromPsi(psiw, psie, b, Smax)
    vapor <- vaporFromPsi(psiw, theta, Tk, Smax)
    iter <- iter + 1
  }
  # State of the soil at the end of the time step, with the lower boundary
  soilwatermod$psiw <- c(psiw, psiB)
  soilwatermod$theta <- c(theta, if (soilp$FreeDrain) theta[n] else Smax[n])
  soilwatermod$vapor[idx] <- vapor
  soilwatermod$k <- c(kh + kv, khB + kvB)
  soilwatermod$Evapmmhr <- Evapmmhr
  soilwatermod$surplus <- surplusRate * dT
  soilwatermod$drainage <- -ff[n] * dT
  soilwatermod$success <- (massBalance <= tolerance)
  soilwatermod$iterations <- iter
  soilwatermod$error <- massBalance
  return(soilwatermod)
}
