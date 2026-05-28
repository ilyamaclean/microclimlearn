#' @title Calculates solar position
#' @description Calculates the solar zenith and azimuth angles
#' @param tme POSIXlt object of date and time in UTC
#' @param lat latitude (decimal degrees)
#' @param long longitude(decimal degrees, -ve west of Greenwich meridian)
#' @return a list comprising
#' \describe{
#'  \item{zen}{the solar zenith angle (degrees from vertical)}
#'  \item{azi}{the solar azimuth angle (degrees from north)}
#' }
#' @rdname sunposition
#' @export
sunposition <-function(tme, lat, long) {
  day<-tme$mday
  month<-tme$mon + 1
  year <- tme$year + 1900
  hour <- tme$hour + tme$min / 60 + tme$sec / 3600
  # Calculate Astronomical Julian day
  dd <- day + 0.5 # decimal day
  ma <- month + (month < 3) * 12 # adjusted month
  ya <- year + (month < 3) * -1 # adjusted year
  jd <- trunc(365.25 * (ya + 4716)) + trunc(30.6001 * (ma + 1)) + dd - 1524.5
  B <- (2 -trunc(ya / 100) + trunc(trunc(ya / 100) / 4))
  jd <- jd + (jd > 2299160) * B
  # Calculate solar time
  m<- 6.24004077 + 0.01720197 * (jd - 2451545)
  eot <- -7.659 * sin(m) + 9.863 * sin(2 * m + 3.5932)
  st <- hour + (4 * long + eot) / 60
  # Calculate solar zenith
  latr <- lat * pi/180 # radians
  tt <- 0.261799 * (st - 12)
  dec <- (pi * 23.5/180) * cos(2 * pi * ((jd - 159.5) / 365.25))
  coh <- sin(dec) * sin(latr) + cos(dec) * cos(latr) * cos(tt)
  z <- acos(coh) * (180 / pi)
  # Calculate solar azimuth
  sh <- sin(dec) * sin(latr) + cos(dec) * cos(latr) * cos(tt)
  hh <- (atan(sh / sqrt(1 - sh^2)))
  sazi <- cos(dec) * sin(tt) / cos(hh)
  cazi <- (sin(latr) * cos(dec) * cos(tt) - cos(latr) * sin(dec))/
    sqrt((cos(dec) * sin(tt))^2 + (sin(latr) * cos(dec) * cos(tt) -
                                     cos(latr) * sin(dec))^2)
  sqt <- 1 - sazi^2
  sqt[sqt < 0] <- 0
  azi <- 180 + (180 * atan(sazi / sqrt(sqt))) / pi
  azi[cazi < 0 & sazi < 0] <- 180 - azi[cazi < 0 & sazi < 0]
  azi[cazi < 0 & sazi >= 0] <- 540 - azi[cazi < 0 & sazi >= 0]
  solar <- list(zen = z, azi = azi)
  return(solar)
}
#' @title Calculates solar index
#' @description Calculates fraction of direct beam radiation incident on
#' an inclined surface
#' @param slope slope (decimal degrees from horizontal)
#' @param aspect aspect (decimal degrees from north)
#' @param solp a list of solar zenith and azimuths as returned by [sunposition()]
#' @return fraction of direct beam radiation incident on
#' an inclined surface
#' @rdname solarindex
#' @export
solarindex <- function(slope, aspect, solp) {
  i <- with(solp, cos(zen * pi/180) * cos(slope * pi/180)
            + sin(zen * pi/180) * sin(slope * pi/180) * cos((azi - aspect) * pi/180))
  i[i < 0] <- 0
  i[solp$zen > 90] <- 0
  i
}
#' @title Calculates canopy extinction coefficient
#' @description Calculates the canopy extinction coefficient for inclined and flat canopies
#' @param solp a list of solar zenith and azimuths as returned by [sunposition()]
#' @param x ratio of vertical to horizontal projections of leaf foliage
#' @param si the fraction of direct ebam radiation incident on an inclined surface
#' as returned by [solarindex()]
#' @return a list of the following
#' \describe{
#'  \item{k}{canopy extinction coefficient for vegetation above a flat surface}
#'  \item{k0}{canopy extinction coefficient when solar zenith = 0}
#'  \item{kd}{canopy extinction coefficient for vegetation above an inclined surface}
#' }
#' @rdname cank
#' @export
cank <- function(solp, x, si) {
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
  # k over inclined slopes
  kd <- k * cos(zen * pi/180) / si
  kd[kd > 6000] <- 6000
  return(kd)
}
#' @title Scattering directionality factor
#' @description Calculates the canopy scattering directionality factor for use in two-stream models
#' using one of two methods
#' @param lref leaf reflectance (0-1)
#' @param ltra leaf transmittance. Note that `ltra + lref` must be < 1.
#' @param method one of "Sellers" - leaf directionality factor approximated by
#' calculating mean leaf angle or "Pinty" - leaf directionality factor calculated
#' by integrating the inclination density distribution function.
#' @return J - the scattering directionality factor
#' @rdname scattercoef
#' @export
scattercoef <- function(lref, ltra, x, method = "Sellers") {
  if (method == "Pinty") {
    A<-x+1.774*(x+1.182)^(-0.733)
    if (length(x) >1) {
      J<-A
      s <- which(x <= 1)
      J[s] <- ((x[s]^2 * atan(sqrt(1 - x[s]^2) / x[s])) / (A[s] * (1 - x[s]^2)^(3 / 2)))-
        (x[s]^3 / (A[s] * (1 - x[s]^2) * ((1 - x[s]^2) + x[s]^2)))
      y <- 1/x[-s]
      J[-s] <- -(y^2 * sqrt(1 - y^2) * log((sqrt(1 - y^2) - y^2 +1) / (sqrt(1 - y^2) + y^2 - 1))
                 + 2 * y^2 - 2) / (2 * A[-s] * y * (y^2 - 1)^2)
      J[x==1] <- 1 / 3
    } else {
      if (x <1 ) J <- ((x^2 * atan(sqrt(1 - x^2) / x)) / (A * (1 - x^2)^(3 / 2))) -
          (x^3 / (A * (1 - x^2) * ((1 - x^2) + x^2)))
      if (x>1) {
        y<-1/x
        J <- -(y^2 * sqrt(1 - y^2) * log((sqrt(1 - y^2) - y^2 + 1) / (sqrt(1 - y^2) + y^2 - 1))
               + 2 * y^2 - 2) / (2 * A * y * (y^2 - 1)^2)
      }
      if (x == 1) J <- 1 / 3
    }
  } else {
    if (x == 1) {
      J <- 1/3
    } else {
      mla <- 9.65 * (3 + x)^-1.65 # mean leaf angle
      if (mla > pi/2) mla = pi/2
      J <- cos(mla)^2
    }
  }
  return(J)
}
#' @title Calculate canopy albedo
#' @description Calculates canopy albedos for direct and diffuse radiation using
#' a two-stream model.
#' @param pai total one side area of canopy elements per unit ground area
#' @param lref leaf reflectance (0-1)
#' @param ltra leaf transmittance. Note that `ltra + lref` must be < 1.
#' @param x ratio of vertical to horizontal projections of leaf foliage
#' @param gref ground reflectance
#' @param slope slope of ground surface (decimal degrees from horizontal)
#' @param aspect aspect of ground surface (decimal degrees from north)
#' @param solp a list of solar zenith and azimuths as returned by [sunposition()]
#' @param scattermethod one of "Sellers" - leaf directionality factor approximated by
#' calculating mean leaf angle or "Pinty" - leaf directionality factor calculated
#' by integrating the inclination density distribution function.
#' @param leafh optional logical indicating whether `pai` is the total one-side plant
#' area per unit horizontal area (TRUE) or per sloped ground surface area (FALSE)
#' @return a list of the following
#' \describe{
#'  \item{albd}{combined ground and canopy albedo for diffuse radiation}
#'  \item{albc}{canopy only albedo for direct radiation}
#'  \item{albg}{ground only albedo - same as gref}
#'  \item{albb}{combined ground and canopy albedo for direct radiation}
#'  \item{wgtg}{relative weighting of ground contribution to total albedo for direct radiation}
#'  \item{wgti}{relative weighting of ground contribution to total albedo for diffuse radiation}
#'  \item{K}{canopy extinction coefficient for vegetation above a horizontal surface}
#' }
#' @rdname albedo
#' @export
albedo <- function(pai, lref, ltra, x, gref, slope, aspect, solp,
                   scattermethod = "Sellers", leafh = FALSE) {
  if (leafh) {
    pai <- pai / cos(slope * pi/180)
  }
  zen <- solp$zen
  azi <- solp$azi
  # Calculate canopy extinction coefficient
  si <- solarindex(slope, aspect, solp)
  K <- cank(solp, x, si)
  # Calculate two-stream base parameters
  om <- lref + ltra # single scattering albedo of individual canopy elements
  a <- 1 - om # single scattering absorbtance of individual canopy elements on first interaction.
  del <- lref - ltra # scattering asymmetry
  J <- scattercoef(lref, ltra, x, scattermethod)
  gma <- 0.5 * (om + J * del)
  ss <- 0.5 * (om + J * del / K) * K # back-scatter coefficient
  sstr <- om * K - ss  # scatter coefficient
  # Calculate two-stream base parameters
  h <- sqrt(a^2 + 2 * a * gma)
  S1 <- exp(-h * pai)
  S2 <- exp(-K * pai)
  u1 <- a + gma * (1 - 1 / gref)
  D1 <- (a + gma + h) * (u1 - h) * 1 / S1 - (a + gma - h) * (u1 + h) * S1
  # Calculate two-stream diffuse parameters
  p1 <- (gma / (D1 * S1)) * (u1 - h)
  p2 <- (-gma * S1 / D1) * (u1 + h)
  # Calculate direct parameters
  sig <- K^2 + gma^2 - (a + gma)^2
  p5 <- -ss * (a + gma - K) - gma * sstr
  v1 <- ss - (p5 * (a + gma + K)) / sig
  v2 <- ss - gma - (p5 / sig) * (u1 + K)
  p6 <- (1 / D1) * ((v1 / S1) * (u1 - h) - (a + gma - h) * S2 * v2)
  p7 <- (-1 / D1) * ((v1 * S1) * (u1 + h) - (a + gma + h) * S2 * v2)
  # Calculate albedos
  albd <- (p1 + p2) # white sky albedo
  albb <- (p5 / sig + p6 + p7) # black sky albedo
  # Calculate weights
  wgtg <- exp(-K * pai) # ground weighting (direct)
  wgti <- exp(-pai) # ground weighting (diffuse)
  albc <- (albb - wgtg * gref) / (1 - wgtg) # canopy albedo
  # Calculate leaf absorption coefficient
  sc <- K * si
  return(list(albd = albd, albc = albc, albg = gref, albb = albb, wgtg = wgtg, wgti = wgti, sc = sc))
}
#' @title Calculates radiation absorbed by canopy
#' @description Calculates combined longwave and shortwave radiation absorption.
#' @param climdata a data.frame of weather variables - see e.g. example dataset micropoint::climdata
#' @param vegp a list of vegetation parameters - see e.g. example dataset micropoint::vegparams
#' @param groundp a list of ground parameters - see e.g. micropoint::groundparams
#' @param lat latitude (decimal degrees), negative south of the equator
#' @param long longitude (decimal degrees), negative west of the Greenwich meridian
#' @return radiation absorption by canopy (W/m^2)
#' @rdname calcRabs
#' @export
calcRabs <- function(climdata, vegp, groundp, lat, long) {
  # Calculate solar position
  tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  solp <- sunposition(tme, lat, long)
  albs <- albedo(vegp$pai, vegp$lref, vegp$ltra, vegp$x, groundp$gref,
                 groundp$slope, groundp$aspect, solp)
  # Calculate shortwave radiation absorption
  Rb0 <- (climdata$swdown - climdata$difrad) / cos(solp$zen * pi/180)
  Rb0[Rb0 > 1352] <- 1352
  si <- solarindex(groundp$slope, groundp$aspect, solp) # terrain solar coefficient
  Rabs_sw <- (1 - albs$albd) * climdata$difrad +
    (1 - albs$albg) * albs$wgtg * Rb0 * si +
    (1 - albs$albc) * (1 - albs$wgtg) * Rb0 * albs$sc
  # Calculate longwave radiation absorption
  Rabs_lw <- climdata$lwdown * (albs$wgti * groundp$em + (1 - albs$wgti) * vegp$em)
  return(Rabs_sw +  Rabs_lw)
}
#' @title Two-stream radiation model
#' @description Calculates upward and downward fluxes of shortwave radiation within canopies
#' accomodating inclined surfaces and non-random canopy gaps
#' @param vegp a list of vegetation parameters - see e.g. example dataset micropoint::vegparams
#' @param groundp a list of ground parameters - see e.g. micropoint::groundparams
#' @param paii vector of plant area index values per unit ground area ordered
#' sequentially from canopy bottom to top. Note sum(paii) should equal vegp$pai
#' @param swdown a single numeric value of flux density of downward shortwave
#' radiation on a flat horizontal surface
#' @param difrad a single numeric value of flux density of downward diffuse radiation
#' @param lat latitude (decimal degrees), negative south of the equator
#' @param long longitude (decimal degrees), negative west of the Greenwich meridian
#' @param solp a list of the solar zenith and azimith angle
#' @return a list of radiation fluxes in the same units as `swrad` and `difrad`
#' \describe{
#'  \item{Rdirdown}{Flux density of downward direct radiation perpendicular to the solar beam}
#'  \item{Rdifdown}{Flux density of downward diffuse radiation}
#'  \item{Rswup}{Flux density of upward shortwave radiation - assumed entirely diffuse}
#'  \item{Rleafabs}{Flux density of shortwave radiation absorbed by leaves - averaged over both sides}
#'  \item{Rsun}{Flux density of shortwave radiation absorbed by sunlit leaves}
#'  \item{Rshade}{Flux density of shortwave radiation absorbed by shaded leaves}
#'  \item{sunlitfrac}{Fraction of leaves that are sunlit}
#'  \item{Rgroundabs}{Flux density of shortwave radiation absorbed by ground surface}
#' }
#' @rdname twostream
#' @export
twostream <- function(vegp, groundp, paii, swdown, difrad, lat, long, solp, PAR = FALSE) {
  if (swdown > 0) {
    if (sum(paii) != vegp$pai) {
      vegp$pai <- sum(paii)
      warning("vegp$pai != sum(paii). Adjusting total pai")
    }
    # Cumulative leaf areas
    pia <- rev(cumsum(rev(paii)))
    # Canopy extinction coefficient
    si <- solarindex(groundp$slope, groundp$aspect, solp)
    cosz <- cos(solp$zen * pi /180)
    kd <- cank(solp, vegp$x, si)
    k <-  cank(solp, vegp$x, cosz)
    # Calculate two-stream base parames
    if (PAR) {
      om <- vegp$lrefp + vegp$ltrap
      del <- vegp$lrefp - vegp$ltrap # scattering asymmetry
    } else {
      om <- vegp$lref + vegp$ltra # single scattering albedo of individual canopy elements
      del <- vegp$lref - vegp$ltra # scattering asymmetry
    }
    a <- 1 - om # single scattering absorbtance of individual canopy elements on first interaction.
    J <- with(vegp,scattercoef(lref, ltra, x))
    gma <- 0.5 * (om + J * del)
    ss <- 0.5 * (om + J * del / kd) * kd # back-scatter coefficient
    sstr <- om * kd - ss  # scatter coefficient
    # Calculate two-stream base parameters
    h <- sqrt(a^2 + 2 * a * gma)
    S1 <- exp(-h * vegp$pai)
    S2 <- exp(-kd * vegp$pai)
    u1 <- a + gma * (1 - 1 / groundp$gref)
    u2 <- a + gma * (1 - groundp$gref)
    D1 <- (a + gma + h) * (u1 - h) * 1 / S1 - (a + gma - h) * (u1 + h) * S1
    D2 <- (u2 + h) * 1 / S1 - (u2 - h) * S1
    # Calculate two-stream diffuse parameters
    p1 <- (gma / (D1 * S1)) * (u1 - h)
    p2 <- (-gma * S1 / D1) * (u1 + h)
    p3 <- (1 / (D2 * S1)) * (u2 + h)
    p4 <- (-S1 / D2) * (u2 - h)
    # Calculate direct parameters
    sig <- kd^2 + gma^2 - (a + gma)^2
    p5 <- -ss * (a + gma - kd) - gma * sstr
    v1 <- ss - (p5 * (a + gma + kd)) / sig
    v2 <- ss - gma - (p5 / sig) * (u1 + kd)
    p6 <- (1 / D1) * ((v1 / S1) * (u1 - h) - (a + gma - h) * S2 * v2)
    p7 <- (-1 / D1) * ((v1 * S1) * (u1 + h) - (a + gma + h) * S2 * v2)
    p8 <- sstr * (a + gma + kd) - gma * ss
    v3 <- (sstr + gma * groundp$gref - (p8 / -sig) * (u2 - kd)) * S2
    p9 <- (-1 / D2) * ((p8 / (-sig * S1)) * (u2 + h) + v3)
    p10 <- (1 / D2) * (((p8 * S1) / -sig) * (u2 - h) + v3)
    # Calculate normalised fluxes
    Rdu <- p1 * exp(-h * pia) + p2 * exp(h * pia)  # upward diffuse
    Rdd <- p3 * exp(-h * pia) + p4 * exp(h * pia)  # downwrad diffuse
    Rbu <- (p5 / sig) * exp(-kd * pia) + p6 * exp(-h * pia) + p7 * exp(h * pia) # Contribution of direct to upward diffuse
    Rbd <- (p8 / -sig) * exp(-kd * pia) + p9 * exp(-h * pia) + p10 * exp(h * pia) # Contribution of direct to downward diffuse
    Rbb <- exp(-kd * pia)
    # Calculate radiation streams
    Rbeam <- (swdown - difrad) / cos (solp$zen * pi/180)
    Rbeam[Rbeam > 1352] <- 1352
    Rdirdown <- Rbb * Rbeam
    Rdifdown <- Rdd * difrad + Rbd * Rbeam
    Rswup <- Rdu * difrad + Rbu * Rbeam
    # Calculate radiation absorbed by leaf
    RswL_abs <- a * 0.5 * (k * cosz * Rdirdown + Rdifdown + Rswup)  # average absorbed
    RswL_sun <- a * 0.5 * (k * cosz * Rbeam + Rdifdown + Rswup)
    RswL_shade <- a * 0.5 * (Rdifdown + Rswup)
    # Ground absorbed
    Rddg <- p3 * exp(-h * vegp$pai) + p4 * exp(h * vegp$pai)  # downwrad diffuse at ground
    Rbdg <- (p8 / -sig) * exp(-kd * vegp$pai) + p9 * exp(-h * vegp$pai) + p10 * exp(h * vegp$pai) # Contribution of direct to downward diffuse
    Rdifdowng <- Rddg * difrad + Rbdg * Rbeam
    Rdirdowng <- Rbeam * exp(-kd * vegp$pai)
    RswG_abs <- (1 - groundp$gref) * (si * Rdirdowng + Rdifdowng)
  } else {
    n <- length(pia)
    RswL_abs <- RswL_sun <- RswL_shade <- Rbb <- Rdirdown <- Rdifdown <- Rswup <- rep(0, n)
    RswG_abs <- 0
  }
  return(list(Rdirdown =  Rdirdown,  Rdifdown =  Rdifdown,  Rswup = Rswup, Rleafabs = RswL_abs,
              Rsun = RswL_sun, Rshade = RswL_shade, Rgroundabs = RswG_abs, sunlitfrac = Rbb))
}
#' @title Below-canopy longwave radiation model
#' @description Calculates upward and downward fluxes of longwave radiation within
#' canopies along with longwave radiation absorbed by leaves
#' @param  paii vector of plant area index values per unit ground area ordered
#' sequentially from canopy bottom to top.
#' @param lwdown flux density of downward longwave radiation form sky (W/m^2)
#' @param tleaf vector of leaf temperatures (deg C)
#' @param tground ground surface temperature (deg C)
#' @param vegem emissivity of leaves
#' @param groundem emissivity of ground surface
#' @return a list of radiation fluxes
#' \describe{
#'  \item{Rlwdown}{Flux density of downward longwave radiation (W/m^2)}
#'  \item{Rlwdown}{Flux density of upward longwave radiation (W/m^2)}
#'  \item{RlwLabs}{Flux density of longwave radiation absorbed by leaves - averaged over both sides}
#' }
#' @rdname longwavebelow
#' @export
longwavebelow <- function(paii, lwdown, tleaf, tground, vegem = 0.97, groundem = 0.97) {
  # Calculate pai above and below
  pait <- sum(paii)
  paia <- rev(cumsum(rev(paii)))
  paib <- pait - paia
  # For each element i, calculate weight to every other element
  # weight is paii of j * transmission between i & j
  pij <- abs(outer(paia, paia, "-"))
  wgts <- exp(-pij) * matrix(paii, nrow = length(paii), ncol = length(paii), byrow = TRUE)
  diag(wgts) <- 0
  # Calculate transmission from ground and canopy
  trh <- exp(-paia)
  trg <- exp(-paib)
  # Calculate total canopy weight
  wgt <- rowSums(wgts)
  # Check how close result is to two
  wsum <- wgt + trh + trg
  # Calculate and apply adjustment factor
  mu <- 1 + (2 - 2 * (wsum / 2)) / wgt
  wgts <- wgts * mu # mu recycled
  # Calculate weight of canopy elements for ground absorbed
  wgtg <- trg * paii
  gwsum <- sum(wgtg)
  trsky = exp(-pait)
  mug <- 1 + (1 - (gwsum + trsky)) / gwsum
  wgtg <- wgtg * mug
  # Calculate longwave fluxes
  groundLW <- groundem * 5.67e-8 * (tground + 273.15)^4
  leafLW <- vegem * 5.67e-8 * (tleaf + 273.15)^4
  lwsky <- lwdown * trh  # longwave radiation downward from sky
  lwgro <- groundLW * trg # longwave radiation upward from ground
  # Downward component (j >= i)
  upper_mask <- upper.tri(wgts, diag = TRUE)
  lwcand <- rowSums(wgts * leafLW[col(wgts)] * upper_mask)
  # Upward component (j <= i)
  lower_mask <- lower.tri(wgts, diag = TRUE)
  lwcanu <- rowSums(wgts * leafLW[col(wgts)] * lower_mask)
  # Longwave radiation fluxes
  Rlwdown <- lwsky + lwcand
  Rlwup <-   lwgro + lwcanu
  RlwLabs <- 0.5 * vegem * (Rlwdown + Rlwup)
  return(list(Rlwdown = Rlwdown,
              Rlwup =   Rlwup,
              RlwLabs = RlwLabs))
}

#' Clear-sky solar radiation
#'
#' Calculates clear-sky direct solar radiation from time, location, air
#' temperature, relative humidity and pressure.
#'
#' @param climdata Data frame containing at least \code{obs_time},
#'   \code{temp}, \code{relhum} and \code{pres}.
#' @param lat Latitude (decimal degrees)
#' @param long Longitude (decimal degrees)
#'
#' @return Numeric vector of clear-sky radiation (W m^-2)
#' @export
clearskyrad <- function(climdata, lat, long) {

  tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  temp <- climdata$temp
  relhum <- climdata$relhum
  pres <- climdata$pres
  solp <- sunposition(tme, lat = lat, long = long)
  zenr <- solp$zen * pi / 180
  csr <- rep(0, length(zenr))
  sel <- zenr <= pi / 2
  if (any(sel)) {
    m <- 35 * cos(zenr[sel]) / sqrt(1224 * cos(zenr[sel])^2 + 1)
    TrTpg <- 1.021 - 0.084 * sqrt(m * 0.00949 * pres[sel] + 0.051)
    xx <- log(relhum[sel] / 100) + (17.27 * temp[sel]) / (237.3 + temp[sel])
    Td <- (237.3 * xx) / (17.27 - xx)
    u <- exp(0.1133 - log(3.78) + 0.0393 * Td)
    Tw <- 1 - 0.077 * (u * m)^0.3
    Ta <- 0.935 * m
    od <- TrTpg * Tw * Ta
    csr[sel] <- 1352.778 * cos(zenr[sel]) * od
  }
  csr
}


#' Diffuse proportion of shortwave radiation
#'
#' Calculates the diffuse fraction of incoming shortwave radiation from total
#' downward shortwave flux and solar geometry.
#'
#' @param climdata Data frame containing \code{obs_time} and \code{swdown}.
#' @param lat Latitude (decimal degrees)
#' @param long Longitude (decimal degrees)
#'
#' @return Numeric vector giving the diffuse proportion of shortwave radiation
#' @export
difprop <- function(climdata, lat, long) {

  tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  solp <- sunposition(tme, lat = lat, long = long)
  zenr <- solp$zen * pi / 180
  zd <- solp$zen
  n <- nrow(climdata)
  dp <- rep(1, n)
  sel <- zenr < pi / 2
  if (any(sel)) {
    k1 <- 0.83 - 0.56 * exp(-0.06 * (90 - zd[sel]))
    si <- cos(zenr[sel])
    k <- climdata$swdown[sel] / (1352 * si)
    k <- pmin(pmax(k, 0), k1)
    rho <- k / k1
    sigma3 <- ifelse(
      rho > 1.04,
      0.12 + 0.65 * (rho - 1.04),
      0.021 + 0.397 * rho - 0.231 * rho^2 -
        0.13 * exp(-((rho - 0.931) / 0.134)^2 * 0.834)
    )
    k2 <- 0.95 * k1
    d1 <- rep(1, sum(sel))
    s2 <- zd[sel] < 88.6
    d1[s2] <- 0.07 + 0.046 * zd[sel][s2] / (93 - zd[sel][s2])
    K <- 0.5 * (1 + sin(pi * (k - 0.22) / (k1 - 0.22) - pi / 2))
    d2 <- 1 - ((1 - d1) * (0.11 * sqrt(K) + 0.15 * K + 0.74 * K^2))
    d3 <- (d2 * k2) * (1 - k) / (k * (1 - k2))
    alpha <- (1 / cos(zenr[sel]))^0.6
    kbmax <- 0.81^alpha
    kmax <- (kbmax + d2 * k2 / (1 - k2)) / (1 + d2 * k2 / (1 - k2))
    dmax <- (d2 * k2) * (1 - kmax) / (kmax * (1 - k2))
    dpp <- 1 - kmax * (1 - dmax) / k
    dpp[k <= kmax] <- d3[k <= kmax]
    dpp[k <= k2] <- d2[k <= k2]
    dpp[k <= 0.22] <- 1
    kX <- 0.56 - 0.32 * exp(-0.06 * (90 - zd[sel]))
    kL <- (k - 0.14) / (kX - 0.14)
    kR <- (k - kX) / 0.71
    delta <- rep(0, length(k))
    s3 <- k >= 0.14 & k < kX
    delta[s3] <- -3 * kL[s3]^2 * (1 - kL[s3]) * sigma3[s3]^1.3
    s4 <- k >= kX & k < (kX + 0.71)
    delta[s4] <- 3 * kR[s4] * (1 - kR[s4])^2 * sigma3[s4]^0.6
    dpp[sigma3 > 0.01] <- dpp[sigma3 > 0.01] + delta[sigma3 > 0.01]
    dp[sel] <- dpp
  }
  dp
}
