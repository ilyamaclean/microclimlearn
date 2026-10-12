#' Plant functional type parameter values
#'
#' Vegetation parameter values for main plant function types
#'
#' @format ## `PFTtable`
#' A data frame with 23 rows and 13 columns:
#' \describe{
#'   \item{varname}{Short-hand name of vegetation parameter}
#'   \item{Description}{Long-hand description of vegetation parameter}
#'   \item{units}{units used}
#'   \item{BET.Tr}{Parameter values for Broadleaf evergreen tropical trees}
#'   \item{BET.Te}{Parameter values for Broadleaf evergreen temperate trees}
#'   \item{BDT}{Parameter values for Broadleaf deciduous trees}
#'   \item{NET}{Parameter values for Needle-leaf evergreen trees}
#'   \item{NDT}{Parameter values for Needle-leaf deciduous trees}
#'   \item{C3}{Parameter values for C3 grasses}
#'   \item{C4}{Parameter values for C4 grasses}
#'   \item{ESh}{Parameter values for Evergreen shrubs}
#'   \item{DSh}{Parameter values for Deciduous shrubs}
#'   \item{multiplier}{Multiplier used to convert from values in table to values used in model}
#'   ...
#' }
#' @source Photosynthesis parameters are from Harper et al. (2016), Table 2; plant hydraulic parameters
#'   are from Eller et al. (2020), Table 2; root depths are from Schenk & Jackson (2002), Table 4.
"PFTtable"

#' A data frame of hourly weather
#'
#' A data frame of hourly weather in 2017 at Caerthillean Cove, Lizard, Cornwall (49.96807N, 5.215668W)
#'
#' @format a data frame with the following elements:
#' \describe{
#'  \item{obs_time}{POSIXlt object of dates and times}
#'  \item{temp}{temperature (degrees C)}
#'  \item{relhum}{relative humidity (percentage)}
#'  \item{pres}{atmospheric press (kPa)}
#'  \item{swdown}{Total downward shortwave radiation (W / m^2)}
#'  \item{difrad}{Total downward diffuse radiation (W / m^2)}
#'  \item{lwdown}{Total downward longwave radiation (W / m^2)}
#'  \item{windspeed}{Wind speed (m/s)}
#'  \item{winddir}{Wind direction (decimal degrees)}
#'  \item{precip}{precipitation (mm)}
#' }
"climdata"

#' Plant area index profile of a broadleaved woodland measured by lidar
#'
#' The vertical distribution of plant area in Wytham Woods, Oxfordshire, UK, in
#' full leaf. It was measured by a fixed, upward-looking lidar sensor (sensor 302
#' of the StrucNet network) and is the mean of the daily profiles from 8 to 22
#' June 2022. The source gives plant area per unit volume (m^2/m^3) at height
#' steps of 0.5 m, which is multiplied here by layer depth to give the plant area
#' index of each layer. Layers are ordered from the ground upward, as the
#' functions of this package that take `paii` require. Canopy height is 24.5 m
#' (49 layers of 0.5 m) and total plant area index is 3.38. The lowest two
#' layers lie below the sensor and hold no plant area.
#'
#' @format a data frame with 49 rows and the following elements:
#' \describe{
#'  \item{z}{height above ground of the middle of each layer (m)}
#'  \item{paii}{plant area index of each layer (m^2 of plant area per m^2 of ground)}
#' }
#' @source Calders K, Brede B, Newnham G, et al. (2023) Vertical plant profiles
#'   for Dassenbos (NL, 2014-2018, TLS); Wytham Woods (UK, 2022, LEAF) & Northern
#'   Australia (2021-2022, LEAF). Zenodo. \doi{10.5281/zenodo.7778743}. Licensed
#'   under CC BY 4.0.
#' @references Calders K, Brede B, Newnham G, et al. (2023) StrucNet: a global
#'   network for automated vegetation structure monitoring. Remote Sensing in
#'   Ecology and Conservation 9: 587-598.
"wythampai"
