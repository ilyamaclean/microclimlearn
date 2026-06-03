#' Plant functional type parameter values
#'
#' Vegetation parameter values for main plant function types
#'
#' @format ## `PFTtable`
#' A data frame with 26 rows and 13 columns:
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
#'   \item{Esh}{Parameter values for Evergreen shrubs}
#'   \item{Dsh}{Parameter values for Deciduous shrubs}
#'   \item{multiplier}{Multiplier used to convert form values in table to values used in model}
#'   ...
#' }
#' @source <https://gmd.copernicus.org/articles/11/2857/2018/gmd-11-2857-2018.pdf>
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
#' "weather"
