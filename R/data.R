#' French departement
#'
#' Geographical data of the French administrative units "département"
#' This data set excludes overseas departement as well as the two departements of the island of Corsica.
#'
#' @format An \code{\link[sf]{sf}} object, with lattitude and longitude
#' @source The data set comes from \url{github.com/gregoiredavid/france-geojson}, under the "Licence Ouverte / Open Licence v2.0" licence.
#' @examples
#' \dontrun{
#' ## Using "plot"
#' plot(departement)
#' ## Using "ggplot"
#' library(ggplot)
#' ggplot(departement) + geom_sf()
#' }
"departement"
#' IRIS of Paris
#'
#' Geographical data of the French statistical units "IRIS" (\emph{grouped islets for statistical information}) forming the city of Paris.
#' Five iris units (forming two parcs: that of Boulogne and Vincennes) of dispropotionate size have been removed from the data set.
#'
#' @format An \link[sf]{sf} object, with lattitude and longitude
#' @source The data set comes from the French National Statistics Institute and the National Geography Institute. It is accessible at \url{www.data.gouv.fr/fr/datasets/contours-iris-insee-ign/}, under the "Licence Ouverte / Open Licence v2.0" licence.
#' @examples
#' \dontrun{
#' library(graphseg)
#' ## Using "plot"
#' plot(paris["geometry"])
#' ## Using "ggplot"
#' library(ggplot)
#' ggplot(paris) + geom_sf()
#' }
"paris"

