#' French departement
#'
#' Geographical data of the French administrative units "département"
#' This data set excludes overseas departement as well as the two departements of the island of Corsica.
#'
#' @format An \code{\link[sf]{sf}} object, with lattitude and longitude
#' @source The data set comes from \url{https://github.com/gregoiredavid/france-geojson/}, under the "Licence Ouverte / Open Licence v2.0" licence.
#' @examples
#' data(departement)
#' plot(departement["geometry"])
"departement"
#' IRIS of Paris
#'
#' Geographical data of the French statistical units "IRIS" (\emph{grouped islets for statistical information}) forming the city of Paris.
#' Five iris units (forming two parcs: that of Boulogne and Vincennes) of dispropotionate size have been removed from the data set.
#'
#' @format An \link[sf]{sf} object, with lattitude and longitude
#' @source The data set comes from the French National Statistics Institute and
#' the National Geography Institute. It is accessible at
#' \url{https://www.data.gouv.fr/fr/datasets/contours-iris-insee-ign/},
#' under the "Licence Ouverte / Open Licence v2.0" licence.
#' @examples
#' data(paris)
#' plot(paris["geometry"])
"paris"

#' Administrative areas of the Netherlands aroung the city of Utrecht
#'
#' Geographical data of the Dutch administrative units "District".
#' @format
#' \describe{
#'   \item{utrecht_district}{An \code{\link[sf]{sf}} object, with lattitude and longitude. With
#'   7 variables and 650 zones.}
#'   \item{graph_utrecht_district}{The adjacency graph as an \code{\link{igraph}} object.}
#' }
#' @source The data set comes from \url{https://geodata.nationaalgeoregister.nl/}.
#' @name utrecht_district
#' @examples
#' data(utrecht_district); data(graph_utrecht_district)
#' coord <- sf::st_coordinates(sf::st_centroid(utrecht_district))
#' adj_municip <- as(as(igraph::as_adjacency_matrix(graph_utrecht_district, type = "both"),
#'                      "symmetricMatrix"),
#'                   "TsparseMatrix")
#' edge_list <- data.frame(adj_municip@i + 1, adj_municip@j + 1)
#' segment_df <- cbind(coord[edge_list[, 1], ], coord[edge_list[, 2], ])
#' ptmat <- as.matrix(segment_df[, 1:4])[2:nrow(segment_df), ]
#' linesegs <- lapply(split(ptmat, 1:nrow(ptmat)), function(x) {
#'     x <- matrix(x, nrow = 2, byrow = TRUE)
#'     x <- sf::st_linestring(x)})
#' final_sf <- sf::st_sf(sf::st_sfc(linesegs), 'ID' = 1:length(sf::st_sfc(linesegs)))
#' op <- par(mar = rep(0, 4))
#' plot(sf::st_geometry(utrecht_district), lwd = 0.6, border = "grey")
#' plot(sf::st_geometry(final_sf), lwd = 0.5, add = TRUE)
#' plot(sf::st_centroid(utrecht_district), add = TRUE, col = "black", pch = 20,
#'      cex = 0.5)
#' par(op)
#' @keywords datasets
NULL
#' @rdname utrecht_district
"utrecht_district"
#' @rdname utrecht_district
"graph_utrecht_district"
