# old_zones <- list(1, 20, 381, 400)
# while (length(setdiff(unlist(zones), seq(1, p))) > 0) {
#
#   zones <- adjacent_vertices(graph, old_zones) %>%
#     mapply(union, ., old_zones, SIMPLIFY = FALSE)
#   zones
#   old_zones <- zones
# }
# zones