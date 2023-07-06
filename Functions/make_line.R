make_line <- function(start_x, start_y, end_x, end_y) {
  sf::st_linestring(matrix(c(start_x, end_x, start_y, end_y), 2, 2))
}