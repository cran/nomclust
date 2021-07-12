as.hclust.nomclust <- function(x, ...) {
    output <- as.hclust(as.agnes(x, ...))
  }
