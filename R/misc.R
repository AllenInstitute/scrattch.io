collapse_along <- function(x,
                           collapse = "/") {

  out <- vector(length = length(x))

   for(i in 1:length(x)) {
    out[i] <- paste(x[1:i], collapse = collapse)
  }

  out
}
