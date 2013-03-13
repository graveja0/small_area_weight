defactor <- function(f) as.numeric(levels(f)[as.integer(f)])

# Version for character data
defactor.c <- function(f) {
  l <- as.character(levels(f))
  out <- as.numeric(l)[f]
  return(out)
}
