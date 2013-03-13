comma <- function(number,r=0) {
  out <- round(number)
  c <- nchar(out)
  out2 <- ""
  for (i in 1:c) {
    j <- (c+1)-i
    out2 <- paste(out2,ifelse(round((i-1)/3)==((i-1)/3) & j!= c,",",""),substr(out,j,j),sep="")
  }
  c2 <- floor(nchar(out)/3)+c     
  out3 <- ""  
  for (i in 1:(c2)) {
    j <- (c2+1)-i
    out3 <- paste(out3,substr(out2,j,j),sep="")
  }
  
  return(out3)
}
