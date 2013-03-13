recode_factors <- function(data,marginal) {
  x.vars <- c()
  for (c in 1:length(marginal)) {
    lev <- levels(factor(data[,marginal[c]])); lev
    # Factors with 2 levels are dummies and don't need to be edited
    if (length(lev)<=2) {
      x.vars <- c(x.vars,marginal[c])
    }
    if (length(lev)>2) {
      for( x in 1:length(lev)) {
        data[,paste(marginal[c],x,sep="")] <- ifelse(data[,marginal[c]]==lev[x],1,0)
        if (x>1) {
          x.vars <- c(x.vars,paste(marginal[[c]],x,sep=""))
        }
      }    
    }
  }
  out <- list(marginal,x.vars,data)
}
