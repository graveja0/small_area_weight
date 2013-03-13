create_interaction <- function(indata,var1,var2,stem1=var1,stem2=var2,dummy=TRUE) {
  lev1 <- levels(factor(indata[,var1])); lev1
  lev2 <- levels(factor(indata[,var2])); lev2
  nvar <- ncol(indata); nvar
  i <- 1
  for (x in 1:length(lev1)) {
    for (y in 1:length(lev2)) {
      if (dummy) {
        indata[,paste(stem1,lev1[x],"_",stem2,lev2[y],sep="")] <- ifelse(indata[,var1]==lev1[x] & indata[,var2]==lev2[y],1,0    )
      }
      if (!dummy) {
        indata[indata[,var1]==lev1[x] & indata[,var2]==lev2[y],paste(stem1,"_",stem2,sep="")] <- i
      }
      i <- i +1
    }
  }
  if (dummy) {
    begin <- nvar+1
    end <- ncol(indata)
    name <- names(indata[,begin:end])
    out <- indata[,begin:end]
  }
  if (!dummy) {
    out <- indata[,paste(stem1,"_",stem2,sep="")]
  }
  return(out)
}
