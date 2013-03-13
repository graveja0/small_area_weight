D.s.fast = function(dt,xx,ww) {
  D = matrix(0,nrow=length(xx),ncol=length(xx))
  for (r in 1:length(xx)) {
    for (c in 1:length(xx)) {
      D[r,c] <- sum(dt[[ww]]*dt[[xx[r]]]*dt[[xx[c]]])
    }
  }
  return(D)
}

saw = function(
  data =  findat,
  area = "area.id",
  unit = "id",
  ww = "ww.rk",
  xx = c("cons","x"),
  xx.control = NA,
  intercept = FALSE,
#   restrict = NA,
#   restrict.prefix = NA,
  startval = "random",
  maxstep = 1e7,
  maxiter = 3,
  tolerance = 1,
  trace = TRUE,
  restriction.matrix = NA,
  restriction.prefix = "gamma",
  seed,
  breakpoint=5
  ) {

  require(data.table)
  
  
  # Store as a Data Table if Not Already
  if (!is.factor(data[,area])) data[,area] = factor(data[,area])
  if (!is.data.table(data)) data = data.table(data)
  # Extract out various elements
  I = length(xx)
  X = data[,xx,with=FALSE]
  W = data[,ww,with=FALSE]
  A = data[,area,with=FALSE]
  ID = data[,unit,with=FALSE]
  setnames(ID,c("id"))
  # Temporarily assign generic names.
  setnames(X,1:length(xx),paste("x",1:length(xx),sep=""))
  xx.orig = xx
  xx = paste("x",1:length(xx),sep="")
  setnames(W,1,"w")
  setnames(A,1,"area")
  # Cr  
  A.names = names(table(A[,area]))
  A.names.num = as.numeric(A.names)
  S = length(A.names)
  N = nrow(data)
 
  if (trace) {
    cat("\nReweighting Statistics:\n")
    cat(paste("  Number of Small Areas:",S,"\n",sep=" "))
    cat(paste("  Covariates:",length(xx),"\n",sep=" "))
  }

  # Create an Unrestricted Restriction Matrix if One Was Not Specified
  if (missing(restriction.matrix)) {
    if (trace) cat("  No restriction matrix found.  Borrowing will be unrestricted.\n")
    R = data.frame(cbind(area=A.names,matrix(1,nrow=length(A.names),ncol=length(A.names))))
    names(R) = c("area",paste("gamma",A.names,sep=""))
    R.prefix = "gamma"
    R = data.table(R)
    setkey(R,"area")
    restriction.prefix = "gamma"
  } 
  if (!missing(restriction.matrix)) {
    if (!is.factor(restriction.matrix[,area])) restriction.matrix[,area] = factor(restriction.matrix[,area])
    if (trace) cat("  Restriction Matrix Found\n")
    if (is.data.table(restriction.matrix)==FALSE) {
      R = data.table(restriction.matrix) 
    } else R = restriction.matrix
    setnames(R,c("area",paste("gamma",A.names,sep="")))
    setkey(R,"area")
    R.prefix = "gamma"
    R = data.table(R)
    setkey(R,"area")
    restriction.prefix="gamma"
  } 
  #Generate a weighted version of X
  Xw = X
  for (c in 1:I) Xw[,c] = W[,w]*X[[xx[c]]]
  Xw = cbind(A,Xw)
  # Check step shortening inputs
  if(length(maxstep) == 1 ) maxstep=rep(maxstep,I)
  if(length(maxstep) != I) stop("maxstep is wrong length")
  if(any(maxstep<=0)) stop("zero or negative values in maxstep")
  # Set Initial Parameters  
  if (missing(seed)) set.seed(12340) 
  if (startval=="random") startval = rnorm(I*S)
  beta <- matrix(startval,nrow=I,ncol=S)
  if (intercept) beta[1,] <- rep(1,S)
  rownames(beta) <- paste("x",1:length(xx),sep="")
  colnames(beta) <- as.character(A.names.num)
  d.s <- matrix(tolerance+1,nrow=S,ncol=I)
  bx <- matrix(nrow=N,ncol=S)
  # Control Totals 
  if (missing(xx.control)) {
    if (trace) cat("  Missing Control Totals So Using Direct Totals \n")
    X.Cx = Xw[,lapply(.SD,sum),by=area] 
    xx.control = X.Cx
  } else { 
      if (trace)   cat("  Control Totals Detected \n")
      X.Cx = data.table(data.frame(xx.control)[,c(area,xx.orig)])
      setnames(X.Cx,c("area",xx))
  }
  if (trace) cat("Estimating Weights and Applying Restrictions\n\n")
  max = c()
  nonconverge= 0
  for (k in 1:maxiter) {
    # Update Betas
    if (k>1) beta <- newbeta
    #cat(paste(head(w.h)))

    for (s in 1:S) bx[,s] = tcrossprod(as.matrix(X),t(as.matrix(beta[,s])))

    
    # delta.h
    # Applying Restrictions
    exp.bx <- data.table(exp(bx))
    setnames(exp.bx,1:dim(exp.bx)[2],paste("exp.bx",A.names,sep=""))
    exp.bx = cbind(exp.bx,A,ID)
    setkey(exp.bx,area)
    levels(R$area) = levels(exp.bx$area)
    
    exp.bx2 = exp.bx[R] #JG EDIT
    setkey(exp.bx2,id)
    MyMatrix = as.matrix(exp.bx2[,which(grepl("exp.bx",colnames(exp.bx2))),with=FALSE])
    MultBy = as.matrix(exp.bx2[,which(grepl(restriction.prefix,colnames(exp.bx2))),with=FALSE])
    MultBy = apply(MultBy,2,function(i)as.numeric(i))
    exp.bx = data.table(t(t(MyMatrix*MultBy)))
    Sum.exp.bx <- apply(data.frame(exp.bx),1,sum)  
    delta.h <- log( W / Sum.exp.bx)
#     cat("\n\n\n")    
#     cat(paste(head(delta.h)))
#     cat("\n\n\n")
    # Poisson Model 
    ln.wh = data.table(bx)+delta.h$w
    w.h = data.table(exp(ln.wh))
    setnames(w.h,1:dim(w.h)[2],paste("w",A.names,sep=""))
    w.h = cbind(w.h,A,ID)
    setkey(w.h,area)
    w.h = w.h[R]
    setkey(w.h,id)
    w.h.1 = as.matrix(w.h[,which(grepl("w",colnames(w.h))),with=FALSE])
    w.h.r = as.matrix(w.h[,which(grepl(restriction.prefix,colnames(w.h))),with=FALSE])
    w.h.r = apply(w.h.r,2,function(i) as.numeric(i))
    w.h =  data.table(t(t(w.h.1*w.h.r)))


    ###############
    # Step 2
    ##############
    # Update Betas
    
    newbeta <- matrix(0,nrow=I,ncol=S)      
    d.s <-  matrix(0,nrow=S,ncol=I)
    T.s <- matrix(0,nrow=S,ncol=I)
    rownames(d.s) <- A.names
    colnames(d.s)  <- xx
    T.s = data.table(cbind(A.names,T.s))
    setnames(T.s,1,"area")
    setkey(T.s,area)
    
    if (trace) cat(paste("Updating Betas (S=",S, ") s=",sep=""))
    maxi = 0
    maxi.area = c(1)
    pb <- txtProgressBar(min = 0, max = S, style = 3)
    for (s in 1:S) {
      #if (trace & s!=1 & s%%10!=0) cat(".")
      #if (trace & (s==1 | s%%10==0)) cat(s)  
      #if (trace & (s!=1 & s %% 100 ==0)) cat("\n")
      setTxtProgressBar(pb, s)
      Xw2 = X*w.h[[s]]
      temp = data.table(w.h[[s]],X)
      D.s = D.s.fast(dt=temp,xx=xx,ww="V1")
      tot = apply(Xw2,2,sum)
      diff <- t(data.frame(X.Cx[s])[,-1]-tot)
      step <- solve(D.s,diff)
      step = ifelse(abs(step)<=maxstep,step,sign(step)*maxstep)
      newbeta[,s] <- beta[,s] + step 
      if (max(maxi,abs(diff))==max(abs(diff))) {
       maxi.area = A.names[s]
       maxi.Diff = diff
       maxi.n = s
       maxi.d = tot
      }
      maxi <- max(maxi,abs(diff))
    }  
    max[k] = maxi
    if (k>5) {
      if (max[k]>max[k-1]) {
        nonconverge=nonconverge+1
      } else {
        nonconverge = 0
      }
      if (nonconverge==breakpoint) {
        cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nNot converging. Maximum difference found in area",maxi.area,"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",sep=" ")
        print(maxi.Diff)
        cat("Control Total:\n")
        print(X.Cx[maxi.n])
        cat("Indirect Total:\n")
        print(maxi.d)
        cat("Betas:\n")
        print(beta[,maxi.n])
        beta[,maxi.n] = rnorm(I)
        #maxstep = 0.1
        return(data.table(data,w.h))
        break
      }
    }
    if (trace) cat("\n\n")
    if (maxi<tolerance)       break
    
    if (trace) {    
      if (maxi>1) maxpaste = comma(maxi) else maxpaste = round(maxi,5)
      if (k==1) cat("k =",k,": Largest Difference = ",maxpaste,"\n")
      if (k>1) cat("k =",k,": Largest Difference = ",maxpaste,"\n")
      cat("\n")
    } 
    
  }
  output <- data.table(data,w.h)
  return(output)
  
  
} # END SAW FUNCTION


