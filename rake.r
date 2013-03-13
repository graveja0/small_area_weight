##################################################################
#  R function for loglinear raking of a dataset
#   rake() does loglinear weighting to controls
#   rake.trim() does the same subject to constraints on weights relative
#           to mean weight.
#   19 Sep 2008: modified to control largest step (in rake() ONLY
##################################################################

rake = function(xx,ww,xx.control,beta,iters=10,tolerance = 1, intercept=T,
                maxstep=1e5,trace=T,tracefile="",check.rank=F,warn=T)
{  ## fit a weighting model of the form w.new = w*exp(xx%*%beta)
  ## subject to constraint sum(w.new * xx) = xx.control
  # ARGUMENTS
  #	xx = matrix of variables to be controlled by weighting
  #	ww = vector of base weights (assumed to be all 1 if missing)
  #	xx.control = vector of control totals corresponding to columns of xx
  #			(see also intercept argument)
  #	beta = initial guess for coefficients of loglinear weighting model
  #	iters = maximum number of iterations allowed
  #	tolerance = criterion for stopping iteration, maximum difference
  #		between weighted sum and control totals
  #	intercept = logical, should a column of 1's be prepended to xx?
  #		if yes, first element of xx.control is assumed to correspond
  #		to weighted sum of this column
  #	maxstep = upper bound on absolute change in components of beta
  #		on any iteration (either a single number or vector like beta)
  #	trace = logical, should tracing information be printed with cat()
  #	check.rank = logical or numeric, should xx.control be checked for
  #		nonsingularity and redundant columns dropped accordingly?
  #	warn = logical, emit message if columns dropped?
  # VALUE: object with fitted exponential weighting model
  # 	$maxdiff = maximum error of fitted - control weights
  #	$beta = coefficients of weighting model
  #	$w.rake = exponential factor in weights
  #	$w.both = combined weight ww*w.rake
  ### set up and check arguments
  if(intercept) xx=cbind(1,xx)
  nvars=ncol(xx)
  if(missing(ww)) ww = rep(1,nrow(xx)) else {
    ww = c(ww)
    if(length(ww) != nrow(xx)) stop(paste("ww is wrong length,",length(ww),
                                          "instead of",nrow(xx)))
  }
  if(missing(beta)) {
    beta = rep(0,nvars)
    if(intercept) beta[1]=log(xx.control[1]/sum(ww))
  } else if(length(beta) != nvars) stop("starting beta is wrong length")
  if(length(maxstep) == 1 ) maxstep=rep(maxstep,nvars)
  if(length(maxstep) != nvars) stop("maxstep is wrong length")
  if(any(maxstep<=0)) stop("zero or negative values in maxstep")
  if(length(xx.control) != nvars) stop(paste(
    "xx.control is wrong length,",length(xx.control),"instead of",nvars))
  ### singularity check
  if(check.rank) {
    tol=if(is.logical(check.rank)) 1e-6 else check.rank
    qrx=qr(xx)[c("rank","pivot")]
    dropped=setdiff(1:nvars,qrx$pivot[1:qrx$rank])
    if(length(dropped)) {
      xx.control=xx.control[-dropped]
      xx=xx[,-dropped]
      beta=beta[-dropped]
      maxstep=maxstep[-dropped]
      if (warn) cat("raking, dropped=",dropped,", dim(xx)=",dim(xx),",
		 length(xx.control)=",length(xx.control),"\n")
     }
  }
  ### begin main loop
  for(iter in 1:iters) {
    #browser()
    w.rake = c(exp(xx %*% beta))
    w.both = w.rake * ww
    xx.sum = c(w.both %*% xx)
    maxdiff = max(abs(xx.sum-xx.control))
    if(maxdiff<tolerance) break
    if(trace) cat(file=tracefile, "After iteration",
                  iter,"largest difference =",format(maxdiff,digits=4),"\n")
    deriv = t(xx) %*% (w.both *xx)
    step = solve(deriv,xx.sum-xx.control)
    step = ifelse(abs(step)<=maxstep,step,sign(step)*maxstep)
    beta = beta - step
  }
  if(trace) cat(file=tracefile, 
                "Largest difference =",format(maxdiff,digits=4),"FINISHED\n")
  list(maxdiff=maxdiff,beta=beta,w.rake=w.rake,w.both=w.both,xx.sums=xx.sum,
       dropped=if(check.rank) dropped)
}

rake.trim = function(xx,ww,xx.control,beta,iters=10,tolerance = 1, intercept=T,
                     w.max,w.max.try=.95*w.max, w.min=0,w.min.try=1.1*w.min, 
                     iters.trim=10, w.trim=1)
{	## fit a weighting model of the form w.new = w*exp(xx%*%beta)
  ## subject to constraint sum(w.new * xx) = xx.control
  # ARGUMENTS
  #	xx, ww, xx.control, beta, iters, tolerance, intercept:
  #		values passed to rake()
  #	w.max = maximum allowed combined weight as multiple of mean weight
  #	w.max.try = targeted value of maximum wt before re-iteration of raking
  #	w.min = minimum allowed combined weight as multiple of mean weight
  #	w.min.try = targeted value of minimum wt before re-iteration of raking
  #	iters.trim = maximum iterations of trimmming
  #	w.trim = initial estimates of trimming factor applied to base weights
  # VALUE: object with fitted exponential weighting model
  # 	$maxdiff = maximum error of sum(fitted) - control weights
  #	$beta = coefficients of weighting model
  #	$w.rake = exponential factor in weights
  #	$w.trim = trimming factor applied to base weights
  #	$w.both = combined weight ww*w.trim*w.rake
  ### set up and check arguments
  if(intercept) xx=cbind(1,xx)
  nvars=ncol(xx)
  if(missing(ww)) ww = rep(1,nrow(xx)) else {
    ww = c(ww)
    if(length(ww) != nrow(xx)) stop(paste("ww is wrong length,",length(ww),
                                          "instead of",nrow(xx)))
  }
  if(missing(beta)) {
    beta = rep(0,nvars)
    if(intercept) beta[1]=log(xx.control[1]/sum(ww))
  } else if(length(beta) != nvars) stop("starting beta is wrong length")
  if(length(xx.control) != nvars) stop(paste(
    "xx.control is wrong length,",length(xx.control),"instead of",nvars))
  for(iter in 1:iters.trim) {
    rake.out = rake(xx,ww*w.trim,xx.control,beta,iters,tolerance,intercept=F)
    w.both.ratio = rake.out$w.both/mean(rake.out$w.both)
    cat("TRIMMING ITERATION ", iter,
        ": CV(w.trim)=",format(sqrt(var(w.trim))/mean(w.trim),digits=3),
        ", CV(w)=",format(sqrt(var(w.both.ratio)),digits=3),
        ", range(w/mean)=(",format(min(w.both.ratio),digits=3),
        ", ",format(max(w.both.ratio),digits=3),")",
        ",\n       % trimmed (below,above)=(",
        paste(format(100*c(sum(w.trim>1),sum(w.trim<1))/length(ww),digits=3),
              collapse=","),
        ")\n",sep="")
    beta = rake.out$beta
    if(max(w.both.ratio) <= w.max && min(w.both.ratio)>=w.min) break
    w.trim = w.trim*ifelse(w.both.ratio>w.max,w.trim*(w.max.try/w.both.ratio),1)
    if(w.min) w.trim = w.trim*
      ifelse(w.both.ratio<w.min,w.trim*(w.min.try/w.both.ratio),1)
  }
  rake.out$w.trim=w.trim
  rake.out$summary=round(apply(cbind(base=ww,rake=rake.out$w.rake,
                                     trim=rake.out$w.trim,final=rake.out$w.both),2,
                               function(wt) { 
                                 wt=wt/mean(wt)
                                 c(quantile(wt,c(0,.05,.25,.5,.75,.9,.95,.99,1)),CV=sqrt(var(wt))) }),3)
  rake.out$trimmed = c(below=sum(w.trim>1),above=sum(w.trim<1))/length(ww)
  rake.out
}
