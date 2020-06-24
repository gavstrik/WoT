# Various functions for the analysis

gmean <- function(x, w = NULL, na.rm=TRUE){
  # Calculate (weighted) geometric mean of input
  y = log(x)
  if(is.null(w)){
    z = exp(mean(y,na.rm=na.rm))
  } else{
    z = exp(sum(w*log(x),na.rm = TRUE)/sum(w,na.rm = TRUE))
  }

  return(z)
}




getParameters <- function(model, se.fit=TRUE){
  # Extract GMM parameters

  pnames    = names(getpars(model))               # use this to find parameters
  delta.idx = grepl("pr",pnames)                  # state probabilities
  mu.idx    = grepl("Intercept",pnames)           # intercept (group bias)
  sig.idx   = grepl("sd",pnames)                  # standard deviations
  B.idx     = (!sig.idx & !mu.idx & !delta.idx)   # all the rest (might be more than one type of parameter)
  nstates   = sum(delta.idx)                      # number of states in model
  nvar      = sum(B.idx)/nstates                  # number of independent/regressor variables in model

  parnames  = character(length(pnames))           # parameter names, with identified states
  parnames[mu.idx] = paste0("mu",1:nstates)
  parnames[sig.idx] = paste0("sig",1:nstates)
  parnames[delta.idx] = paste0("pr",1:nstates)
  parnames[B.idx] = as.vector(outer(paste0("B",1:nvar),1:nstates,paste0))

  tmp       = data.frame(par=depmixS4::getpars(model),"ci.lo"=NA,"ci.hi"=NA)
  covmat = TRUE
  tryCatch({ vcov(model) },
           warning = function(w) {
             covmat <<- FALSE
             # warning("Hessian is singular")
           },
           error   = function(e) {
             covmat <<- FALSE
             # warning("Hessian is singular" )
           })

  if(se.fit & covmat){
    tmp = as.data.frame(suppressWarnings(confint(model)[,-2]))
    colnames(tmp) = c("par","ci.lo","ci.hi")
  }
  tmp$approx=FALSE

  if( any( is.na(tmp$ci.lo) | is.na(tmp$ci.hi) ) ){
    idx = which( is.na(tmp$ci.lo) | is.na(tmp$ci.hi) )
    tmp$approx[idx] = TRUE
  }

  rownames(tmp) = parnames

  delta = tmp[delta.idx,]
  mu    = tmp[mu.idx,]
  sig   = tmp[sig.idx,]
  B     = tmp[B.idx,]

  # if(se.fit & covmat){
  rownames(delta) = parnames[delta.idx]
  rownames(mu) = parnames[mu.idx]
  rownames(sig) = parnames[sig.idx]
  rownames(B) = parnames[B.idx]
  # } else{
  #   names(delta) = parnames[delta.idx]
  #   names(mu) = parnames[mu.idx]
  #   names(sig) = parnames[sig.idx]
  #   names(B) = parnames[B.idx]
  # }

  out = list(delta=delta,mu=mu,B=B, sig=sig, pars=tmp)
  return(out)
}


get_res <- function(M,dat){
  # Make residuals from model, based on weighted states and fitted pars

  nB = length(M$B$par)/M$states
  if(nB==1){
    B  = M$B$par
    Y = dat$y
    X = dat$x1
    W = as.matrix(M$sweight)
    A = M$mu$par
    res = Y-W%*%B*X-W%*%A
  }
  S   = M$sig$par
  res = res/sqrt(W^2%*%S^2)

  return(res)
}


fit_gmm <- function(formula, in.dat, opt=NULL,trystates=2:3, maxiter=10){
  # Fit a Gaussian mixture model
  
  Mstat = list()
  M.fit = list()
  # h = unique(in.dat$method)
  # n = unique(in.dat$v)
  # d = unique(in.dat$d)
  
  for(i in 1:length(trystates)){
    
    M = mix(formula, data = in.dat, nstates = trystates[i], family = gaussian())
    
    iter=0
    ok = TRUE
    while(iter<maxiter & ok){
      tryCatch( expr    = {
        capture.output({tmp.fit  = fit(M)})
        aic      = AIC(tmp.fit,k=2)
        bic      = AIC(tmp.fit,k=log(nrow(in.dat)))
        tmp.crit = c(aic,bic, logLik(tmp.fit))
      },
      warning = function(w) {
        tmp.crit <<- c(NA,NA,NA)
        # warning("Problem fitting ",trystates[i]," states for type ",h," and ",d," dots with ",n," views" )
      },
      error   = function(e) {
        tmp.crit <<- c(NA,NA,NA)
        # warning("Error fitting ",trystates[i]," states for type ",h," and ",d," dots with ",n," views" )
      })
      if(all(!is.na(tmp.crit))){
        ok=FALSE
      }
      iter=iter+1
    }
    
    M.fit[[i]] = tmp.fit
    Mstat$crit = rbind(Mstat$crit,c(trystates[i],tmp.crit))
    
  }
  
  
  colnames(Mstat$crit) = c("nstates","aic","bic","logLik")
  Mstat$crit = as.data.frame(Mstat$crit)
  
  tmp.aic = which(Mstat$crit$aic == min(Mstat$crit$aic,na.rm = TRUE))
  tmp.bic = which(Mstat$crit$bic == min(Mstat$crit$bic,na.rm = TRUE))
  if(all(is.na(Mstat$crit$aic))){
    tmp.opt = c(NA,NA)
  } else{
    tmp.opt = as.data.frame(t( Mstat$crit$nstates[c(tmp.aic,tmp.bic)] ))
  }
  names(tmp.opt) = c("aic","bic")
  Mstat$optimal = tmp.opt
  
  if(all(is.na(tmp.opt))){
    M.out = NULL
  } else{
    if(is.null(opt)){
      M.out = list(fit = M.fit, states=NULL, optimal=tmp.opt, crit=Mstat$crit)
    } else if(tolower(opt)=='bic'){
      M.out = list(fit = M.fit[[tmp.bic]], states = tmp.opt$bic, optimal = tmp.opt, crit=Mstat$crit)
    } else if(tolower(opt)=='aic'){
      M.out = list(fit = M.fit[[tmp.aic]], states = tmp.opt$aic, optimal = tmp.opt, crit=Mstat$crit)
    }
  }
  
  return(M.out)
}
