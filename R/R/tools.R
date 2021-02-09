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
    tmp = as.data.frame(suppressWarnings(depmixS4::confint(model)[,-2]))
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


# Approx. SE when not available (NOT USED)
.approx_se<- function(M){
  # Not yet functional...

  pidx = unique(which(M$weighted$pars$approx))
  for(i in pidx){
    s = as.numeric(rev(unlist(strsplit(rownames(M$weighted$pars[pidx[i],]),split = "")))[1])
    w = M$weighted$sweight[,s]
    x = M$data$guess
    y = M$data$wgmean
    tmp = lm(y ~ x, weights = w)
    tmp.se = confint(tmp)
    mu.idx = grepl("Intercept",rownames(tmp.se))
    B.idx  = !mu.idx
    M$mu[sidx[i],2:3] = tmp.se[mu.idx,]
    M.tmp$B[sidx[i],2:3] = tmp.se[B.idx,]
  }

  idx = unique(which(M.tmp$pars$approx)) # parameters with NA std.errors
  for(i in 1:length(sidx)){
    s = as.numeric(rev(unlist(strsplit(rownames(M.tmp$pars[idx[i],]),split = "")))[1])
    w = M.tmp$sweight[,s]
    x = M.tmp$data$guess
    y = M.tmp$data$gmean.weighted
    tmp = lm(y ~ x, weights = w)
    tmp.se = confint(tmp)
    mu.idx = grepl("Intercept",rownames(tmp.se))
    B.idx  = !mu.idx

    M.tmp$mu[sidx[i],2:3] = tmp.se[mu.idx,]
    M.tmp$B[sidx[i],2:3] = tmp.se[B.idx,]

    tmp.names = rownames(M.tmp$pars)
    if(log(sidx[i],10)<1){
      tmp.idx = substr(tmp.names,nchar(tmp.names),nchar(tmp.names))==sidx[i]
    } else{
      tmp.idx = substr(tmp.names,nchar(tmp.names)-1,nchar(tmp.names))==sidx[i]
    }
    M.tmp$pars[grepl("mu",tmp.names) & tmp.idx,2:3] = tmp.se[mu.idx]
    M.tmp$pars[grepl("B",tmp.names) & tmp.idx,2:3] = tmp.se[B.idx]
  }

}
