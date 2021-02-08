fit_gmm <- function(formula, in.dat, opt=NULL,trystates=2:3, maxiter=10){
  # Fit a Gaussian mixture model
  # @params:
  #   formula     standard R notation as for lm type functions
  #   in.dat      dataset
  #   opt         if NULL, all models for given states are returned, if 'aic' or 'bic', model with lowest score is returned
  #   trystates   number of states to fit, if more than 1 input as vector
  #   maxiter     if a given attempt at fitting doesnt converge, try this number of times to redo the procedure

  Mstat = list()
  M.fit = list()



  # Fitting models with number of states given in 'trystates' and formula as input.
  for(i in 1:length(trystates)){

    M = mix(formula, data = in.dat, nstates = trystates[i], family = gaussian())

    iter=0
    ok = FALSE
    while(iter<maxiter & !ok){
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
      if(all(!is.na(tmp.crit))){ # if algorithm converged, exit and continue
        ok=TRUE
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

  # If 'opt' is given return only the model with lowest score (either AIC or BIC)
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
