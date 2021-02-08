#########################################################################################################################
#   Script for the quantile regression analysis of dots.rda used in the paper:
#     "The Effects of Bounded Social Information on Sequential Decision Making" - Engelhardt et al (2021)
#
#   Author:   Jacob Stærk-Østergaard,
#   email:    ostergaard@sund.ku.dk
#
#########################################################################################################################

library(WoT)
# misc::clean_up() # cleaning the workspace

#########################################################################################################################
#     Setup: declare functions and create variables
#########################################################################################################################
library(quantreg)

# Aux functions
get_coef <- function(x, a=.05){
  est = x$coefficients[,1]
  hi  = x$coefficients[,1]+qnorm(1-a/2)*x$coefficients[,2]
  lo  = x$coefficients[,1]-qnorm(1-a/2)*x$coefficients[,2]
  out = data.frame(tau=x$tau, est = est, lo=lo, hi=hi)
  return(out)
}
qr_est <- function(qfit){
  npar = 8
  ntau = length(taus)
  qres = array(dim=c(4,ntau,npar), dimnames = list(c("tau", "est", "lo", "hi")))
  for(i in 1:length(qfit)){
    tmp = get_coef( qfit[[i]] )
    qres[,i,] = t(tmp)
  }
  return(qres)
}
lm_est <- function(qfit){
  tmp = qfit$coefficients[,1:2]
  out = data.frame(est = tmp[,1], lo = tmp[,1]-1.96*tmp[,2], hi = tmp[,1]+1.96*tmp[,2])
  return(out)
}



# Data is loaded automatically, otherwise use:    load("data/dots.rda")

dots$logerr = log(dots$guess/dots$d)  # create logerr variable to be analyzed
dots$logd2  = log(dots$d)-log(55)     # adjust the number dots so that the d=55 group intercept parameter corresponds to the group quantile
dots$v      = as.factor(dots$v)       # interpret the number of previous estimates visible as a factor variable

#########################################################################################################################
#  Fit regression models (standard and quantile)
#########################################################################################################################

# Set quantiles for regression
taus       = seq(.1,.9,.1)
taus       = c(.05,taus,.95)

# Fit models: both standard regression and quantile regression for each type of social information
fits = list(lm = list(historical  = lm(logerr ~ logd2*v, data=dots, subset=(method=="history")),
                      manipulated = lm(logerr ~ logd2*v, data=dots, subset=(method!="history" | v==0))),
            qr = list(historical  = rq(logerr ~ logd2*v, data=dots, subset=(method=="history"), tau=taus),
                      manipulated = rq(logerr ~ logd2*v, data=dots, subset=(method!="history" | v==0), tau=taus)))

# Make summaries of all models
sums = list(lm = list(historical  = summary(fits$lm$historical),
                      manipulated = summary(fits$lm$manipulated)),
            qr = list(historical  = summary(fits$qr$historical, covariance=TRUE),
                      manipulated = summary(fits$qr$manipulated, covariance=TRUE)))


# Get estimates with confidence bounds
est = list(lm = list(historical  = lm_est(sums$lm$historical),
                     manipulated = lm_est(sums$lm$manipulated)),
           qr = list(historical  = qr_est(sums$qr$historical),
                     manipulated = qr_est(sums$qr$manipulated)))

qr = list(est = est, fit=fits, sum = sums)

fn = "data/qr_est.rda"
save(qr, file=fn)

