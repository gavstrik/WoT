#########################################################################################################################
#   Aux script for the quantile regression analysis of dots.rda used in the paper:
#     "The Effects of Bounded Social Information on Sequential Decision Making" - Engelhardt et al (2021)
#
#   Author:   Jacob Stærk-Østergaard,
#   email:    ostergaard@sund.ku.dk
#
#########################################################################################################################

library(WoT)
# misc::clean_up() # cleaning the workspace

load("data/qr_est.rda")
ndots = c(55,148,403,1097)


#########################################################################################################################
# Control group estimates
#########################################################################################################################

# qr$sum$qr$historical
rownames(qr$fit$qr$historical$coefficients)
tmp = qr$est$qr$historical[-1,6,]
colnames(tmp) = rownames(qr$fit$qr$historical$coefficients)

v0 = tmp[,1:2]
B   = outer(v0[,2],log(ndots/55),"*")
v0 = v0[,1]+B
colnames(v0) = paste0("d=",ndots)
round(v0,4)
round(exp(v0)*matrix(ndots,nr=3,nc=4,byrow=TRUE))
round(exp(v0),4)
v0 = list(log_ratio = v0, ratio = exp(v0), dots = round(exp(v0)*matrix(ndots,nr=3,nc=4,byrow=TRUE)))
v0 = lapply(v0,t)


#########################################################################################################################
# v=1,3,9 groups: (refit model at tau=.5 with v=1,3,9 as reference groups)
#########################################################################################################################

dots$logerr = log(dots$guess/dots$d)  # create logerr variable to be analyzed
dots$logd2  = log(dots$d)-log(55)     # adjust the number dots so that the d=55 group intercept parameter corresponds to the group quantile
dots$v      = as.factor(dots$v)       # interpret the number of previous estimates visible as a factor variable

dots$v1     = relevel(dots$v,ref="1")
dots$v3     = relevel(dots$v,ref="3")
dots$v9     = relevel(dots$v,ref="9")

get_coef <- function(x, a=.05){
  est = x$coefficients[,1]
  hi  = x$coefficients[,1]+qnorm(1-a/2)*x$coefficients[,2]
  lo  = x$coefficients[,1]-qnorm(1-a/2)*x$coefficients[,2]
  out = data.frame(tau=x$tau, est = est, lo=lo, hi=hi)
  return(out)
}

# Set quantiles for regression
tau = .5
qr1 = list(historical  = rq(logerr ~ logd2*v1, data=dots, subset=(method=="history"), tau=tau),
           manipulated = rq(logerr ~ logd2*v1, data=dots, subset=(method!="history" | v==0), tau=tau))
qr3 = list(historical  = rq(logerr ~ logd2*v3, data=dots, subset=(method=="history"), tau=tau),
           manipulated = rq(logerr ~ logd2*v3, data=dots, subset=(method!="history" | v==0), tau=tau))
qr9 = list(historical  = rq(logerr ~ logd2*v9, data=dots, subset=(method=="history"), tau=tau),
           manipulated = rq(logerr ~ logd2*v9, data=dots, subset=(method!="history" | v==0), tau=tau))

make_tables <- function(model){
  tmp1 = get_coef(summary(model, covariance=TRUE))
  tmp2 = t(tmp1[1:2,-1])
  tmp3 = outer(tmp2[,2],log(ndots/55),"*")
  tmp2  = tmp2[,1]+tmp3
  colnames(tmp2) = paste0("d=",ndots)
  out = list(log_ratio = tmp2, ratio = exp(tmp2), dots = round(exp(tmp2)*matrix(ndots,nr=3,nc=4,byrow=TRUE)))
  out = lapply(out,t)
  return(out)
}


## Naming: vxy, where x refers to dots and y to historical or manipulated
res = list(v1h = make_tables(qr1$historical),
           v1m = make_tables(qr1$manipulated),
           v3h = make_tables(qr3$historical),
           v3m = make_tables(qr3$manipulated),
           v9h = make_tables(qr9$historical),
           v9m = make_tables(qr9$manipulated))

res = lapply(res, function(x)  matrix(unlist(x), nrow=4))

xtable::xtable(matrix(unlist(v0),nrow=4), digits=4)
lapply(res, xtable::xtable)


