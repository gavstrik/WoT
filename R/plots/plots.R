#########################################################################################################################
#   Script for plots used in the paper:
#     "The Effects of Bounded Social Information on Sequential Decision Making" - Engelhardt et al (2021)
#
#   Author:   Jacob Stærk-Østergaard,
#   email:    ostergaard@sund.ku.dk
#
#   Input a filename (fn) to save plots as pdf
#########################################################################################################################

library(WoT)
# misc::clean_up() # cleaning the workspace

which_model <- function(type,d,v){
  # aux function to map model i to thread type, dots and views
  if(type=='pristine'){
    idx = which(sessions$method=='history' & sessions$d == d & sessions$v == v)
  } else if(type=='manipulated'){
    idx = which(sessions$method=='max' & sessions$d == d & sessions$v == v)
  }

  i = which(names(model) == sessions$session[idx])
  return(i)
}

#########################################################################################################################
#     1. Quantile regression plots
#########################################################################################################################
load("data/qr_est.rda")

plot_qr_est <- function(est, type=NULL, fn = NULL){
  # Plot quantile regression results
  ntau = dim(est)[2]
  npar = dim(est)[3]

  # Plot results
  if(!is.null(fn)) pdf(file=fn, width=10, height=6)

  par(mfrow=c(2,4), mar=c(2,2,1,1), oma=c(2,1.5,6,1))
  layout(matrix(c(1,2,3,6,4,7,5,8), nr=2))
  for(i in 1:npar){
    x  = est[1,,i]
    y  = est[2,,i]
    lo = est[3,,i]
    hi = est[4,,i]

    # plot(0,0, type='n', xlim=c(0,1), ylim=c(range(qres[-1,i,])))
    yLim = ifelse(i %in% c(2,6,7,8), 1, 1)*c(-1,1)
    plot(0,0, type='n', xlim=c(0,1), ylim=yLim, ann=FALSE)
    polygon(c(x,rev(x)), c(lo, rev(hi)), border=NA, col=add_alpha('black',.25))
    abline(h=0, v=0.5, lty=3, col=add_alpha('black',.75))
    lines(x,y)

    if(i==1){
      mtext(text = "Intercept", side=2, line=2)
      mtext(text = "v=0", side=3, line=1)
      # legend("topleft", c("Quantile","Linear"), col=add_alpha(c("black","red"),.75), lty=1, lwd=2, bty='n')
    } else if(i==2){
      mtext(text = "Slope", side=2, line=2)
    } else if(i==3){
      mtext(text = "v=1", side=3, line=1)
    } else if(i==4){
      mtext(text = "v=3", side=3, line=1)
    } else if(i==5){
      mtext(text = "v=9", side=3, line=1)
    }
  }
  mtext("Reference group", side=3, line = 1.7, outer=TRUE, adj = .1)
  mtext("Contrast groups", side=3, line = 1.7, outer=TRUE, adj = .64)
  mtext("Quantiles", side=1, line = 1, outer=TRUE)
  if(!is.null(type)){
    mtext(type, side=3, line = 4, outer=TRUE, cex=1.5)
  }


  if(!is.null(fn)) dev.off()
}
plot_qrlines <- function(qfit, dots, type=NULL, fn = NULL){

  plot_qrline <- function(v,x, yLim = c(0,4), logscale=FALSE){
    taus = as.numeric(substr(colnames(v),6, nchar(colnames(v))))
    ntau = ncol(v)
    plot(0,0, xlim=range(exp(x+log(55))), ylim=yLim, type='n', bty='n')
    for(i in 1:ntau){
      if(taus[i]==.5){
        cols = add_alpha('tomato2',.95)
        Lty  = 1
        Lwd  = 1.5
      } else{
        cols = add_alpha('tomato2',.95)
        Lty  = 2
        Lwd  = 1.5
      }
      if(logscale){
        abline(h=0, lty=3, col=add_alpha('black',.5))
        lines(exp(x+log(55)), exp(v[1,i]+v[2,i]*x), col=cols, lty=Lty, lwd=Lwd)
        text(max(exp(x+log(55))),exp(v[1,i]+v[2,i]*max(x)), labels = taus[i], pos = 4, cex=.75, offset = .1)

      } else{
        abline(h=1, lty=3, col=add_alpha('black',.5))
        lines(exp(x+log(55)), exp(v[1,i]+v[2,i]*x), col=cols, lty=Lty, lwd=Lwd)
        text(max(exp(x+log(55))),exp(v[1,i]+v[2,i]*max(x)), labels = taus[i], pos = 4, cex=.75, offset = .1)

      }
    }
  }


  # Plot quantile regression lines
  v0 = coef(qfit)[1:2,]
  v1 = coef(qfit)[1:2,]+coef(qfit)[c(3,6),]
  v3 = coef(qfit)[1:2,]+coef(qfit)[c(4,7),]
  v9 = coef(qfit)[1:2,]+coef(qfit)[c(5,8),]

  # Plot results
  if(!is.null(fn)) pdf(file=fn, width=10, height=6)
  par(mfrow=c(2,2), mar=c(3,3,1,1), oma=c(1.5,1.5,3.5,0))
  x = seq(0,3.1,.1)

  plot_qrline(v0,x); mtext("v=0", side=3, line=0)
  tmp = subset(dots, method=="history" & v==0)
  tmp$logerr = log(tmp$guess/tmp$d)
  points(tmp$d,exp(tmp$logerr), pch=16, col=add_alpha('black',.25), cex=.5)

  plot_qrline(v1,x); mtext("v=1", side=3, line=0)
  tmp = subset(dots, method=="history" & v==1)
  tmp$logerr = log(tmp$guess/tmp$d)
  points(tmp$d,exp(tmp$logerr), pch=16, col=add_alpha('black',.25), cex=.5)

  plot_qrline(v3,x); mtext("v=3", side=3, line=0)
  tmp = subset(dots, method=="history" & v==3)
  tmp$logerr = log(tmp$guess/tmp$d)
  points(tmp$d,exp(tmp$logerr), pch=16, col=add_alpha('black',.25), cex=.5)

  plot_qrline(v9,x); mtext("v=9", side=3, line=0)
  tmp = subset(dots, method=="history" & v==9)
  tmp$logerr = log(tmp$guess/tmp$d)
  points(tmp$d,exp(tmp$logerr), pch=16, col=add_alpha('black',.25), cex=.5)
  mtext("Number of dots visible", side=1, outer=TRUE, line=0)
  mtext("Relative error", side=2, outer=TRUE, line=0)
  if(!is.null(type)){
    mtext(type, side=3, line = 2, outer=TRUE, cex=1.5)
  }
  if(!is.null(fn)) dev.off()
}

plot_qr_est(qr$est$qr$historical,  type='Pristine',  fn = NULL)
plot_qr_est(qr$est$qr$historical,  type='Pristine',  fn = "plots/historical_qr.pdf")
plot_qr_est(qr$est$qr$manipulated, type='Manipulated', fn = "plots/manipulated_qr.pdf")

plot_qrlines(qr$fit$qr$historical,  dots, type='Pristine',  fn = NULL)
plot_qrlines(qr$fit$qr$historical,  dots, type='Pristine',  fn = "plots/historical_qrlines.pdf")
plot_qrlines(qr$fit$qr$manipulated, dots, type='Manipulated', fn = "plots/manipulated_qrlines.pdf")

#########################################################################################################################
#     2. QQ-plots for GMMs
#########################################################################################################################
load("data/models.rda")
  plot_gmm_res <- function(i, info=TRUE){
    M.tmp = model[[i]]$weighted
    dat   = model[[i]]$data
    res   = get_res(M.tmp,dat)
    cols = add_alpha('black',.5)

    d = unique(dat$d)
    v = unique(dat$v)
    h = unique(dat$method)
    s = model[[i]]$session
    n = nrow(dat)
    # if(i %in% c(1,2,3,5,14,15,16,17,18,23, 24)){
    #   cols = add_alpha('red',.75)
    # }
    qqnorm(res, pch=16, col=cols, bty='n', xlim=c(-4,4), ylim=c(-4,4), main=paste("Session no.",i,"id:",s), cex.main=.9)

    abline(0,1, lty=3)

    if(info){
      if(h == "history"){
        h = "pristine"
      }
      legend("topleft", c(paste("History:  ",h),
                          paste("Dots:      ",d),
                          paste("Views:    ",v),
                          paste("N obs:    ",n)), cex=.25, bty='n')
    }
  }

# Plot all models
par(mfrow=c(6,5), mar=c(2,2,1,1), oma=c(0,0,0,0))
for(i in 1:length(model))
  plot_gmm_res(i)

layout(1)
plot_gmm_res(which_model("pristine",55,9))


#########################################################################################################################
#     3. Weighted beta plots for GMMs
#########################################################################################################################
load("data/models.rda")
plot_beta <- function(M.tmp, probs=NULL,xlim=NULL, noCov=FALSE, ann=FALSE){

  wB   = as.matrix(M.tmp$sweight)%*%M.tmp$B$par
  cols = wbeta_colors(wB)


  # Plot weighted betas
  V = vcov(M.tmp$fit)$vcov
  rownames(V) = rownames(M.tmp$pars)#[!M.tmp$pars$approx]
  colnames(V) = rownames(M.tmp$pars)#[!M.tmp$pars$approx]

  W = as.matrix(M.tmp$sweight)
  idx = which(grepl("B",rownames(V)))

  VB = V[idx,idx]
  B = M.tmp$B$par

  wB = W%*%B
  wV = diag(W%*%VB%*%t(W))

  idx = order(wB)
  y = seq(0,1,length=length(wB))
  x = wB[idx]
  n = length(x)

  if(all(wV>=0) & !noCov){
    wB.hi = wB+1.96*sqrt(wV)
    wB.lo = wB-1.96*sqrt(wV)
    x.hi = wB.hi[idx]
    x.lo = wB.lo[idx]
  } else {
    x.hi = rep(NA,n)
    x.lo = rep(NA,n)
    x.lo = x.hi = x

  }

  if(length(cols)==n){
    cols = cols[idx]
  }
  if(is.null(xlim)){
    plot(c(min(x.lo),max(x.hi)), c(0,1), type='n', bty='n', lwd=2, ann=FALSE, yaxt='n')
  } else {
    plot(c(min(x.lo),max(x.hi)), c(0,1), type='n', bty='n', lwd=2, xlim=xlim, ann=FALSE, yaxt='n')
  }

  axis(side = 2, at = seq(0,1,.25), las=1, cex=1.5)

  idx.hi = order(x.hi)
  idx.lo = order(x.lo)

  segments(x.hi[idx.hi],y,x.lo[idx.lo],y, col=add_alpha(cols,.5),lwd=2)
  segments(x[-n],y[-n],x[-1],y[-1], col=add_alpha(cols,.85), lwd=2)

  # mtext("CDF",side = 2, line=2.5, cex=1.5)
  if(ann){
    mtext(expression(tilde(beta)),side = 1, line=3, cex=.75)
  }

  # plot(wB[idx],y, type='l', bty='n', lwd=2, xlim=c(0,1))
  # polygon(c(wB.hi[idx],rev(wB.lo[idx])), c(y,rev(y)), border=NA, col=add_alpha('grey',.5))
  # points(x,y, pch=16, col=rcols[idx], cex=.5)


  if(!is.null(probs)){
    for(p in probs){
      idx1=max(which(y<= p))
      idx2=min(which(y > p))
      y.tmp = (y[idx1]+y[idx2])/2
      x.tmp = (x[idx1]+x[idx2])/2
      segments(par()$usr[1],y.tmp,x.tmp,y.tmp, lty=3)
      segments(x.tmp,par()$usr[3],x.tmp,y.tmp, lty=3)

      # text(par()$usr[1],y.tmp,round(y.tmp,2),cex=0.75, adj=c(-.2,-.3))
      # text(x.tmp,par()$usr[3],round(x.tmp,2),cex=0.75, adj=c(-.3,-.4))

      # text((par()$usr[1]+x.tmp)/2,y.tmp,round(y.tmp,2),cex=0.75, adj=c(-.2,-.3))
      # text(x.tmp,(par()$usr[3]+y.tmp)/2,round(x.tmp,2),cex=0.75, adj=c(-.3,-.4))

      # text((par()$usr[1]+x.tmp)/2,y.tmp,bquote("p="~.(round(y.tmp,2))),cex=0.75, adj=c(-.2,-.3))
      # text(x.tmp,(par()$usr[3]+y.tmp)/2,bquote(beta~"="~.(round(x.tmp,2))),cex=0.75, adj=c(-.3,-.4))

      # text(x.tmp-.03,y.tmp,bquote("p="~.(round(y.tmp,2))),cex=0.75, adj=c(-.2,-.3))
      # text(x.tmp,y.tmp-.03,bquote(beta%~~%.(round(x.tmp,3))),cex=0.75, adj=c(-.3,-.4))
      text(x.tmp,y.tmp-.03,bquote(tilde(beta)~"="~.(round(x.tmp,3))),cex=1.2, adj=c(-.3,-.4))

    }
  }
  # mtext(paste0("Type: ",unique(M.tmp$data$hist), ", ",unique(M.tmp$data$dots)," dots and ", unique(M.tmp$data$nobs), " views") , side=3, line=1, font=2)

}

# Plot all models
par(mfrow=c(6,5), mar=c(2,2,1,1), oma=c(0,0,0,0))
for(i in 1:length(model))
  plot_beta(model[[i]]$weighted, xlim=c(-.5,1.5))


# Plot one model
layout(1)
plot_gmm_res(which_model("pristine",55,9))


#########################################################################################################################
#     4. Plot threads
#########################################################################################################################
load("data/models.rda")
plot_thread <- function(i,Xlim=c(-2,5), fn=NULL){
  Xlim = c(-2,5)

  M.tmp = model[[i]]$weighted
  dat   = model[[i]]$data
  wB = as.matrix(M.tmp$sweight)%*%M.tmp$B$par

  Bcols = wbeta_colors(wB)

  h = unique(model[[i]]$data$method)
  d = unique(model[[i]]$data$d)
  v = unique(model[[i]]$data$v)
  n = nrow(model[[i]]$data)
  s = model[[i]]$session

  if(!is.null(fn)) pdf(fn, width=12, height=10)
  plot(dat$y, col=add_alpha(Bcols,.8), pch=16,cex.axis=1, bty='n', ylim=Xlim, ann=FALSE)
  lines(dat$x1, lty=3, col=add_alpha('black',.5))
  mtext(expression("log-ratio ("~y[i]~")"), side = 2, line = 1.8, cex=.8, las=0)
  mtext("Observation no.", side = 1, line = 2.25, cex=.8)
  abline(h=0, col=add_alpha('black',.3))

  mtext(paste("Session no.",i,"id:",s), 3, line=0.5, cex=.8)
  if(h=='history'){
    h='pristine'
  }
  legend("topleft", c(paste("History:  ",h),
                      paste("Dots:      ",d),
                      paste("Views:    ",v),
                      paste("N obs:    ",n)), cex=.75, bty='n')
  if(!is.null(fn)) dev.off()
}

# Plot all models
par(mfrow=c(6,5), mar=c(2,2,1,1), oma=c(0,0,0,0))
for(i in 1:length(model))
  plot_thread(i)

# Plot one model
layout(1)
plot_thread(which_model("pristine",55,9))


#########################################################################################################################
#     5. Plot social info vs log-errors
#########################################################################################################################
load("data/models.rda")

plot_info <- function(i, Mlim=c(-2,5), Xlim=c(-2,5), fn=NULL){
  M.tmp = model[[i]]$weighted
  dat   = model[[i]]$data

  wB = as.matrix(M.tmp$sweight)%*%M.tmp$B$par

  Bcols = wbeta_colors(wB)
  h = unique(model[[i]]$data$method)
  d = unique(model[[i]]$data$d)
  v = unique(model[[i]]$data$v)
  n = nrow(model[[i]]$data)
  s = model[[i]]$session

  if(!is.null(fn)) pdf(fn, width=12, height=10)
      plot(dat$x1, dat$y, col=add_alpha(Bcols,.8), pch=16, cex.axis=1, bty='n', ann=FALSE, xlim=Mlim, ylim=Mlim)
      for(i in 1:M.tmp$states){
        a = M.tmp$mu$par[i]
        b = M.tmp$B$par[i]
        # abline(a,b, col=add_alpha('black',.5), lwd=3)
      }
      abline(v=0,h=0,a=0,b=1, lty=3)
      mtext(expression("log-ratio ("~y[i]~")"), side = 2, line = 1.8, cex=.8, las=0)
      mtext(expression("Social information ("~x[i]~")"), side = 1, line = 2, cex=.8)
      mtext(paste("Session no.",i,"id:",s), 3, line=0.5, cex=.8)
      if(h=='history'){
        h='pristine'
      }
      legend("topleft", c(paste("History:  ",h),
                          paste("Dots:      ",d),
                          paste("Views:    ",v),
                          paste("N obs:    ",n)), cex=.75, bty='n')
    if(!is.null(fn)) dev.off()

}


# Plot all models
par(mfrow=c(6,5), mar=c(2,2,1,1), oma=c(0,0,0,0))
for(i in 1:length(model))
  plot_info(i)

# Plot one model
layout(1)
plot_info(which_model("pristine",55,9))


#########################################################################################################################
#     6. Plot beta vs social information
#########################################################################################################################
load("data/models.rda")

plot_socialbeta <- function(type, fn=NULL){
  X = Y = cols = numeric(0)
  X.max = Y.max = cols.max = numeric(0)
  X.prs = Y.prs = cols.prs = numeric(0)
  for(i in 1:length(model)){
    M.tmp = model[[i]]$weighted
    dat   = model[[i]]$data
    wB    = as.matrix(M.tmp$sweight)%*%M.tmp$B$par
    Bcols = wbeta_colors(wB)

    X = c(X,dat$x1)
    Y = c(Y,wB)
    cols = c(cols, Bcols)
    if(all(dat$method=='max')){
      X.max = c(X.max,dat$x1)
      Y.max = c(Y.max,wB)
      cols.max = c(cols.max, Bcols)
    } else {
      X.prs = c(X.prs,dat$x1)
      Y.prs = c(Y.prs,wB)
      cols.prs = c(cols.prs, Bcols)
    }
  }
  tmp = list(x.all=X,
             y.all=Y,
             x.pristine=X.prs,
             y.pristine=Y.prs,
             x.manipulated =X.max,
             y.manipulated=Y.max,
             c.all = cols,
             c.pristine = cols.prs,
             c.manipulated = cols.max)

  if(!is.null(fn)) pdf(fn, width=12, height=10)
      if(type=='all'){
        plot(tmp$x.all,tmp$y.all, col=add_alpha(tmp$c.all,.8), pch=16, cex.axis=1, bty='n', ann=FALSE, xlim=c(-4,6), ylim=c(-1,1.5))
      } else if(type=='pristine'){
        plot(X.prs,Y.prs, col=add_alpha(cols.prs,.8), pch=16, cex.axis=1, bty='n', ann=FALSE, xlim=c(-4,6), ylim=c(-1,1.5))
      } else if(type=='manipulated'){
        plot(X.max,Y.max, col=add_alpha(cols.max,.8), pch=16, cex.axis=1, bty='n', ann=FALSE, xlim=c(-4,6), ylim=c(-1,1.5))
      }
      mtext("Social info", 1, 2)
      mtext(expression(beta[w]), 2, 2)
  if(!is.null(fn)) dev.off()
}

par(mfrow=c(3,1), mar=c(3,3,1,1))
plot_socialbeta('all')
plot_socialbeta('pristine')
plot_socialbeta('manipulated')



#########################################################################################################################
#     7. Plot beta scales
#########################################################################################################################

plot_beta_scale <- function(align, fn=NULL){
  cols   = colorRampPalette(c("red","chartreuse3","blue"), space="rgb", interpolate="linear")(200)
  midx   = seq(0,.8, length=200)
  cidx   = c(seq(-.5,min(midx),length=50),midx,seq(max(midx),1.5,length=50))
  Bscale = c(rep(cols[1],50),cols,rep(cols[length(cols)],50))

  if(!is.null(fn)) pdf(file = fn, width = 2, height=8)
      par(mfrow=c(1,1), mar=c(2,0,0,0), oma=c(0,0,0,0), cex.axis=1.2)
      if(align=="vert"){
        # Vertical beta-scale

        plot(rep(0, length(Bscale)),  cidx, col=add_alpha(Bscale,1), pch=15, cex=3, lwd=1, bty='n', ann=FALSE, axes=FALSE)
        axis(side = 4, at = pretty(cidx), cex=1.8, las=1, line = -3)
        mtext(expression(tilde(beta)~"scale"), side = 1, line=1, las=1, cex=1.8)
      } else{
        # Horizonal beta-scale
        plot(cidx,rep(0, length(Bscale)), col=add_alpha(Bscale,1), pch=15, cex=3, lwd=1, bty='n', ann=FALSE, axes=FALSE)
        axis(side = 1, at = pretty(cidx), cex=1.8, las=1, line = -3)
        mtext(expression(tilde(beta)~"scale"), side = 1, line=1, las=1, cex=1.8)
      }
  if(!is.null(fn)) dev.off()
}
plot_beta_scale('vert')
plot_beta_scale('hori')



