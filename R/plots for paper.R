wdir = dirname(rstudioapi::getActiveDocumentContext()$path)
wdir = paste0(wdir,"/")

source(paste0(wdir,"tools.R"))
source(paste0(wdir,"plot tools.R"))


fn = paste0(wdir,"data/models.Rda")
load(fn)

fn = paste0(wdir,"data/sessions.Rda")
load(fn)

fn   = paste0(wdir,"data/dots.Rda")
load(fn)


plot_empty <- function(){
  plot(0,0, ann=FALSE, axes=FALSE, type='n')
}

plot_res <- function(i, info=TRUE){
  M.tmp = model[[i]]$weighted
  dat   = model[[i]]$data
  res   = get_res(M.tmp,dat)
  cols = add.alpha('black',.5)

  d = unique(dat$d)
  v = unique(dat$v)
  h = unique(dat$method)
  s = model[[i]]$session
  n = nrow(dat)
  qqnorm(res, pch=16, col=cols, bty='n', xlim=c(-4,4), ylim=c(-4,4), main="", ann=FALSE, cex.main=.9)

  abline(0,1, lty=3)

  if(info){
    legend("topleft", c(paste("History:  ",h),
                        paste("Dots:      ",d),
                        paste("Views:    ",v),
                        paste("N obs:    ",n)), cex=1, bty='n')
  }
}

plot_betadist <- function(i, info=TRUE){

  s = sessions$session[sessions$v>0]
  idx = which(names(model) == s[i])
  M.tmp = model[[idx]]$weighted
  Bcols = beta_colors(M.tmp)
  dat   = model[[idx]]$data
  h = unique(model[[idx]]$data$method)
  d = unique(model[[idx]]$data$d)
  v = unique(model[[idx]]$data$v)
  n = nrow(model[[idx]]$data)

  plot_beta(M.tmp,Bcols,noCov = FALSE, xlim=c(-.5,1.5), ann=FALSE)

  if(info){
    legend("topleft", c(paste("History:  ",h),
                        paste("Dots:      ",d),
                        paste("Views:    ",v),
                        paste("N obs:    ",n)), cex=1, bty='n')
  }
}

plot_info <- function(i, info=TRUE){
  Mlim = c(-2,5)
  Xlim = c(-2,5)

  s = sessions$session[sessions$v>0]
  idx = which(names(model) == s[i])
  M.tmp = model[[idx]]$weighted
  Bcols = beta_colors(M.tmp)
  dat   = model[[idx]]$data
  h = unique(model[[idx]]$data$method)
  d = unique(model[[idx]]$data$d)
  v = unique(model[[idx]]$data$v)
  n = nrow(model[[idx]]$data)


  plot(dat$x1, dat$y, col=add.alpha(Bcols,.8), pch=16, cex.axis=1, bty='n', ann=FALSE, xlim=Mlim, ylim=Mlim)
  # for(i in 1:M.tmp$states){
  #   a = M.tmp$mu$par[i]
  #   b = M.tmp$B$par[i]
  #   abline(a,b, col=add.alpha('black',.5), lwd=3)
  # }
  abline(v=0,h=0,a=0,b=1, lty=3)
  # mtext(expression("log-ratio ("~y[i]~")"), side = 2, line = 1.8, cex=1, las=0)
  # mtext(expression("Social information ("~s[i]~")"), side = 1, line = 2, cex=1)

  if(info){
    legend("topleft", c(paste("History:  ",h),
                        paste("Dots:      ",d),
                        paste("Views:    ",v),
                        paste("N obs:    ",n)), cex=1, bty='n')
  }
}

plot_thread <- function(i, info=TRUE){
  s = sessions$session[sessions$v>0]
  idx = which(names(model) == s[i])
  M.tmp = model[[idx]]$weighted
  Bcols = beta_colors(M.tmp)
  dat   = model[[idx]]$data
  h = unique(model[[idx]]$data$method)
  d = unique(model[[idx]]$data$d)
  v = unique(model[[idx]]$data$v)
  n = nrow(model[[idx]]$data)

  Mlim = c(-2,5)
  Xlim = c(-2,5)

  plot(dat$y, col=add.alpha(Bcols,.8), pch=16,cex.axis=1, bty='n', ylim=Xlim, ann=FALSE, xlim=c(0,500))
  lines(dat$x1, lty=3, col=add.alpha('black',.5))
  # mtext(expression("log-ratio ("~y[i]~")"), side = 2, line = 1.8, cex=1, las=0)
  # mtext("Observation no.", side = 1, line = 2.25, cex=1)
  abline(h=0, col=add.alpha('black',.3))

  if(info){
    legend("topleft", c(paste("History:  ",h),
                        paste("Dots:      ",d),
                        paste("Views:    ",v),
                        paste("N obs:    ",n)), cex=1, bty='n')
  }
}

plot_betascale <- function(layout="vertical", info=TRUE, line=-1){

  cols   = colorRampPalette(c("red","chartreuse3","blue"), space="rgb", interpolate="linear")(200)
  midx   = seq(0,.8, length=200)
  cidx   = c(seq(-.5,min(midx),length=50),midx,seq(max(midx),1.5,length=50))
  Bscale = c(rep(cols[1],50),cols,rep(cols[length(cols)],50))


  if(layout == "vertical"){
    plot(rep(0, length(Bscale)),  cidx, col=add.alpha(Bscale,1), pch=15, cex=3, lwd=1, bty='n', ann=FALSE, axes=FALSE)
    axis(side = 4, at = pretty(cidx), cex=1.8, las=1, line = line)

    if(info){
      mtext(expression(beta^w~"scale"), side = 4, line=1, las=0, cex=1.4)
    }
  } else{
    plot(cidx,rep(0, length(Bscale)), col=add.alpha(Bscale,1), pch=15, cex=3, lwd=1, bty='n', ann=FALSE, axes=FALSE)
    axis(side = 1, at = pretty(cidx), cex=1.8, las=1, line = line)
    if(info){
      mtext(expression(beta^w~"scale"), side = 1, line=1, las=1, cex=1.4)
    }

  }
}



# Sessions to plot
# sessions
s = sessions$session[sessions$v>0]


wdir
# Plot Fig. 6-9

fn = paste0(wdir,"plots/qqplots.pdf")
pdf(fn, width=10, height=12)
    layout(matrix(c(1:30), nr=6, byrow=TRUE))
    par(mar=c(3,3,0,0), oma=c(3,3,0,0))
    for(i in 1:length(s)){
      plot_res(i)
    }
    mtext("Theoretical Quantiles", side = 1, line=1, outer=TRUE)
    mtext("Sample Quantiles", side = 2, line=1, outer=TRUE)
dev.off()

fn = paste0(wdir,"plots/allbetas.pdf")
pdf(fn, width=10, height=13)
    layout(matrix(c(1:30,30,31,31,31,30), nr=7, byrow=TRUE), heights = c(rep(4,6),1))
    par(mar=c(3,3,0,0), oma=c(4,3,0,0.25))
    for(i in 1:length(s)){
      plot_betadist(i)
    }

    plot_empty()
    par(mar=c(0,0,2,0))
    plot_betascale("horizontal", FALSE,line = 0)
    mtext("CDF",side = 2, line=1, cex=1, outer=TRUE)
    mtext(expression(beta^w),side = 1, line=-3.25, cex=1, outer=TRUE)

dev.off()

fn = paste0(wdir,"plots/allthreads.pdf")
pdf(fn, width=10, height=13)
    layout(matrix(c(1:30,30,31,31,31,30), nr=7, byrow=TRUE), heights = c(rep(4,6),1))
    par(mar=c(3,3,0,0), oma=c(4,3,0,0.25))
    for(i in 1:length(s)){
      plot_thread(i)
    }

    plot_empty()
    par(mar=c(0,0,2,0))
    plot_betascale("horizontal", FALSE, line=0)
    mtext(expression("log-ratio ("~y[i]~")"),side = 2, line=1, cex=1, outer=TRUE)
    mtext("Observation no.",side = 1, line=-3.25, cex=1, outer=TRUE)
dev.off()


fn = paste0(wdir,"plots/allinfo.pdf")
pdf(fn, width=10, height=13)
    layout(matrix(c(1:30,30,31,31,31,30), nr=7, byrow=TRUE), heights = c(rep(4,6),1))
    par(mar=c(3,3,0,0), oma=c(4,3,0,0.25))
    for(i in 1:length(s)){
      plot_info(i)
    }

    plot_empty()
    par(mar=c(0,0,2,0))
    plot_betascale("horizontal", FALSE, line=0)
    mtext(expression("log-ratio ("~y[i]~")"),side = 2, line=1, cex=1, outer=TRUE)
    mtext(expression("Social information ("~s[i]~")"),side = 1, line=-3.25, cex=1, outer=TRUE)
dev.off()



# Plot Fig. 10

plot_pair <- function(h,d,v, infox=FALSE, infoy=FALSE){
  i = which(s==sessions$session[which(sessions$method==h & sessions$d==d & sessions$v==v)])
  plot_thread(i)
  if(infoy){
    mtext(expression("log-ratio ("~y[i]~")"),side = 2, line=2, cex=1)
  }
  if(infox){
    mtext("Observation no.",side = 1, line=2.5, cex=1)
  }


  plot_betadist(i, info=FALSE)
  if(infox){
    mtext(expression(beta^w),side = 1, line=2.5, cex=1)
  }

}


fn = paste0(wdir,"plots/thread_evolution.pdf")
pdf(fn, width=10, height=12)
    layout( matrix( c( 1, 1, 1, 2, 2, 3, 3, 3, 4, 4,
                       5, 5, 5, 6, 6, 7, 7, 7, 8, 8,
                       9, 9, 9,10,10,11,11,11,12,12,
                      13,13,13,14,14,15,15,15,16,16,
                      17,17,17,18,18,19,19,19,20,20,
                      21,21,22,22,22,22,22,22,21,21), nc=10, byrow=TRUE), heights = c(rep(4,5),1))
    par(mar=c(3,3,0,0), oma=c(4,3,0,0.25))
    plot_pair("history",d=55,v=1, infoy = TRUE)
    plot_pair("max",d=55,v=1)

    plot_pair("history",d=55,v=3, infoy=TRUE)
    plot_pair("max",d=55,v=3)

    plot_pair("history",d=148,v=3, infoy=TRUE)
    plot_pair("max",d=148,v=3)

    plot_pair("history",d=403,v=9, infoy=TRUE)
    plot_pair("max",d=403,v=9)

    plot_pair("history",d=1097,v=9, infoy=TRUE, infox=TRUE)
    plot_pair("max",d=1097,v=9, infox=TRUE)

    par(mar=c(0,0,3,0))

    plot_empty()
    plot_betascale("horizontal",FALSE, line=0)
dev.off()



# Plot Fig. 11



# Plot Fig. 3

plot_triple <- function(h,d,v, infox=FALSE, infoy=FALSE){
  i = which(s%in%sessions$session[which(sessions$method==h & sessions$d==d & sessions$v==v)])
  if(length(i)>1){
    i = i[1]
  }
  plot_thread(i)
  if(infoy){
    mtext(expression("log-ratio ("~y[i]~")"),side = 2, line=2, cex=1)
  }
  if(infox){
    mtext("Observation no.",side = 1, line=2.5, cex=1)
  }

  plot_info(i, info=FALSE)
  if(infoy){
    mtext(expression("log-ratio ("~y[i]~")"), side = 2, line = 2, cex=1, las=0)
  }
  if(infox){
    mtext(expression("Social information ("~s[i]~")"), side = 1, line = 2.5, cex=1)
  }

  plot_betadist(i, info=FALSE)
  if(infox){
    mtext(expression(beta^w),side = 1, line=2.5, cex=1)
  }

}

fn = paste0(wdir,"plots/threads.pdf")
pdf(fn, width=10, height=5)
    layout(matrix(c(1,1,4,4,7,
                    2,3,5,6,7), nc=5, byrow=TRUE), widths = c(4,4,4,4,2))
    par(mar=c(4,3,0,0), oma=c(1,1,0,0.25))
    plot_triple("history",1097,1, infox = TRUE, infoy = TRUE)
    plot_triple("max",1097,9, infox = TRUE, infoy = TRUE)
    par(mar=c(4,0,0,0))
    plot_betascale(line = -3, info = FALSE)
dev.off()




# Plot Fig. 5

dat = data.frame(session = dots$session, type = dots$method, d = dots$d, v = as.factor(dots$v), y = log(dots$guess/dots$d))
dat = aggregate(y ~ type+d+v, data=dat, median)

dat.h = subset(dat,type=='history')
dat.m = subset(dat,type=='max' | v==0)

M.h = glm(y ~ log(d)*v, data=dat.h)
M.m = glm(y ~ log(d)*v, data=dat.m)

res.h = residuals(M.h)
res.m = residuals(M.m)

res.h = (res.h-mean(res.h))/sd(res.h)
res.m = (res.m-mean(res.m))/sd(res.m)


plot_res <- function(res){
  qqnorm(res, pch=16, ann=FALSE, xlim=c(-2,2), ylim=c(-2,2))
  abline(0,1, lty=3)
  mtext("Theoretical",1, line=2)
  mtext(expression("Empirical"~hat(epsilon)),2, line=2)

  plot(res, pch=16, ylim=c(-2,2), xaxt='n', ann=FALSE)
  axis(side = 1, at = 1:4, labels=NA)
  axis(side = 1, at = 1:4+4, labels=NA)
  axis(side = 1, at = 1:4+8, labels=NA)
  axis(side = 1, at = 1:4+12, labels=NA)
  text(x = 1:16,par()$usr[3]*1.2, rep(c(55,148,403,1097),4), srt=60, cex=.75, xpd=TRUE)
  axis(side = 1, at = seq(2.5,16,4), labels=c(0,1,3,9), tick = FALSE, line = .75, cex=.75)
  mtext("Dots",side = 1, line = 1, adj = -0.35)
  mtext("Views",side = 1, line = 2, adj = -0.35)
  mtext(expression("Empirical"~hat(epsilon)),2, line=2)
}


fn = paste0(wdir,"plots/qqplots_median.pdf")
pdf(fn, width = 6, height = 6)
par(mfrow=c(2,2), mar=c(5,5,0,1), oma=c(0,0,0,0), bty='n')
  plot_res(res.h)
  plot_res(res.m)
dev.off()


# Plot Fig. 2


x    = seq(50,1100,5)
v = c(0,1,3,9)
cols = RColorBrewer::brewer.pal(4,"Set1")

plot_medians <- function(M, dat, ylim, info=TRUE){
  plot(0,0, xlim=log(c(55,1100)), ylim=ylim, type='n', ann=FALSE, las=1, xaxt='n', cex.axis=1, yaxt='n', bty='n')
  axis(side = 1, at = seq(4,7), cex.axis=1)
  axis(side = 2, at = seq(min(ylim),max(ylim),.25), labels = seq(min(ylim),max(ylim),.25), cex.axis=1, las=1)

  abline(h=0, lty=3)
  mtext("log(d)", 1, line=2, cex=1)
  mtext(expression("median log-ratio"), 2, line=3, cex=1)

  for(i in 1:length(v)){
    tmp.fit = predict(M, newdata = expand.grid(v=as.factor(v[i]), d=x), se.fit = TRUE)
    fit.lo = tmp.fit$fit-1.96*tmp.fit$se.fit
    fit.hi = tmp.fit$fit+1.96*tmp.fit$se.fit
    polygon(log(c(x,rev(x))), c(fit.lo,rev(fit.hi)), border=NA, col=add.alpha(cols[i],.6))
    lines(log(x),tmp.fit$fit, col=cols[i], lwd=2)
  }
  for(i in 1:length(v)){
    points(log(jitter(dat$d[dat$v==v[i]], .1)),jitter(dat$y[dat$v==v[i]],.1), col=cols[i], pch=16, cex=1.2)
  }

  if(info){
    legend("topright", legend=paste0("v=",rev(v)), col=rev(cols), pch=16, bty='n', cex=1, ncol = 2)
  }

}

fn = paste0(wdir,"plots/medians.pdf")
pdf(fn, width=10, height=4)
    par(mfrow=c(1,2), mar=c(3,4,0,1))
    plot_medians(M.h, dat.h, ylim=c(-.75,.5))
    plot_medians(M.m, dat.m, ylim=c(-.75,1.25), FALSE)
dev.off()








layout(1)


