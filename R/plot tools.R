wbeta_colors <- function(wB, neutral = 0.4){
  # Make colors of weighted betas, with neutral at 0.4
  cols = colorRampPalette(c("red","chartreuse3","blue"), space="rgb", interpolate="linear")(200)

  midx = seq(0,2*neutral, length=length(cols))

  out = numeric(length(wB))
  for(i in 1:length(wB)){
    out[i] = cols[(midx-wB[i])^2 == min((midx-wB[i])^2)]
  }
  return(out)
}

beta_colors <- function(M, neutral = 0.4){
  # Make colors of weighted betas, with neutral at 0.4
  W = as.matrix(M$sweight)
  B = M$B$par
  
  wB = W%*%B
  # idx = order(wB)
  # y = seq(0,1,length=length(wB))
  # x = wB[idx]
  # n = length(x)
  # 
  # cols = colorRampPalette(c("red","chartreuse3","blue"), space="rgb", interpolate="linear")(200)
  # 
  # midx = seq(0,2*neutral, length=length(cols))
  # 
  # out = numeric(length(wB))
  # for(i in 1:length(wB)){
  #   out[i] = cols[(midx-wB[i])^2 == min((midx-wB[i])^2)]
  # }
  out = wbeta_colors(wB)
  return(out)
}


plot_beta <- function(M.tmp, cols, probs=NULL,xlim=NULL, noCov=FALSE, ann=FALSE){
  # Plot weighted betas
  V = vcov(M.tmp$fit)$vcov
  rownames(V) = rownames(M.tmp$pars)[!M.tmp$pars$approx]
  colnames(V) = rownames(M.tmp$pars)[!M.tmp$pars$approx]


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

  segments(x.hi[idx.hi],y,x.lo[idx.lo],y, col=add.alpha(cols,.5),lwd=2)
  segments(x[-n],y[-n],x[-1],y[-1], col=add.alpha(cols,.85), lwd=2)

  # mtext("CDF",side = 2, line=2.5, cex=1.5)
  if(ann){
    mtext(expression(tilde(beta)),side = 1, line=3, cex=.75)
    }

  # plot(wB[idx],y, type='l', bty='n', lwd=2, xlim=c(0,1))
  # polygon(c(wB.hi[idx],rev(wB.lo[idx])), c(y,rev(y)), border=NA, col=add.alpha('grey',.5))
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



