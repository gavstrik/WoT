# library(misc)
library(MASS)
# library(depmixS4)
library(WoT)


# usefull function
clean_up <- function(...){
  for(i in dev.list()) dev.off(which=i)
  keep = as.list(match.call())
  rm(list=setdiff(ls(envir = .GlobalEnv), keep), envir = .GlobalEnv)
  # assign("savePDF", FALSE, envir = .GlobalEnv)
  assign("icloud", "~/Library/Mobile Documents/com~apple~CloudDocs/", envir = .GlobalEnv)
  cat("\014")  # Clear console
}

# usefull functions
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

# set_cobalt_plot_bg <- function(getcol=FALSE){
#  par(bg="#002240", col.lab="#FFFFFF", col.axis="#FFFFFF", col.main="#FFFFFF", col.sub="#FFFFFF", col="#FFFFFF", fg="#FFFFFF")
#  if(getcol){
#    print("#002240")
#  }
#}

#set_normal_plot_bg <- function(){
#  par(bg="#FFFFFF", col.lab="#000000", col.axis="#000000", col.main="#000000", col.sub="#000000", col="#000000", fg="#000000")
  

dat = WoT::newdat
dat$info = apply(newdat[,5:13], 1, gmean)

# Remove first guesses in threads with observations, since there was no prior info
dat = dat[!(dat$nobs!=0 & is.na(dat$info)),]
# dat = dat[dat$dots!=1233,] # no ox data

dat$X = log(dat$guess/dat$dots)    # response is the log of the relative error
dat$M = log(dat$info/dat$dots)     # social information is the log of the geo. average og previous estimates
dat   = dat[which(!is.na(dat$X)),] # remove any NA's in the response

nstates = 3
nobs    = 9
dots    = 1097
hist    = 'h'
hmm.dat = dat[dat$hist==hist & dat$nobs==nobs & dat$dots==dots, ]
hmm.dat = hmm.dat[hmm.dat$X < quantile(hmm.dat$X,.95) & hmm.dat$X > quantile(hmm.dat$X,.02),] 


nrow(hmm.dat)
M = get_hmm(hmm.dat, nstates, optim_state = TRUE, use_info = TRUE)


layout(1:3, heights = c(3,1,1))
    par(bty='n')
    set_cobalt_plot_bg()
    # scols = c("dodgerblue","chartreuse3","darkorange", "magenta")
    scols = RColorBrewer::brewer.pal(max(M$nstates,3),"Set1")
    
    # Colors as interpolated values between the two most probable states
      # tmp = M$states[,-1]
      # w   = apply(tmp,1,function(x) nmax(x,1)/sum(nmax(x,1:2)) )
      # c1  = scols[unlist(apply(tmp,1,function(x) which(x==nmax(x,1)) ))]
      # c2  = scols[unlist(apply(tmp,1,function(x) which(x==nmax(x,2)) ))]
      # cols = numeric(0)
      # for(i in 1:length(w)){
      #   ctmp    = colorRampPalette(c(c1[i],c2[i]))(1000)
      #   cols[i] = ctmp[round(w[i]*1000)]
      # }
      # cols = add.alpha(cols,.75)
    
    
    # Colors as intensity of most probable state
      pchs = rep(16,length(M$states$state))
      w    = apply(M$states[,-1],1,max)
      cols = scols[M$states$state]
      for(i in 1:length(cols)){
        cols[i] = add.alpha(cols[i],w[i])  
      }
      
      # idx = which(apply(M$states[,-1],1,min)>.20)
      # cols[idx] = 'red'
    
    # Colors as most probable state, when highly confident, grey when not
      # w    = apply(tmp,1,max)
      # cols = scols[M$states$state]
      # conf = 0.75
      # 
      # pchs = rep(17,length(M$states$state))
      
    
    # for(i in 1:length(cols)){
    #   cols[i] = add.alpha(cols[i],.75)  
    #   
    #   if(w[i]>.5){
    #     # cols[i] = add.alpha(cols[i],.75)  
    #     pchs[i] = 18
    #   } else if(w[i]>.75){
    #     # cols[i] = add.alpha(par()$fg,.75)
    #     pchs[i] = 15
    #   } else if(w[i]>.9){
    #     # cols[i] = add.alpha(par()$fg,.75)
    #     pchs[i] = 16
    #   }
    # }
    # 
    
   
par(mar=c(2,5,1,1))    
    
    
    # Plot the observed log relative errors
      plot(hmm.dat$X, col=add.alpha(par()$fg,.5), type='m', lty=3, ylim=c(-2,4))
      points(hmm.dat$X, col=cols, pch=pchs, cex=2)
      abline(h=0, lty=3)
    
    # Plot the estimated states
      plot(map_states(M$states$state,order(M$mu)), col=add.alpha(par()$fg,.5), type='l', lty=3, yaxt='n', cex=2)
      points(map_states(M$states$state,order(M$mu)), col=scols[M$states$state], pch=pchs)
      axis(2, at=1:M$nstates, las=1)#, labels = c("lo",'mid',"hi"))
    
    # Plot the data distribution with normal components and weighted average 
      hist(hmm.dat$X, breaks=50, prob=TRUE, border=NA, col=add.alpha(par()$fg,.25), xlim=c(-3,3), ann=FALSE)
      for(i in 1:M$nstates){
        curve(M$delta[i]*dnorm(x,M$mu[i],M$sig[i]),add=TRUE, col=scols[i], lty=1, lwd=2)
      }
      x = seq(-10,10,.01)
      y = numeric(length(x))
      for(i in 1:length(M$mu)){
        y = y+M$delta[i]*dnorm(x,M$mu[i],M$sig[i])  
      }
      lines(x,y, col=par()$fg,lwd=2,lty=1)
    
      # m = median(hmm.dat$X)
      # s = diff(quantile(hmm.dat$X, c(.25,.75)))/2
      # curve(dcauchy(x,m,s), add=TRUE, col='black',lwd=2, lty=3)
      # 
      # m = median(hmm.dat$X)
      # s = diff(quantile(hmm.dat$X, c(.25,.75)))/2
      # curve(dnorm(x,m,s), add=TRUE, col='black',lwd=2, lty=2)
      
      # legend("topright", c("Gaussian","Cauchy","Mixed Gaussian", "Components"), lty=c(1,1,1,3), lwd=c(2,2,2,1), col=c("orange","dodgerblue","red","black"), bty='n')
    
      
    idx = order(M$mu)
    round(M$info[idx],3)
    round(M$mu[idx],3)
    round(M$sig[idx],3)
    round(M$aic[idx,],3)
    
    round(table(M$states$state)/length(M$states$state),3)[idx]
    
    round(M$delta[idx],3)
    round(M$Gam[idx,idx],3)
    
      mtext(paste("Type:", hist, "Views:",nobs,"Dots:",dots), 3, outer=TRUE, line=-2)

      