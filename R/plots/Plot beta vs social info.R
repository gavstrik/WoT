library(WoT)
# misc::clean_up()
load("data/models.rda")
hidx = names(model) %in% sessions$session[which(sessions$method=='history')]
midx = which(!hidx)
hidx = which(hidx)

x = y = numeric(0)
for(i in hidx){
  M.tmp = model[[i]]$weighted
  wB    = as.matrix(M.tmp$sweight)%*%M.tmp$B$par
  dat   = model[[i]]$data
  Si    = log(dat$wgmean/dat$d)
  x = c(x,Si)
  y = c(y,wB)
}
cols = wbeta_colors(y)
# par(mfrow=c(2,1),mar=c(3.5,3.5,1,1))

fn = "plots/historical_beta_social_info.pdf"
pdf(file=fn, width=10, height=6)
  par(mfrow=c(1,1),mar=c(3.5,3.5,1,1))
  plot(x,y,col = add_alpha(cols,.75), pch=16, bty='n', ann=FALSE, ylim=c(-.5,1), xlim=c(-4,4))
  mtext(expression("Social information ("*S[i]*")"), side=1, line=2.25)
  mtext(expression("Estimated score ("*beta[w]*")"), side=2, line=2.25)
dev.off()

x = y = numeric(0)
for(i in midx){
  M.tmp = model[[i]]$weighted
  wB    = as.matrix(M.tmp$sweight)%*%M.tmp$B$par
  dat   = model[[i]]$data
  Si    = log(dat$wgmean/dat$d)
  x = c(x,Si)
  y = c(y,wB)
}
cols = wbeta_colors(y)
fn = "plots/manipulated_beta_social_info.pdf"
pdf(file=fn, width=10, height=6)
    par(mfrow=c(1,1),mar=c(3.5,3.5,1,1))
    plot(x,y,col = add_alpha(cols,.75), pch=16, bty='n', ann=FALSE, ylim=c(-.5,1), xlim=c(-4,4))
    mtext(expression("Social information ("*S[i]*")"), side=1, line=2.25)
    mtext(expression("Estimated score ("*beta[w]*")"), side=2, line=2.25)
dev.off()
