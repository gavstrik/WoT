#########################################################################################################################
#   Script for the Gaussian mixture model analysis of dots.rda used in the paper:
#     "The Effects of Bounded Social Information on Sequential Decision Making" - Engelhardt et al (2021)
#
#   Author:   Jacob Stærk-Østergaard,
#   email:    ostergaard@sund.ku.dk
#
#########################################################################################################################

library(WoT)
# misc::clean_up() # cleaning the workspace

#########################################################################################################################
#     1. Load data and create variables
#########################################################################################################################

# Find the geometric (including weighted) mean of available social information
dots <- read.csv("data/dots.csv", colClasses = c("character",
                                                 "numeric",
                                                 "numeric",
                                                 "character",
                                                 "numeric",
                                                 "character",
                                                 "numeric") , stringsAsFactors=FALSE)

dots$hist = substr(dots$hist, 2, nchar(dots$hist)-1)
tmp = lapply(strsplit(dots$hist,", "), as.numeric)
names(tmp) = dots$d

# Make a function to weigh previous estimates according to the control groups (v=0)
v0_df <- function(s){
  tmp = subset(dots, session == s)
  d   = unique(tmp$d)
  out = approxfun(density(tmp$guess,from=0, to=100*d, bw = 10*d) ) # wide bandwidth ensures almost equal weights to "reasonable" estimates
}

f = list(d55   = v0_df("fscmakcz"), # v=0 d=55 session
         d148  = v0_df("dwnjf9mb"), # v=0 d=148 session
         d403  = v0_df("hqx0v7t5"), # v=0 d=403 session
         d1097 = v0_df("hal5jdl0")) # v=0 d=1097 session


w = list()
for(i in 1:length(tmp)){
  d = as.numeric(names(tmp[i]))
  h = tmp[[i]]
  g = f[[which(d == c(55,148,403,1097))]]
  w[[i]] = g(h)/sum(g(h), na.rm=TRUE)
  # k = g(h)
  # if(any(is.na(k))){
  #   k[which(is.na(k))] = 1e-9
  # }
  # w[[i]] = k/sum(k)

  dots$gmean[i] = gmean(h)
  dots$wgmean[i] = gmean(h,w[[i]])
}

cat("\nTotal number of datapoints:", nrow(dots),
    "\nNumber of unique sessions: ", length(unique(dots$session)))

fn = "data/dots.rda"
save(dots, file=fn)

#########################################################################################################################
#     2. Remove outliers
#########################################################################################################################

      # Find the 1%, 2.5%, 97.5% and 99% quantiles of each session as well as .1d and 10d extremes

        q01  = aggregate(guess~session, data=dots, function(x) quantile(x,.01))
        q025 = aggregate(guess~session, data=dots, function(x) quantile(x,.025))
        q975 = aggregate(guess~session, data=dots, function(x) quantile(x,.975))
        q99  = aggregate(guess~session, data=dots, function(x) quantile(x,.99))
        names(q01)[2] = "q01"
        names(q025)[2] = "q025"
        names(q975)[2] = "q975"
        names(q99)[2] = "q99"

        dots = dplyr::left_join(dots,q01)
        dots = dplyr::left_join(dots,q025)
        dots = dplyr::left_join(dots,q975)
        dots = dplyr::left_join(dots,q99)

        dots$d01 = dots$d/10
        dots$d10 = dots$d*10

      # Find the number of data points left within the three defined boundaries
        idx1 = dots$guess>=dots$q01 & dots$guess<=dots$q99
        idx2 = dots$guess>=dots$q025 & dots$guess<=dots$q975
        idx3 = dots$guess>=dots$d01 & dots$guess<=dots$d10

        tmp  = aggregate(guess ~ method+d+v+session, data=dots, FUN = length)
        tmp1 = aggregate(guess ~ method+d+v+session, data=dots[idx1,], FUN = length)
        tmp2 = aggregate(guess ~ method+d+v+session, data=dots[idx2,], FUN = length)
        tmp3 = aggregate(guess ~ method+d+v+session, data=dots[idx3,], FUN = length)

      # Compare with full number of observations
        names(tmp)[5] = "all"
        tmp1 = dplyr::left_join(tmp,tmp1)
        tmp2 = dplyr::left_join(tmp,tmp2)
        tmp3 = dplyr::left_join(tmp,tmp3)

        tmp1 = round((tmp1$all-tmp1$guess)/tmp1$all,2)*100
        tmp2 = round((tmp2$all-tmp2$guess)/tmp2$all,2)*100
        tmp3 = round((tmp3$all-tmp3$guess)/tmp3$all,2)*100

        tmp = cbind(tmp,tmp1,tmp2,tmp3)
    # How many % of observations are removed with either cutoff method
        tmp = tmp[order(tmp$method, tmp$d, tmp$v),]
        names(tmp) = c("method","d","v","session","obs","Pct removed at 2%","Pct removed at 5%","Pct removed by a factor 10")
        print(tmp)
        sessions = tmp

        # Conclusions:  The factor 10 wrong removes far too many observations in the manipulated series in some cases >20%.
        #               We choose to remove 5% extreme guesses, since a) we are still left with the majority of data (95%)
        #               and removal is fairly consistent at 5% for series with more that 100 obs.

    # Dots for analysis is where guess is within quantiles [.01, .99]:
        cat("\nPct used for analysis:", nrow(dots[idx2,])/nrow(dots)) # Should be approx. .95
        dots = dots[idx2,]

    # Fit GMMs to dots data where gmean or wgmean is not missing, i.e. neither first guess or if v=0 or where social information is extreme beyond reason
        idx = which(!is.na(dots$gmean) & !is.na(dots$wgmean))
        dots = dots[idx,]


        fitsessions = subset(tmp, v!=0 & obs>25) # fit sessions with v>0 and more than 25 observations
        # fitsessions

#########################################################################################################################
#     3. Fit GMMs
#########################################################################################################################

      # Variables to include

        dots$y  = log(dots$guess/dots$d)
        dots$x0 = log(dots$gmean/dots$d)
        dots$x1 = log(dots$wgmean/dots$d)


        h <- function(i){
          fit.dat = subset(dots, session==fitsessions$session[i])
          full = fit_gmm(y~x0, fit.dat, opt = 'bic', trystates = 2:5, maxiter = 10)
          weighted = fit_gmm(y~x1, fit.dat, opt = 'bic', trystates = 2:5, maxiter = 10)
          M.out = list(session = fitsessions$session[i], full=full, weighted=weighted)
          return(M.out)
        }

        set.seed(12)
        # depmixS4::em.control(maxit=5e3)

        # Fit using parallel threads fo speed up (takes approx 100 seconds, 30 seconds with paralell).

        system.time({
          model = lapply(1:nrow(fitsessions),h)
        })
        # If parallel is installed:
        # system.time({
        #   model = parallel::mclapply(1:nrow(fitsessions),f,mc.cores = parallel::detectCores())
        # })

        names(model) = unlist(lapply(model, function(x) x$session))


#########################################################################################################################
#     4. Check GMM convergence
#########################################################################################################################
          # sessions[sessions$session %in% names(model)[idx],]
          # plot_session(3, rng=1)  # Session 3 shows no significant relation between guess and info, hence 2 states is sufficient
          # plot_session(23, rng=3) # Session 23 shows no significant relation between guess and info, hence 2 states is sufficient
          # Check criteria:
          # tmp = lapply(model, function(x) x$weighted$fit)
          lapply(model, function(x) x$weighted$fit)
          lapply(model, function(x) x$weighted$crit)



          # for(i in idx){
          #   fit.dat = subset(dots, session==fitsessions$session[i])
          #   full = fit_gmm(y~x0, fit.dat, opt = 'bic', trystates = 2, maxiter = 10)
          #   weighted = fit_gmm(y~x1, fit.dat, opt = 'bic', trystates = 2, maxiter = 10)
          #   model[[i]] = list(session = fitsessions$session[i], full=full, weighted=weighted)
          #   fit.problems$weighted[i] = any(is.na(model[[i]]$weighted$crit$logLik))
          # }

sessions$fitted = sessions$session %in% fitsessions$session
# cat("\nAny remaining problems regarding fit? ",any(fit.problems$weighted),"\nFixed sessions: \n",paste(fitsessions$session[idx],collapse = "\n"), sep="")


#########################################################################################################################
#     5. Extract parameters and summaries
#########################################################################################################################

        for(i in 1:nrow(fitsessions)){
          M.tmp = model[[i]]

          M.tmp$full = append(M.tmp$full, WoT::getParameters(M.tmp$full$fit, se.fit = TRUE))  # append parameter estimates and std.errs
          tmp2       = depmixS4::posterior(M.tmp$full$fit)                                    # append estimated states
          M.tmp$full = append(M.tmp$full,list(state=tmp2$state, sweight = tmp2[,-1]))

          M.tmp$weighted = append(M.tmp$weighted, WoT::getParameters(M.tmp$weighted$fit, se.fit = TRUE))       # append parameter estimates and std.errs
          tmp2       = depmixS4::posterior(M.tmp$weighted$fit)                                    # append estimated states
          M.tmp$weighted = append(M.tmp$weighted,list(state=tmp2$state, sweight = tmp2[,-1]))

          tmp   = list(data=dots[dots$session==model[[i]]$session,])                                             # append dataset
          M.tmp = append(M.tmp,tmp)

          model[[i]] = M.tmp
        }

    # Check whether standard errors are present or must be approximated
        # cat("\014")
        se.problems = data.frame(session=names(model), full=NA, weighted=NA)
        for(i in 1:nrow(fitsessions)){
          se.problems$full[i] = any(model[[i]]$full$pars$approx)
          se.problems$weighted[i] = any(model[[i]]$weighted$pars$approx)
        }
  cat("\nTroubled SEs for geometric mean:         ",which(se.problems$full),
      "\nTroubled SEs for weighted geometric mean:",which(se.problems$weighted))

        idx = which(se.problems$weighted)
        # fitsessions[idx,]

        for(i in idx){
          pidx = which(model[[i]]$weighted$pars$approx)
          cat("\n",rownames(model[[i]]$weighted$pars[pidx,]))
        }

    miss = lapply(model[idx], function(x) {
                                    pidx = x$weighted$pars$approx
                                    return(rownames(x$weighted$pars[pidx,]))
                                  })

    cat("\nMissing SEs:", unique(unlist(miss)))

    # No betas are missing SEs, hence we continue without any further issues
        fn = "data/models.rda"
        save(model, file=fn)

        fn = "data/sessions.rda"
        save(sessions, file=fn)

