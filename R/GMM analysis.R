library(depmixS4)

###### 1. Read data and prepare for analysis ##################################################################################################################################################################

        wdir = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/WoT/"

        wdir = dirname(rstudioapi::getActiveDocumentContext()$path)
        wdir = paste0(wdir,"/")


        source(paste0(wdir,"tools.R"))


        dots <- read.csv(paste0(wdir,"data/dots.csv"), colClasses = c("numeric",
                                                                      "character",
                                                                      "character",
                                                                      "numeric",
                                                                      "numeric",
                                                                      "character",
                                                                      "character",
                                                                      "numeric",
                                                                      "character",
                                                                      "numeric") , stringsAsFactors=FALSE)

        # Remove brackets [] from hist variable and make a list vector of history as numerical values
          dots$hist = substr(dots$hist, 2, nchar(dots$hist)-1)
          hist = lapply(strsplit(dots$hist,", "), as.numeric)
          cat("Total number of datapoints:", nrow(dots))

        # Find geometric mean using all available info
          dots$gmean = unlist(lapply(hist, gmean))

        # Make weighted estimates of geometric means, based on baseline densities for v=0
              dots$gmean.weighted = NA
              make_df <- function(s){
                # Make a density estimate function for session 's'
                tmp   = subset(dots, session == s)
                f.out = approxfun(density(tmp$guess,from=0, to=max(dots$guess) ) )

                return(f.out)
              }
              # tmp = aggregate(guess ~ method+d+v+session, data=dots[dots$v==0,], FUN = length) # Sessions with 0 views
              tmp = aggregate(guess ~ method+d+v+session, data=dots[dots$v!=0,], FUN = length) # Sessions with >0 views

              for(i in 1:nrow(tmp)){
                idx = which(dots$session == tmp$session[i])
                d   = unique(dots$d[idx])

                if(d == 55){
                  df = make_df("fscmakcz") # v=0 d=55 session
                } else if(d == 148){
                  df = make_df("dwnjf9mb") # v=0 d=148 session
                } else if(d == 403){
                  df = make_df("hqx0v7t5") # v=0 d=403 session
                } else if(d == 1097){
                  df = make_df("hal5jdl0") # v=0 d=1097 session
                } else{
                  warning("Session does not exist")
                }

                tmp.hist = unlist(lapply(hist[idx], function(x){
                  w = df(x)
                  if(sum(w)==0){
                    w = rep(1,length(w))
                  }
                  w = w/sum(w)
                  gmean(x, w)
                }))
                dots$gmean.weighted[idx] = tmp.hist
              }

      # check whether any weighted gmeans are missing ()
        id = which(is.na(dots$gmean.weighted) & dots$v!=0 & dots$decision_order>1)
        dots[id,]
      # check whether any gmeans are missing ()
        id = which(is.na(dots$gmean) & dots$v!=0 & dots$decision_order>1)
        dots[id,]

      # Set NaN values to NA
        dots$gmean[is.nan(dots$gmean)] = NA
        dots$gmean.weighted[is.nan(dots$gmean.weighted)] = NA

      # Keep relevant columns
        dots = dots[,-c(1,2,9)]

      # Store prepared data
        save(dots, file=paste0(wdir,"data/dots.Rda"))



###### 2. Fit models ##########################################################################################################################################################################################


      # Remove extreme guesses

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
    # How many % of observations are remove with either cutoff method
        tmp = tmp[order(tmp$method, tmp$d, tmp$v),]
        print(tmp)
        sessions = tmp

        # Conclusions:  The factor 10 wrong removes far too many observations in the manipulated series in some cases >20%.
        #               We choose to remove 5% extreme guesses, since a) we are still left with the majority of data (95%)
        #               and removal is fairly consistent at 5% for series with more that 100 obs.


    # Dots for analysis is where guess is within quantiles [.025, .975]:
        nrow(dots[idx2,])/nrow(dots) # Should be approx. .95
        dots = dots[idx2,]



    # Fit GMMs to dots data
        idx = !is.na(dots$gmean)
        dots = dots[idx,]

        fitsessions = subset(tmp, v!=0)
        fitsessions

        dots$y  = log(dots$guess/dots$d)
        dots$x0 = log(dots$gmean/dots$d)
        dots$x1 = log(dots$gmean.weighted/dots$d)

        f <- function(i, states=2:5){
          fit.dat = subset(dots, session==fitsessions$session[i])
          full = fit_gmm(y~x0, fit.dat, opt = 'bic', trystates = states, maxiter = 10)
          weighted = fit_gmm(y~x1, fit.dat, opt = 'bic', trystates = states, maxiter = 10)
          M.out = list(session = fitsessions$session[i], full=full, weighted=weighted)
          return(M.out)
        }

        set.seed(1234)
        # Fit using parallel threads fo speed up.
        model = parallel::mclapply(1:nrow(fitsessions),f,mc.cores = 4)
        names(model) = unlist(lapply(model, function(x) x$session))


        # Check results for convergence
          # cat("\014")
          fit.problems = data.frame(session=names(model), full=NA, weighted=NA)
          for(i in 1:29){
                # cat("------------------------------   ",i,"full       ----------------------------------------------------\n")
                # print(model[[i]]$full$crit)
                # cat("------------------------------   ",i,"weighted   ----------------------------------------------------\n")
                # print(model[[i]]$weighted$crit)
                # cat("\n")

                fit.problems$full[i] = any(is.na(model[[i]]$full$crit$logLik))
                fit.problems$weighted[i] = any(is.na(model[[i]]$weighted$crit$logLik))
            # Sys.sleep(1)
          }
          which(fit.problems$full)
          which(fit.problems$weighted)

          idx = which(fit.problems$weighted) # Problem threads, either few observations or simple relation between social info and estimate

          for(i in idx){
            model[[i]] = f(i,2)
          }

        fn = paste0(wdir,"data/models.Rda")
        save(model, file=fn)


        fn = paste0(wdir,"data/sessions.Rda")
        save(sessions, file=fn)



###### 3. Extract parameters and summaries ####################################################################################################################################################################

        # Check whether standard errors are present
        # cat("\014")
        se.problems = data.frame(session=names(model), full=NA, weighted=NA)
        for(i in 1:29){
          # cat("------------------------------   ",i,"   ----------------------------------------------------\n")
          # m.tmp = model[[i]]$weighted$fit
          m.tmp = model[[i]]$full$fit
          tmp = WoT::getParameters(m.tmp, se.fit = TRUE)
          se.problems$full[i] = any(tmp$pars$approx)

          m.tmp = model[[i]]$weighted$fit
          tmp = WoT::getParameters(m.tmp, se.fit = TRUE)
          se.problems$weighted[i] = any(tmp$pars$approx)
          # Sys.sleep(1)
        }

        which(se.problems$full)      #
        which(se.problems$weighted)  #

        # Note from above that session 21 ("6c4s02ki") is discarded

        # Model 1 has far to many states...
        m   = model[[1]]$weighted
        tmp = WoT::getParameters(m$fit, se.fit = TRUE)

        # Force model 1 to have only 2 states...
        fit.dat = subset(dots, session==fitsessions$session[1])
        full = fit_gmm(y~x0, fit.dat, opt = 'bic', trystates = 2, maxiter = 10)
        weighted = fit_gmm(y~x1, fit.dat, opt = 'bic', trystates = 2, maxiter = 10)
        model[[1]] = list(session = fitsessions$session[1], full=full, weighted=weighted)

        # Check again:
        m = model[[1]]$weighted
        tmp = WoT::getParameters(m$fit, se.fit = TRUE)
        # Now all se's are available! (except for 21 which is discarded)


        # Approx. SE when not available (NOT USED)



        for(i in 1:29){
          M.tmp = model[[i]]

          M.tmp$full = append(M.tmp$full, WoT::getParameters(M.tmp$full$fit, se.fit = TRUE))  # append parameter estimates and std.errs
          tmp2       = depmixS4::posterior(M.tmp$full$fit)                                    # append estimated states
          M.tmp$full = append(M.tmp$full,list(state=tmp2$state, sweight = tmp2[,-1]))

          M.tmp$weighted = append(M.tmp$weighted, WoT::getParameters(M.tmp$weighted$fit, se.fit = TRUE))       # append parameter estimates and std.errs
          tmp2       = depmixS4::posterior(M.tmp$weighted$fit)                                    # append estimated states
          M.tmp$weighted = append(M.tmp$weighted,list(state=tmp2$state, sweight = tmp2[,-1]))

          tmp   = list(data=dots[dots$session==model[[i]]$session,])                                             # append dataset
          M.tmp = append(M.tmp,tmp)


          if(any(M.tmp$pars$approx)){ # Make approx. confidence bounds, based on weighted regression

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


          model[[i]] = M.tmp
        }

        # i=8
        # # model[[i]]$full$pars
        # model[[i]]$weighted$pars

        fn = paste0(wdir,"data/models.Rda")
        save(model, file=fn)




###### 4. Visualize results - modelfit: residuals and qqplots #################################################################################################################################################


        # See 'Plots for paper'



###### 5. Visualize results - conclusions on beta #############################################################################################################################################################


        # See 'Plots for paper'

