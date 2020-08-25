### Code for implementing the STEMA forecasting procedure
### Author: Hannah Correia, Auburn University, 2018

library(mgcv) # GAM
library(forecast) # ARIMA
library(ggplot2)
library(reshape2)
library(spTimer) # Bayesian space-time forecasing 


###################################################### 
# snrankgam function                                 #
# Written by: Asheber Abebe, Auburn University, USA  #
#             Hannah Correia, Auburn University, USA #
######################################################

# Wilcoxon score function (sign rank)
phi.w <- function(u){
  sqrt(3)*u
}

# Bent score function (Kloke and McKean, 2014)
phi.bent <- function(u){
  ifelse(u < .75, u*4*sqrt(2)/3, sqrt(2))
}

phi.s <- function(u){
  sign(u)
}

wfun <- function(r, wtol, phi){
  rres <- phi(rank(r)/(length(r)+1))
  w <- ifelse(r > wtol, rres/r, 0)
  w/mean(w)
}

snrankgam <- function(mod, data=data, family=gaussian, phi=phi, tol=1e-3, wtol=1e-4, nmax=100){
  mod <- formula(mod)
  z0 <- gam(mod, data=data, family=family)
  y0 <- z0$fitted.values
  y <- z0$y
  r0 <- abs(y-y0)
  dff <- mean(r0)
  w <<- wfun(r0, wtol, phi)
  c <- 1
  pb <- txtProgressBar(min = 1e-10, max = nmax, initial = 0, style = 1)
  while(dff > tol & c < nmax){
    z1 <- gam(mod, data=data, weights=w, family=family)
    y1 <- z1$fitted.values
    r1 <- abs(z1$y - y1)
    w <<- wfun(r1, wtol, phi)
    dff <- sum(abs(y1 - y0))/sum(abs(y0))
    y0 <- y1
    c <- c+1
    setTxtProgressBar(pb, c)
    if(c >= nmax & dff > tol){cat("\n\n Did not converge. Consider increasing tol and/or nmax. \n\n")}
    if(dff < tol){cat("\n\n Converged. \n\n")} 
  }
  return(z1)
  close(pb)
}



#### 'stmodavg' FUNCTION TO DO MODAVG FOR EACH SPECIES ####
stmodavg <- function(dataname) {
  
  dset <- (dataname) # Load dataset
  
  species <- if(unique(dset$Species) == "Sablefish"){"sable"} else 
    if(unique(dset$Species) == "Pacific cod"){"cod"} else 
      if(unique(dset$Species) == "Pacific halibut"){"halibut"} else
        if(unique(dset$Species) == "Giant grenadier"){"grenadier"}
  
  # Set phi to be by species in rankGAMs
  phi <- if(unique(dset$Species) == "Sablefish"){phi.w} else
    if(unique(dset$Species) == "Giant grenadier"){phi.w} else
      if(unique(dset$Species) == "Pacific cod"){phi.bent} else
        if(unique(dset$Species) == "Pacific halibut"){phi.bent}
  
  # Create folder for species results & enter it
  maindir <- getwd()
  dir.create(species)
  setwd(paste0(maindir, "/", species))
  speciesdir <- getwd()
  
  # Loop to go through all possible prediction years available
  # (to replicate results in paper, use the following)
  #predyears <- c((min(dset$Year)+10):max(dset$Year))
  # Or, set specific prediction year (this will do one prediction year)
  predyears <- c((min(dset$Year)+10):(min(dset$Year)+10))
  
  # Save each type of model's predictions for each station in year=ypred
  
  for(ypred in predyears){
    
    #### COMPARED METHODS ####
    # (Uncomment if wanting to replicate original analysis)
    #### Plain time-series (ARIMA) model - no LOO ###

    stations <- unique(as.numeric(dset[which(dset$Year<=ypred),]$Station))
    years <- unique(as.numeric(dset[which(dset$Year<ypred),]$Year))

    # Number of years used for model FITTING (i.e. years before ypred); needed for leave-one-out
    ylen <- length(years)

    # Matrix to hold plain time-series prediction for each station after each leave-one-out sampling;
    # (i.e. col 1 is all station's pred. val. for year ypred after 1st leave-one-out sampling; col 2 is 2nd; etc.)
    pred.arima <- matrix(NA, max(stations)-61, 2); colnames(pred.arima) <- c("station", "pred")
    aic.arima <- matrix(NA, max(stations)-61, 2); colnames(aic.arima) <- c("station", "aic")

    # Dataset w/ all years up to and including ypred
    jack <- dset[which(dset$Year <= ypred), ]

    years.jack <- unique(as.numeric(jack$Year))

    # Have to run through each station when doing time-series
    for(i in stations){

      # Create matrix of a station's CPUE and lagged SST values
      # Have to fill in missing y with NA for CPUE
      fish.ts <- data.frame(matrix(NA, length(years)+1, 7))
      fish.ts[,1] <- c(min(years):(max(years)+1))

      fish.ts[,2] <- jack[which(jack$Station==i & jack$Year<=ypred),]$CPUE;
      fish.ts[,3] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW5;
      fish.ts[,4] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW4;
      fish.ts[,5] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW3;
      fish.ts[,6] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW2;
      fish.ts[,7] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW1

      colnames(fish.ts) <- c("year", "CPUE", "SST_cvW5", "SST_cvW4", "SST_cvW3", "SST_cvW2", "SST_cvW1")

      # SST of all years up to and including ypred
      xreg <- cbind(SST_cvW5 = fish.ts$SST_cvW5, SST_cvW4 = fish.ts$SST_cvW4, SST_cvW3 = fish.ts$SST_cvW3,
                    SST_cvW2 = fish.ts$SST_cvW2, SST_cvW1 = fish.ts$SST_cvW1)

      # CPUE of all years up to and including ypred
      # Have to fill in missing y with NA for CPUE
      CPUE.ts <- ts(fish.ts$CPUE, start = min(years), end = (max(years)+1))

      # Training data of 1990-(ypred-1)
      # SST
      xreg.train <- xreg[1:(length(years)), ]
      # CPUE
      CPUE.train <- CPUE.ts[1:(length(years))]

      # Fit model using training data
      fish_fit <- auto.arima(CPUE.train, xreg = xreg.train)

      # Predict ypred using trained model

      xreg.pred <- rbind(xreg[(length(years)+1):(length(years)+1), ], NA)

      fish_forecast <- forecast(fish_fit, h = 1, xreg = xreg.pred)

      aic.arima[i-61, 1] <- i
      aic.arima[i-61, 2] <- fish_forecast$model$aicc #Use AIC for small samples

      # Save the forecasted value in the correct row of pred.time
      # Column 2 is first leave-one-out, column 3 is second leave-one-out etc.
      pred.arima[i-61, 1] <- i
      pred.arima[i-61, 2] <- fish_forecast$mean[1]

    }

    # Remove NA values created by missing numbers in station list
    aic.arima <- na.omit(aic.arima)
    pred.arima <- na.omit(pred.arima)

    pred.name.arima <- paste("pred.arima", ypred, sep = ".")
    assign(pred.name.arima, data.frame(pred.arima), envir = .GlobalEnv)



    #### Plain time-series (ARIMA) model w/ LOO ###

    stations <- unique(as.numeric(dset[which(dset$Year<=ypred),]$Station))
    years <- unique(as.numeric(dset[which(dset$Year<ypred),]$Year))

    # Number of years used for model FITTING (i.e. years before ypred); needed for leave-one-out
    ylen <- length(years)

    # Matrix to hold plain time-series prediction for each station after each leave-one-out sampling;
    # (i.e. col 1 is all station's pred. val. for year ypred after 1st leave-one-out sampling; col 2 is 2nd; etc.)
    pred.arima <- matrix(NA, max(stations)-61, ylen + 1); colnames(pred.arima) <- c("station", paste("i =", min(years):max(years)))
    aic.arima <- matrix(NA, max(stations)-61, ylen + 1); colnames(aic.arima) <- c("station", paste("i =", min(years):max(years)))

    for(y in years){
      # Create a matrix for the n-1 years obtained from jackknife sampling;
      jack <- dset[which(dset$Year != y & dset$Year <= ypred), ]

      years.jack <- unique(as.numeric(jack$Year))

      # Have to run through each station when doing time-series
      for(i in stations){

        # Create matrix of a station's CPUE and lagged SST values
        # Have to fill in missing y with NA for CPUE
        fish.ts <- data.frame(matrix(NA, length(years)+1, 7))
        fish.ts[,1] <- c(min(years):(max(years)+1))
        if(y == min(years)){
          fish.ts[1, 2] <- NA; fish.ts[2:(length(years)+1), 2] <- jack[which(jack$Station==i & jack$Year<=ypred),]$CPUE;
          fish.ts[1, 3] <- NA; fish.ts[2:(length(years)+1), 3] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW5;
          fish.ts[1, 4] <- NA; fish.ts[2:(length(years)+1), 4] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW4;
          fish.ts[1, 5] <- NA; fish.ts[2:(length(years)+1), 5] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW3;
          fish.ts[1, 6] <- NA; fish.ts[2:(length(years)+1), 6] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW2;
          fish.ts[1, 7] <- NA; fish.ts[2:(length(years)+1), 7] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW1
        } else
          if(y > min(years) & y < (max(years)+1)){
            fish.ts[(y-min(years)+1), 2] <- NA; fish.ts[1:(y-min(years)), 2] <- jack[which(jack$Station==i & jack$Year<=ypred),]$CPUE[1:(y-min(years))]; fish.ts[(y-min(years)+2):(length(years)+1), 2] <- jack[which(jack$Station==i & jack$Year<=ypred),]$CPUE[(y-min(years)+1):(length(jack[which(jack$Station==i & jack$Year<=ypred),]$CPUE))];
            fish.ts[(y-min(years)+1), 3] <- NA; fish.ts[1:(y-min(years)), 3] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW5[1:(y-min(years))]; fish.ts[(y-min(years)+2):(length(years)+1), 3] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW5[(y-min(years)+1):(length(jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW5))];
            fish.ts[(y-min(years)+1), 4] <- NA; fish.ts[1:(y-min(years)), 4] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW4[1:(y-min(years))]; fish.ts[(y-min(years)+2):(length(years)+1), 4] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW4[(y-min(years)+1):(length(jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW4))];
            fish.ts[(y-min(years)+1), 5] <- NA; fish.ts[1:(y-min(years)), 5] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW3[1:(y-min(years))]; fish.ts[(y-min(years)+2):(length(years)+1), 5] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW3[(y-min(years)+1):(length(jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW3))];
            fish.ts[(y-min(years)+1), 6] <- NA; fish.ts[1:(y-min(years)), 6] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW2[1:(y-min(years))]; fish.ts[(y-min(years)+2):(length(years)+1), 6] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW2[(y-min(years)+1):(length(jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW2))];
            fish.ts[(y-min(years)+1), 7] <- NA; fish.ts[1:(y-min(years)), 7] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW1[1:(y-min(years))]; fish.ts[(y-min(years)+2):(length(years)+1), 7] <- jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW1[(y-min(years)+1):(length(jack[which(jack$Station==i & jack$Year<=ypred),]$SST_cvW1))]
          }
        colnames(fish.ts) <- c("year", "CPUE", "SST_cvW5", "SST_cvW4", "SST_cvW3", "SST_cvW2", "SST_cvW1")

        # SST of all years up to and including ypred
        xreg <- cbind(SST_cvW5 = fish.ts$SST_cvW5, SST_cvW4 = fish.ts$SST_cvW4, SST_cvW3 = fish.ts$SST_cvW3,
                      SST_cvW2 = fish.ts$SST_cvW2, SST_cvW1 = fish.ts$SST_cvW1)

        # CPUE of all years up to and including ypred
        # Have to fill in missing y with NA for CPUE
        CPUE.ts <- ts(fish.ts$CPUE, start = min(years), end = (max(years)+1))

        # Training data of 1990-(ypred-1)
        # SST
        xreg.train <- xreg[1:(length(years)), ]
        # CPUE
        CPUE.train <- CPUE.ts[1:(length(years))]

        # Fit model using training data
        fish_fit <- auto.arima(CPUE.train, xreg = xreg.train)

        # Predict ypred using trained model

        xreg.pred <- rbind(xreg[(length(years)+1):(length(years)+1), ], NA)

        fish_forecast <- forecast(fish_fit, h = 1, xreg = xreg.pred)

        aic.arima[i-61, 1] <- i
        aic.arima[i-61, (y-min(years)+2)] <- fish_forecast$model$aicc #Use AIC for small samples

        # Save the forecasted value in the correct row of pred.time
        # Column 2 is first leave-one-out, column 3 is second leave-one-out etc.
        pred.arima[i-61, 1] <- i
        pred.arima[i-61, (y-min(years)+2)] <- fish_forecast$mean[1]

        # NOTE: six warnings produced for sablefish for one ypred run
        # In forecast.Arima(fish_fit, h = 1, xreg = xreg.pred) :
        # Upper prediction intervals are not finite.

      }
    }

    # Remove NA values created by missing numbers in station list
    aic.arima <- na.omit(aic.arima)
    pred.arima <- na.omit(pred.arima)

    # Compute leave-one out values
    pred.est.arima_l1o <- cbind(pred.arima[,1], rowMeans(pred.arima[,2:ncol(pred.arima)]), NA, NA) # Mean
    n <- length(years) #number of leave-one-out estimates created for each station (removed years 1990:1999 one at a time per station => 10 est)
    # Variance of leave-one-out predictions
    for(i in 1:nrow(pred.arima)){pred.est.arima_l1o[i,3] <- ((n-1)^2/n)*var(pred.arima[i, 2:ncol(pred.arima)]) + 1}
    # Std error of leave-one-out predictions
    for(i in 1:nrow(pred.arima)){pred.est.arima_l1o[i,4] <- sqrt(((n-1)^2/n)*var(pred.arima[i, 2:ncol(pred.arima)]) + 1)}

    colnames(pred.est.arima_l1o) <- c("station", "l1o.mean", "l1o.var", "l1o.sterr")

    pred.name.arima <- paste("pred.arima_l1o", ypred, sep = ".")
    assign(pred.name.arima, data.frame(pred.est.arima_l1o), envir = .GlobalEnv)



    #### Naive spatial GAM w/ LOO - no averaging over time performed ###

    stations <- unique(as.numeric(dset[which(dset$Year<=ypred),]$Station))
    years <- unique(as.numeric(dset[which(dset$Year<=ypred),]$Year))

    # Matrix to hold the predicted values for each station after each leave-one-out sampling
    # (i.e. col 1 is all station's pred. val. for year ypred after 1st leave-one-out sampling; col 2 is 2nd; etc.)
    pred.spGAM <- matrix(NA, length(stations), length(stations)+1)
    colnames(pred.spGAM) <- c("station_removed", stations)
    row.names(pred.spGAM) <- stations

    for(s in stations){
      # Create a matrix for the n-1 stations obtained from leave-one-out sampling;
      looset <- dset[which(dset$Station != s), ]

      stations.l1o <- unique(as.numeric(looset$Station))

      # Create a subset of leave-one-out sample with only year we wish to fit & fit GAM to it
      mod <- snrankgam(CPUE ~ s(Longitude, Latitude, k=5) + s(Longitude, Latitude, by = SST_cvW5, k=5) + s(Longitude, Latitude, by = SST_cvW4, k=5) +
                         s(Longitude, Latitude, by = SST_cvW3, k=5) + s(Longitude, Latitude, by = SST_cvW2, k=5) + s(Longitude, Latitude, by = SST_cvW1, k=5),
                       data = looset[which(looset$Year==ypred-1),], family=gaussian, phi = phi)

      # Next year's prediction is previous year's fitted value
      sc <- as.character(s)
      sc.l1o <- as.character(stations.l1o)

      pred.spGAM[row.names(pred.spGAM) %in% sc, 1] <- s
      pred.spGAM[row.names(pred.spGAM) %in% sc, colnames(pred.spGAM) %in% sc] <- NA
      pred.spGAM[row.names(pred.spGAM) %in% sc, colnames(pred.spGAM) %in% sc.l1o] <- mod$fitted.values

    }


    # Spatial leave-one-out calculations: Note that the spatial predictions should be averaged by column to obtain the
    # mean leave-one-out prediction for the station whose number is at the top of the column

    # Compute leave-one-out values
    pred.est.spGAM <- cbind(as.numeric(row.names(pred.spGAM)), colMeans(pred.spGAM[,2:ncol(pred.spGAM)], na.rm = TRUE), NA, NA) # Mean
    n <- length(stations) #number of leave-one-out values created for each station (removed a station one at a time per station => 73 est)
    # Variance of leave-one-out predictions
    for(i in 1:(ncol(pred.spGAM)-1)){pred.est.spGAM[i,3] <- ((n-1)^2/n)*var(pred.spGAM[, i+1], na.rm = TRUE) + 1}
    # Std error of leave-one-out predictions
    for(i in 1:(ncol(pred.spGAM)-1)){pred.est.spGAM[i,4] <- sqrt(((n-1)^2/n)*var(pred.spGAM[, i+1], na.rm = TRUE) + 1)}

    colnames(pred.est.spGAM) <- c("station", "l1o.mean", "l1o.var", "l1o.stderr")

    pred.name.spGAM <- paste("pred.spGAM_l1o", ypred, sep = ".")
    assign(pred.name.spGAM, data.frame(pred.est.spGAM), envir = .GlobalEnv)




    #### Bayesian spTimer no LOO ###

    stations <- unique(as.numeric(dset[which(dset$Year<=ypred),]$Station))
    years <- unique(as.numeric(dset[which(dset$Year<ypred),]$Year))

    # Number of years used for model FITTING (i.e. years before ypred); needed for leave-one-out
    ylen <- length(years)

    # Matrix to hold plain time-series prediction for each station after each leave-one-out sampling;
    # (i.e. col 1 is all station's pred. val. for year ypred after 1st leave-one-out sampling; col 2 is 2nd; etc.)
    pred.bayes <- matrix(NA, length(stations), 2); colnames(pred.bayes) <- c("station", "pred")
    pred.bayes[, 1] <- stations

    # Training data of 1990-(ypred-1)
    data.fit <- dset[which(dset$Year<ypred),]

    # Set ideal acceptance rate to 32%; search for fit (using training data) that obtains closest acceptance rate
    opt.acc <- 32
    prev.diff <- 100
    k <- 1
    gp.gamm <- spT.Gibbs(formula = CPUE ~ SST_cvW5 + SST_cvW4 + SST_cvW3 + SST_cvW2 + SST_cvW1,
                         data = data.fit, model = "GP", coords = ~ Longitude + Latitude, scale.transform = "NONE",
                         spatial.decay = spT.decay(distribution=Gamm(2,1), tuning = 10^k))

    acc.diff <- abs(opt.acc - gp.gamm$accept)
    while(acc.diff <= prev.diff){
      prev.diff <- acc.diff
      k <- k - 1
      gp.old <- gp.gamm
      gp.gamm <- spT.Gibbs(formula = CPUE ~ SST_cvW5 + SST_cvW4 + SST_cvW3 + SST_cvW2 + SST_cvW1,
                           data = data.fit, model = "GP", coords = ~ Longitude + Latitude, scale.transform = "NONE",
                           spatial.decay = spT.decay(distribution=Gamm(2,1), tuning = 10^k))
      acc.diff <- abs(opt.acc - gp.gamm$accept)
    }
    # Save appropriate fit model
    bayes_fit <- gp.old

    # Predict ypred using trained model
    data.pred <- dset[which(dset$Year==ypred),]

    bayes_forecast <- predict(bayes_fit, newdata = data.pred,
                              newcoords = ~ Longitude + Latitude,
                              type = "temporal", foreStep = 1)

    # Save the forecasted value in the correct row of pred.time
    # Column 2 is first leave-one-out, column 3 is second leave-one-out etc.
    pred.bayes[, 2] <- bayes_forecast$Mean

    pred.est.bayes <- pred.bayes
    pred.name.bayes <- paste("pred.bayes", ypred, sep = ".")
    assign(pred.name.bayes, data.frame(pred.est.bayes), envir = .GlobalEnv)



    #### Bayesian spTimer w/ LOO ###

    stations <- unique(as.numeric(dset[which(dset$Year<=ypred),]$Station))
    years <- unique(as.numeric(dset[which(dset$Year<ypred),]$Year))

    # Number of years used for model FITTING (i.e. years before ypred); needed for leave-one-out
    ylen <- length(years)

    # Matrix to hold plain time-series prediction for each station after each leave-one-out sampling;
    # (i.e. col 1 is all station's pred. val. for year ypred after 1st leave-one-out sampling; col 2 is 2nd; etc.)
    pred.bayes <- matrix(NA, length(stations), ylen + 1); colnames(pred.bayes) <- c("station", paste("i =", min(years):max(years)))
    pred.bayes[, 1] <- stations

    for(y in years){
      # Create a matrix for the n-1 years obtained from jackknife sampling;
      jack <- dset[which(dset$Year != y & dset$Year <= ypred), ]

      years.jack <- unique(as.numeric(jack$Year))

      # Does not have to run through each station when doing Bayes spt forecasting

      # Training data of 1990-(ypred-1)
      data.fit <- jack[which(jack$Year<ypred),]

      # Set ideal acceptance rate to 32%; search for fit (using training data) that obtains closest acceptance rate
      opt.acc <- 32
      prev.diff <- 100
      k <- 1
      gp.gamm <- spT.Gibbs(formula = CPUE ~ SST_cvW5 + SST_cvW4 + SST_cvW3 + SST_cvW2 + SST_cvW1,
                           data = data.fit, model = "GP", coords = ~ Longitude + Latitude, scale.transform = "NONE",
                           spatial.decay = spT.decay(distribution=Gamm(2,1), tuning = 10^k))

      acc.diff <- abs(opt.acc - gp.gamm$accept)
      while(acc.diff <= prev.diff){
        prev.diff <- acc.diff
        k <- k - 1
        gp.old <- gp.gamm
        gp.gamm <- spT.Gibbs(formula = CPUE ~ SST_cvW5 + SST_cvW4 + SST_cvW3 + SST_cvW2 + SST_cvW1,
                             data = data.fit, model = "GP", coords = ~ Longitude + Latitude, scale.transform = "NONE",
                             spatial.decay = spT.decay(distribution=Gamm(2,1), tuning = 10^k))
        acc.diff <- abs(opt.acc - gp.gamm$accept)
      }
      # Save appropriate fit model
      bayes_fit <- gp.old

      # Predict ypred using trained model
      data.pred <- jack[which(jack$Year==ypred),]

      bayes_forecast <- predict(bayes_fit, newdata = data.pred,
                                newcoords = ~ Longitude + Latitude,
                                type = "temporal", foreStep = 1)

      # Save the forecasted value in the correct row of pred.time
      # Column 2 is first leave-one-out, column 3 is second leave-one-out etc.
      pred.bayes[, (y-min(years)+2)] <- bayes_forecast$Mean

    }

    # Compute leave-one out values
    pred.est.bayes_l1o <- cbind(pred.bayes[,1], rowMeans(pred.bayes[,2:ncol(pred.bayes)]), NA, NA) # Mean
    n <- length(years) #number of leave-one-out estimates created for each station (removed years 1990:1999 one at a time per station => 10 est)
    # Variance of leave-one-out predictions
    for(i in 1:nrow(pred.bayes)){pred.est.bayes_l1o[i,3] <- ((n-1)^2/n)*var(pred.bayes[i, 2:ncol(pred.bayes)]) + 1}
    # Std error of leave-one-out predictions
    for(i in 1:nrow(pred.bayes)){pred.est.bayes_l1o[i,4] <- sqrt(((n-1)^2/n)*var(pred.bayes[i, 2:ncol(pred.bayes)]) + 1)}

    colnames(pred.est.bayes_l1o) <- c("station", "l1o.mean", "l1o.var", "l1o.sterr")

    pred.name.bayes <- paste("pred.bayes_l1o", ypred, sep = ".")
    assign(pred.name.bayes, data.frame(pred.est.bayes_l1o), envir = .GlobalEnv)



    #### Step (1) of STEMA technique - Temporal varying coefficient model with ARIMA w/ LOO ####

    # Get station list and year list for each ypred
    stations <- unique(as.numeric(dset[which(dset$Year<=ypred),]$Station))
    years <- unique(as.numeric(dset[which(dset$Year<ypred),]$Year))

    # Number of years used for model FITTING (i.e. years before ypred); needed for leave-one-out
    ylen <- length(years)

    # Matrix to hold the predicted values for each station after each leave-one-out sampling
    # (i.e. col 1 is all station's pred. val. for year ypred after 1st leave-one-out sampling; col 2 is 2nd; etc.)
    pred.t <- matrix(NA, max(stations)-61, ylen + 1); colnames(pred.t) <- c("station", paste("i =", min(years):max(years)))
    aic.t <- matrix(NA, max(stations)-61, ylen + 1); colnames(aic.t) <- c("station", paste("i =", min(years):max(years)))

    for(y in years){
      # Create a matrix for the n-1 years for leave-one-out process
      jack <- dset[which(dset$Year != y & dset$Year < ypred), ]

      years.jack <- unique(as.numeric(jack$Year))

      for(i in stations){
        # Fit a varying coefficient model to obtain CPUE fitted values
        mod <- snrankgam(CPUE ~ s(Year, k=3) + s(Year, by = SST_cvW5, k=3),
                         data = jack[which(jack$Station==i & jack$Year<ypred), ],
                         family=gaussian, phi = phi)

        # Create a time-series object from the fitted values of the varying coefficient model
        # thereby gaining some of the smoothing & information of SST
        ts.gam <- data.frame(matrix(NA, length(years), 2))
        ts.gam[, 1] <- years
        if(y == min(years)){ts.gam[1, 2] <- NA; ts.gam[2:(length(years)), 2] <- mod$fitted.values
        } else if(y > min(years) & y < max(years)){ts.gam[(y-min(years)+1), 2] <- NA;
        ts.gam[1:(y-min(years)), 2] <- mod$fitted.values[1:(y-min(years))]; ts.gam[(y-min(years)+2):(length(years)), 2] <- mod$fitted.values[(y-min(years)+1):(length(mod$fitted.values))]
        } else if(y == max(years)){ts.gam[1:(length(years)-1), 2] <- mod$fitted.values; ts.gam[length(years), 2] <- NA}
        # Training data of 1990-(ypred-1)
        train <- ts(ts.gam[,2], start = min(years), end = max(years))
        # Testing data of ypred
        test <- ts(dset[which(dset$Station==i & dset$Year==ypred),]$CPUE, start = ypred)

        # fit 1990-(ypred-1) data using auto.arima (chooses best AR model)
        auto_fit <- auto.arima(train)
        # predict CPUE for ypred using forecast();
        auto_forecast <- forecast(auto_fit, h = 1)
        # accuracy of the time-series model in predicting future years
        #auto_accuracy <- accuracy(auto_forecast, test) # causes error when y=(ypred-1)

        aic.t[i-61, 1] <- i
        aic.t[i-61, (y-min(years)+2)] <- auto_forecast$model$aicc #Use AIC for small samples

        pred.t[i-61, 1] <- i
        pred.t[i-61, (y-min(years)+2)] <- auto_forecast$mean[1]
      }
    }

    # Remove NA values created by missing numbers in station list
    aic.t <- na.omit(aic.t)
    pred.t <- na.omit(pred.t)

    # Compute leave-one out values
    pred.est.t <- cbind(pred.t[,1], rowMeans(pred.t[,2:ncol(pred.t)]), NA, NA) # Mean
    n <- length(years) #number of leave-one-out estimates created for each station (removed years 1990:1999 one at a time per station => 10 est)
    # Variance of leave-one-out predictions
    for(i in 1:nrow(pred.t)){pred.est.t[i,3] <- ((n-1)^2/n)*var(pred.t[i, 2:ncol(pred.t)]) + 1}
    # Std error of leave-one-out predictions
    for(i in 1:nrow(pred.t)){pred.est.t[i,4] <- sqrt(((n-1)^2/n)*var(pred.t[i, 2:ncol(pred.t)]) + 1)}

    colnames(pred.est.t) <- c("station", "l1o.mean", "l1o.var", "l1o.sterr")

    pred.name.t <- paste("pred.t", ypred, sep = ".")
    assign(pred.name.t, data.frame(pred.est.t), envir = .GlobalEnv)


    #### Step (2) of STEMA technique - Spatial varying coefficient model with ARIMA w/ LOO ####

    # Get station list and year list for each ypred
    stations <- unique(as.numeric(dset[which(dset$Year<=ypred),]$Station))
    years <- unique(as.numeric(dset[which(dset$Year<ypred),]$Year))

    # Number of years used for model FITTING (i.e. years before ypred); needed for leave-one-out
    ylen <- length(years)

    # Matrix to hold the predicted values for each station after each leave-one-out sampling
    # (i.e. col 1 is all station's pred. val. for year ypred after 1st leave-one-out sampling; col 2 is 2nd; etc.)
    pred.spl1o <- matrix(NA, length(stations), length(stations)+1)
    colnames(pred.spl1o) <- c("station_removed", stations)
    row.names(pred.spl1o) <- stations

    for(s in stations){
      # Create a matrix for the n-1 stations obtained from leave-one-out sampling;
      looset <- dset[which(dset$Station != s), ]

      stations.l1o <- unique(as.numeric(looset$Station))

      # Need a matrix to hold predictions for n-1 stations for current year i ready for averaging
      pred.gam <- matrix(NA, length(stations.l1o), length(years)); rownames(pred.gam) <- stations.l1o; colnames(pred.gam) <- years

      for(i in years){
        # Create a subset of leave-one-out sample with only year we wish to fit & fit varying coefficient model to it
        mod <- snrankgam(CPUE ~ s(Longitude, Latitude, k=5) + s(Longitude, Latitude, by = SST_cvW5, k=5) +
                           s(Longitude, Latitude, by = SST_cvW4, k=5) +
                           s(Longitude, Latitude, by = SST_cvW3, k=5) +
                           s(Longitude, Latitude, by = SST_cvW2, k=5) +
                           s(Longitude, Latitude, by = SST_cvW1, k=5),
                         data = looset[which(looset$Year==i),], family=gaussian, phi = phi)

        # Save fitted varying coefficient model values into matrix column corresponding to the year of fitting; 1990 is column 1
        pred.gam[, (i-min(years)+1)] <- mod$fitted.values
      }

      # Create matrix to hold predicted values of ARIMA for each station fit for ypred in the leave-one-out step
      pred.sp <- matrix(NA, max(stations.l1o)-61, 2); colnames(pred.sp) <- c("station", "pred.l1o")
      aic.sp <- matrix(NA, max(stations.l1o)-61, 2); colnames(aic.sp) <- c("station", "pred.l1o")

      for(p in stations.l1o){
        # Create a time-series object from the fitted values of the varying coefficient model for each station
        sp.gam <- data.frame(matrix(NA, length(years), 2))
        sp.gam[, 1] <- years
        sp.gam[, 2] <- pred.gam[row.names(pred.gam) %in% p, ]

        # Training data of 1990-(ypred-1)
        train <- ts(sp.gam[,2], start = min(years), end = max(years))
        # Testing data of ypred
        test <- ts(dset[which(dset$Station==p & dset$Year==ypred),]$CPUE, start = ypred)

        # fit 1990-(ypred-1) data using auto.arima (chooses best AR model)
        auto_fit <- auto.arima(train)
        # predict CPUE for ypred using forecast();
        auto_forecast <- forecast(auto_fit, h = 1)
        # accuracy of the time-series model in predicting future years
        #auto_accuracy <- accuracy(auto_forecast, test) # causes error when y=(ypred-1)

        aic.sp[p-61, 1] <- p
        aic.sp[p-61, 2] <- auto_forecast$model$aicc #Use AIC for small samples

        pred.sp[p-61, 1] <- p
        pred.sp[p-61, 2] <- auto_forecast$mean[1]

      }

      # After all stations have been fit, need to save l1o predictions for that l1o ypred into corresponding column of matrix
      # Match each row of 'spl1o' with the corresponding column in 'pred.spl1o' with the correct station's number attached
      sc <- as.character(s)
      sc.l1o <- as.character(stations.l1o)

      pred.spl1o[row.names(pred.spl1o) %in% sc, 1] <- s
      pred.spl1o[row.names(pred.spl1o) %in% sc, colnames(pred.spl1o) %in% sc] <- NA
      pred.spl1o[row.names(pred.spl1o) %in% sc, colnames(pred.spl1o) %in% sc.l1o] <- pred.sp[pred.sp[, 1] %in% sc.l1o, 2]

    }


    # Spatial leave-one-out calculations: Note that the spatial predictions should be averaged by column to obtain the
    # mean leave-one-out prediction for the station whose number is at the top of the column

    # Compute leave-one-out values
    pred.est.sp <- cbind(as.numeric(row.names(pred.spl1o)), colMeans(pred.spl1o[,2:ncol(pred.spl1o)], na.rm = TRUE), NA, NA) # Mean
    n <- length(stations) #number of leave-one-out values created for each station (removed a station one at a time per station => 73 est)
    # Variance of leave-one-out predictions
    for(i in 1:(ncol(pred.spl1o)-1)){pred.est.sp[i,3] <- ((n-1)^2/n)*var(pred.spl1o[, i+1], na.rm = TRUE) + 1}
    # Std error of leave-one-out predictions
    for(i in 1:(ncol(pred.spl1o)-1)){pred.est.sp[i,4] <- sqrt(((n-1)^2/n)*var(pred.spl1o[, i+1], na.rm = TRUE) + 1)}

    colnames(pred.est.sp) <- c("station", "l1o.mean", "l1o.var", "l1o.stderr")

    pred.name.sp <- paste("pred.sp", ypred, sep = ".")
    assign(pred.name.sp, data.frame(pred.est.sp), envir = .GlobalEnv)



    #### Step (3) of STEMA technique - mod averaging ####
    #### Weigh averaged space pred and time pred for ypred

    weighted.pred <- data.frame(matrix(NA, length(stations), 8))
    colnames(weighted.pred) <- c("station", "l1o.sp.mean", "l1o.sp.err", "l1o.t.mean", "l1o.t.err", "sp.weight", "t.weight", "avg.pred")

    weighted.pred[,1] <- pred.est.sp[,1] # station
    weighted.pred[,2] <- pred.est.sp[,2] # l1o.sp.mean
    weighted.pred[,3] <- pred.est.sp[,4] # l1o.sp.err
    weighted.pred[,4] <- pred.est.t[,2] # l1o.t.mean
    weighted.pred[,5] <- pred.est.t[,4] # l1o.t.err

    weighted.pred[,6] <- weighted.pred$l1o.t.err/(weighted.pred$l1o.sp.err + weighted.pred$l1o.t.err) # Spatial weight

    weighted.pred[,7] <- weighted.pred$l1o.sp.err/(weighted.pred$l1o.sp.err + weighted.pred$l1o.t.err) # Temporal weight

    # Final averaged SPT prediction
    weighted.pred[,8] <- (weighted.pred$l1o.sp.mean * weighted.pred$sp.weight) + (weighted.pred$l1o.t.mean * weighted.pred$t.weight)

    weighted.name <- paste("weighted", ypred, sep = ".")
    assign(weighted.name, weighted.pred, envir = .GlobalEnv)

  }
  
  
  #### Save all predictions (or mean predictions) in tables by method ####
  # Put together the ARIMA (no LOO) predictions in a table by year
  myList <- lapply( paste0("pred.arima.", predyears) , get)
  myList2 <- lapply(myList, function(x) x[2])
  myData <- do.call(cbind, myList2 )
  pred.arima.all <- cbind(get(paste0("pred.arima.", predyears[1]))$station, myData)
  
  colnames(pred.arima.all) <- c("station", predyears)
  
  setwd(speciesdir)
  save(pred.arima.all, file = "pred.arima.all.rda")
  write.table(pred.arima.all, file = "pred.arima.all.txt", sep = "\t")
  
  
  
  # Put together the ARIMA (w/ LOO) predictions in a table by year
  myList <- lapply( paste0("pred.arima_l1o.", predyears) , get)
  myList2 <- lapply(myList, function(x) x[2])
  myData <- do.call(cbind, myList2 )
  pred.arima_l1o.all <- cbind(get(paste0("pred.arima_l1o.", predyears[1]))$station, myData)

  colnames(pred.arima_l1o.all) <- c("station", predyears)

  setwd(speciesdir)
  save(pred.arima_l1o.all, file = "pred.arima_l1o.all.rda")
  write.table(pred.arima_l1o.all, file = "pred.arima_l1o.all.txt", sep = "\t")



  # Put together the spatial GAM predictions in a table by year
  myList <- lapply( paste0("pred.spGAM_l1o.", predyears) , get)
  myList2 <- lapply(myList, function(x) x[2])
  myData <- do.call(cbind, myList2 )
  pred.spGAM_l1o.all <- cbind(get(paste0("pred.spGAM_l1o.", predyears[1]))$station, myData)

  colnames(pred.spGAM_l1o.all) <- c("station", predyears)

  setwd(speciesdir)
  save(pred.spGAM_l1o.all, file = "pred.spGAM_l1o.all.rda")
  write.table(pred.spGAM_l1o.all, file = "pred.spGAM_l1o.all.txt", sep = "\t")



  # Put together the Bayes spTimer (no LOO) predictions in a table by year
  myList <- lapply( paste0("pred.bayes.", predyears) , get)
  myList2 <- lapply(myList, function(x) x[2])
  myData <- do.call(cbind, myList2 )
  pred.bayes.all <- cbind(get(paste0("pred.bayes.", predyears[1]))$station, myData)

  colnames(pred.bayes.all) <- c("station", predyears)

  setwd(speciesdir)
  save(pred.bayes.all, file = "pred.bayes.all.rda")
  write.table(pred.bayes.all, file = "pred.bayes.all.txt", sep = "\t")



  # Put together the Bayes spTimer (w/ LOO) predictions in a table by year
  myList <- lapply( paste0("pred.bayes_l1o.", predyears) , get)
  myList2 <- lapply(myList, function(x) x[2])
  myData <- do.call(cbind, myList2 )
  pred.bayes_l1o.all <- cbind(get(paste0("pred.bayes_l1o.", predyears[1]))$station, myData)

  colnames(pred.bayes_l1o.all) <- c("station", predyears)

  setwd(speciesdir)
  save(pred.bayes_l1o.all, file = "pred.bayes_l1o.all.rda")
  write.table(pred.bayes_l1o.all, file = "pred.bayes_l1o.all.txt", sep = "\t")



  # Put together the temporal GAM predictions in a table by year
  myList <- lapply( paste0("pred.t.", predyears) , get)
  myList2 <- lapply(myList, function(x) x[2])
  myData <- do.call(cbind, myList2 )
  pred.t.all <- cbind(get(paste0("pred.t.", predyears[1]))$station, myData)

  colnames(pred.t.all) <- c("station", predyears)

  setwd(speciesdir)
  save(pred.t.all, file = "pred.t.all.rda")
  write.table(pred.t.all, file = "pred.t.all.txt", sep = "\t")



  # Put together the averaged spatial predictions in a table by year
  myList <- lapply( paste0("pred.sp.", predyears) , get)
  myList2 <- lapply(myList, function(x) x[2])
  myData <- do.call(cbind, myList2 )
  pred.sp.all <- cbind(get(paste0("pred.sp.", predyears[1]))$station, myData)

  colnames(pred.sp.all) <- c("station", predyears)

  setwd(speciesdir)
  save(pred.sp.all, file = "pred.sp.all.rda")
  write.table(pred.sp.all, file = "pred.sp.all.txt", sep = "\t")



  # Put together the averaged weighted predictions in a table by year
  myList <- lapply( paste0("weighted.", predyears) , get)
  myList2 <- lapply(myList, function(x) x[8])
  myData <- do.call(cbind, myList2 )
  weighted.pred.all <- cbind(get(paste0("weighted.", predyears[1]))$station, myData)

  colnames(weighted.pred.all) <- c("station", predyears)

  setwd(speciesdir)
  save(weighted.pred.all, file = "weighted.pred.all.rda")
  write.table(weighted.pred.all, file = "weighted.pred.all.txt", sep = "\t")


  # #### Plot each station's CPUE predictions over the ypred years (2000-2012) ####
  # # Each color represents a different method
  # # Faded color lines are individual stations
  # # Bold color lines are mean predictions over all stations
  # 
  # # First create huge dataset with years, stations, pred method, & pred value
  # pred.arima.all.zz <- melt(pred.arima.all, id.vars = c("station")); colnames(pred.arima.all.zz) <- c("station", "year", "pred.CPUE")
  # pred.arima.all.zz$pred.method <- "ARIMA"
  # pred.arima_l1o.all.zz <- melt(pred.arima_l1o.all, id.vars = c("station")); colnames(pred.arima_l1o.all.zz) <- c("station", "year", "pred.CPUE")
  # pred.arima_l1o.all.zz$pred.method <- "ARIMA_l1o"
  # pred.spGAM_l1o.all.zz <- melt(pred.spGAM_l1o.all, id.vars = c("station")); colnames(pred.spGAM_l1o.all.zz) <- c("station", "year", "pred.CPUE")
  # pred.spGAM_l1o.all.zz$pred.method <- "spGAM_l1o"
  # pred.bayes.all.zz <- melt(pred.bayes.all, id.vars = c("station")); colnames(pred.bayes.all.zz) <- c("station", "year", "pred.CPUE")
  # pred.bayes.all.zz$pred.method <- "Bayes"
  # pred.bayes_l1o.all.zz <- melt(pred.bayes_l1o.all, id.vars = c("station")); colnames(pred.bayes_l1o.all.zz) <- c("station", "year", "pred.CPUE")
  # pred.bayes_l1o.all.zz$pred.method <- "Bayes_l1o"
  # pred.t.all.zz <- melt(pred.t.all, id.vars = c("station")); colnames(pred.t.all.zz) <- c("station", "year", "pred.CPUE")
  # pred.t.all.zz$pred.method <- "tGAM-arima"
  # pred.sp.all.zz <- melt(pred.sp.all, id.vars = c("station")); colnames(pred.sp.all.zz) <- c("station", "year", "pred.CPUE")
  # pred.sp.all.zz$pred.method <- "spGAM-arima"
  # weighted.pred.all.zz <- melt(weighted.pred.all, id.vars = c("station")); colnames(weighted.pred.all.zz) <- c("station", "year", "pred.CPUE")
  # weighted.pred.all.zz$pred.method <- "STEMA"
  # true.zz <- cbind(subset(dset, Year %in% predyears)$Station, subset(dset, Year %in% predyears)$Year, subset(dset, Year %in% predyears)$CPUE, "true"); colnames(true.zz) <- c("station", "year", "pred.CPUE", "pred.method")
  # 
  # zz <- rbind(pred.arima.all.zz, pred.arima_l1o.all.zz, pred.spGAM_l1o.all.zz, pred.bayes.all.zz, pred.bayes_l1o.all.zz, pred.t.all.zz, pred.sp.all.zz, weighted.pred.all.zz, true.zz)
  # zz$station <- as.numeric(zz$station); zz$pred.CPUE <- as.numeric(zz$pred.CPUE); zz$pred.method <- as.factor(zz$pred.method)
  # 
  # # Save huge dataset of predictions by method
  # setwd(speciesdir)
  # save(zz, file = "predictions.SST_final.rda")
  # 
  # # Finally, plot the CPUEs by station with separate graphs for each method
  # pdf(file = paste0(species, "-pred.pdf")) # Without loess smooth
  # print( ggplot(zz, aes(year, pred.CPUE, color = pred.method, group = interaction(pred.method, station))) + theme_bw() + geom_line() + lims(y=range(zz$pred.CPUE)) +
  #          labs(title = "Predicted CPUE by Method", x="Year", y="CPUE", color="") +  facet_wrap(~pred.method) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) )
  # dev.off()
  
  
  
  #### Obtain errors for each station, each year, each method ####
  #### Forecast errors are e=(true-predicted); 
  
  # Change back to main directory in preparation for next species
  setwd(maindir)
  
}


#### Load final datasets ####

dset.env <- new.env()

setwd(getwd())

load("stema_data.rda")

# Split data into 4 subsets within new environment:
# Load only one species into new environment at a time for shorter computation time.
dset.env$sable.stema <- stema_data[stema_data$Species=="Sablefish",]
dset.env$cod.stema <- stema_data[stema_data$Species=="Pacific cod",]
dset.env$halibut.stema <- stema_data[stema_data$Species=="Pacific halibut",]
dset.env$grenadier.stema <- stema_data[stema_data$Species=="Giant grenadier",]


# Do stmodavg() procedure with all prediction methods 
# for all ypreds for each species.
# NOTE: This operation takes significant time.
# (Approximately 45 minutes per species per prediction year 
# on a computer with a 3.1 GHz Intel Core i7 processor.)
eapply(dset.env, stmodavg)


