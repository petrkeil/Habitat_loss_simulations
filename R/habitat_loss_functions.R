
require(spatstat)
require(plyr)

# ------------------------------------------------------------------------------
# FUNCTION FITTING THE THOMAS MODEL TO XY POINT DATA

# Arguments: x and y, the coordinates of the spatial position of the individuals.

# The function returns parameters that can be used in the `rThomas` function.
  Thomas <- function(x, y) # x and y coordinates of the points
  {
    if(length(x)==1) # in case there is just a single point
    {
      return(c(kappa=1/(1000*500), sigma=1, mu=1))  
    }
    
    else # fit the Thomas model and return the parameters
    {
      W <- as.owin(list(xrange=c(0,1000), yrange=c(0,500)))
      xy <- as.ppp(cbind(x, y), W)
      fit.xy <- clusterfit(xy, clusters="Thomas")
      return(fit.xy$modelpar)
    }
  }

# ------------------------------------------------------------------------------
# FUNCTION SIMULATING THE THOMAS PATTERN FOR A SET OF SPECIES

# Arguments:
# param.data - data.frame with columns: "Latin", "kappa", "sigma" and "mu"
# indcs - indices of the rows in the param.data that shuld be included
  
  sim.Thomas.com <- function(param.data, indcs = 1:nrow(param.data))
  {
    simcom.list <- list() # empty list to store the species
    W2000 <- as.owin(list(xrange=c(0,2000), yrange=c(0,2000))) # bounding window
    
    for(i in indcs)
    {
      sp.name <- as.character(param.data[i, "Latin"])
      # print(sp.name)
      # drawing from the Thomas model using the empirical parameters
      sp.pred <- rThomas(kappa = param.data[i, "kappa"], 
                         scale = param.data[i, "sigma"], 
                         mu = param.data[i, "mu"], 
                         win = W2000)
      sp.pred.data <- data.frame(x = sp.pred$x, 
                                 y = sp.pred$y,
                                 marks = rep(sp.name, times=npoints(sp.pred)))
      simcom.list[[i]] <- sp.pred.data
      names(simcom.list)[i] <- sp.name
    }
    
    return(simcom.list)
  }
  
  # Example:
  # a <- sim.Thomas.com(Thomas.res, 1:4)


# ------------------------------------------------------------------------------
# FUNCTION CONVERTING THE PRODUCED LIST TO A ONE BIG MULTISPECIES ppp OBJECT
   
  merge.simcom.list <- function(simcom.list)
  {
    simcom.data <- ldply(simcom.list)
    
    W2000 <- as.owin(list(xrange=c(0,2000), yrange=c(0,2000))) # bounding window
 
    ppp.multisp <- ppp(x = simcom.data$x, 
                       y = simcom.data$y, 
                       marks = simcom.data$marks,
                       window = W2000)
    return(ppp.multisp)
  }
  
  # Example:
  # big.ppp <- merge.simcom.list(a)
  
# ------------------------------------------------------------------------------
# FIT AND EXTRAPOLATE THE THOMAS MODEL  
# Function that fits, extrapolates (to a 2x2 km area) and plots the 
# Thomas model for a specific species:
  
# Arguments: 
# spec.name - character, name of the species to be used
# bci - data.frame, with species name in the first column, and x-y coordinates
#       in the second column.

  Thomas.plot <- function(spec.name, bci)
  {
    sp <- bci[bci$Latin==spec.name,]
    sp.fit <- Thomas(sp$PX, sp$PY)
    W2000 <- as.owin(list(xrange=c(0,2000), yrange=c(0,2000)))
      
    sp.pred <- rThomas(sp.fit['kappa'], 
                       sp.fit['sigma'], 
                       sp.fit['mu'], 
                       win=W2000)
    
    plot(sp.pred, pch=19, cex=0.3, cols="black", main=spec.name)
    rect(0,0,1000,500, col="white")
    points(sp$PX, sp$PY, col="red", pch=19, cex=0.3)
  }

# ------------------------------------------------------------------------------
# READ THE HABITAT LOSS RASTER
# and return it as a matrix

  read.loss.matrix <- function(filename)
  {
    loss <- as.matrix(read.csv(filename, header=FALSE))
    storage.mode(loss) <- "numeric"
    return(loss)
  }

# Example:    
 a <- read.loss.matrix("../data/loss_matrices/branching.csv")
  
# ------------------------------------------------------------------------------
# READ THE HABITAT LOSS RASTER
# and return it as spatstat's owin object

  read.loss.owin <- function(filename)
  {
    loss <- as.matrix(read.csv(filename, header=FALSE))
    storage.mode(loss) <- "numeric"
    loss[loss==0] <- NA
    loss <- as.owin(as.im(loss, W=list(xrange=c(0,2000), yrange=c(0,2000))))
    return(loss)
  }
 
 # Example:    
 #  a <- read.loss.owin("../data/loss_matrices/branching.csv")
 # W2000 <- as.owin(list(xrange=c(0,2000), yrange=c(0,2000)))
 # sp.pred <- rThomas(0.00001, 15, 13, win=W2000)
 # plot(sp.pred[a])

# ------------------------------------------------------------------------------
# READ THE REMAINING HABITAT RASTER
# and return it as spatstat's owin object

  read.remain.owin <- function(filename)
  {
    remaining <- as.matrix(read.csv(filename, header=FALSE))
    storage.mode(remaining) <- "numeric"
    remaining[remaining==1] <- NA
    remaining[remaining==0] <- 1
    remaining <- as.owin(as.im(remaining, W=list(xrange=c(0,2000), yrange=c(0,2000))))
    return(remaining)
  }
 
 # Example:    
 # a <- read.remain.owin("../data/loss_matrices/branching3.csv")
 # W2000 <- as.owin(list(xrange=c(0,2000), yrange=c(0,2000)))
 # sp.pred <- rThomas(0.0001, 15, 13, win=W2000)
 # plot(sp.pred[a])

 
 

  
  