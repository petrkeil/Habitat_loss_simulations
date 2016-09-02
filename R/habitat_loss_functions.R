
require(spatstat)


# ------------------------------------------------------------------------------
# FUNCTION FITTING THE THOMAS MODEL TO XY POINT DATA

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
# FIT AND EXTRAPOLATE THE THOMAS MODEL  
# Function that fits, extrapolates (to a 2x2 km area) and plots the 
# Thomas model for a specific species:

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