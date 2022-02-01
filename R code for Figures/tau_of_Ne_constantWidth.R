N1 <- 1.85*10^6
BG <- 2*N1
io <- 1
t0 <- 0
r <- seq(0,0.5,0.000001)
dev.off()
par(mar=c(5,5,1.1,1))

slope120 <- vector()
depth120 <- vector()

library(viridis)
mycolor = viridis(30)

index <- c(200,500,800,1000,1100,1200,1600,2000,2400,3400,4400,6500,8000,10000,15000)
prior <- c(1028,1661,1932, 2009, 2201, 2245, 2294, 2298, 2287, 2244, 2206, 2154, 2131, 2110, 2079)


tau_of_Ne <- vector()
for (il in index)
{
  N <- il
  #tau <- 1:floor(180/30*N) ##### t50
  #tau <- 200*(1:30)
  tau <- (prior[(1:15)[index==il]]-10):(prior[(1:15)[index==il]]+10)
  
  FWHD <- matrix(nrow=length(io),ncol=length(tau))
  rel_depth <- matrix(nrow=length(io),ncol=length(tau))
  
  for (k in 1:length(io))
  {
    for (j in 1:length(tau))
    {
      N0 <- (N1-N)*exp(tau[j]/2/N)+N
      t <- 1:tau[j]
      dem <- 2*N*(1-t*(1-io[k]/2/N)/tau[j])   ### tau < 2Ne
      
      if (tau[j] > 2*N)
      {
        dem <- 2*N*(0.5 + 1/2*exp(-t/N) - (1-io[k]/N)/2*exp(-(tau[j]-t)/N))  ##### tau > 2Ne
      }
      
      dem[tau[j]] <- 1
      
      probaTa <- t
      sum_1minusX <- t
      probaTa[1] <- 1/dem[1]
      sum_1minusX[1] <- 1-dem[1]/2/N
      
      if (tau[j] > 1)
      {
        for (i in 2:tau[j])
        {
          probaTa[i] <- prod(1-1/head(dem,i-1))/dem[i]
          sum_1minusX[i] <- sum(head(1-dem/2/N,i))
        }
      }
      
      # TBhom <- r
      # TBhet <- r
      # 
      # for (i in 1:length(r))
      # {
      #   TBhom[i] <- sum(t*probaTa*exp(-2*r[i]*sum_1minusX))
      #   TBhet[i] <- 1-sum(exp(-2*r[i]*sum_1minusX)*probaTa)
      # }
      
      c <- 1 - sum(t*probaTa)/BG
      slope <- sum(2*(tau[j]+2*N0-t)*sum_1minusX*probaTa)
      a <- slope/BG/c
      
      exp_fit <-  BG*( 1 - c*exp(-a*r))
      
      #slope[k,j] <- (TBhet[2]*(2*N0+tau[j])+TBhom[2]-TBhet[1]*(2*N0+tau[j])-TBhom[1])/(r[2]-r[1])
      #depth[k,j] <- TBhet[1]*(2*N0+tau[j])+TBhom[1]
      FWHD[k,j] <- -100*log(1-2*r[(1:length(r))[exp_fit>=(BG+sum(t*probaTa))/2][1]])
      rel_depth[k,j] <- 1-sum(t*probaTa)/BG
    }
  }
  
  tau_of_Ne <- c(tau_of_Ne,tau[FWHD[1,]<=0.069][1] )
  
}

plot(index,tau_of_Ne,type="l",xlab="Nc",ylab="tm")
points(index,tau_of_Ne,type="p")
lines(index,index/2,col="red")   ## tm = Nc/2

plot(index,tau_of_Ne/2/index,type="l",xlab="Nc",ylab="tm/2Nc")
points(index,tau_of_Ne/2/index)
abline(h=0.25)


####
###
##
#
#abline(v=0.069,lty=2,lwd=2)  ## (real value 0.06862354) and relative depth value 0.9993096


