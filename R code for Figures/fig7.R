N0_sel <- 1500
io <- 1
t0 <- 0
r <- c(0,0.0000001)

slope1 <- vector()
depth1 <- vector()
slope25 <- vector()
depth25 <- vector()

for (il in 1:5)
{
  N <- 20*il
  tau <- 1:floor(100/30*N) ##### t50
  
  N0_surf <- (N0_sel-N)*exp((t0+tau)/2/N)+N
  
  slope <- matrix(nrow=length(io),ncol=length(tau))
  depth <- matrix(nrow=length(io),ncol=length(tau))
  
  for (k in 1:length(io))
  {
    for (j in 1:length(tau))
    {
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
      
      TBhom <- r
      TBhet <- r
      
      for (i in 1:length(r))
      {
        TBhom[i] <- sum(t*probaTa*exp(-2*r[i]*sum_1minusX))
        TBhet[i] <- 1-sum(exp(-2*r[i]*sum_1minusX)*probaTa)
      }
      slope[k,j] <- (TBhet[2]*(2*N0_surf[j]+tau[j])+TBhom[2]-TBhet[1]*(2*N0_surf[j]+tau[j])-TBhom[1])/(r[2]-r[1])
      depth[k,j] <- TBhet[1]*(2*N0_surf[j]+tau[j])+TBhom[1]
    }
  }
  
  rel_depth <- (2*N+(2*N0_surf-2*N)*exp(-(t0+tau)/2/N)-depth[1,])/(2*N+(2*N0_surf-2*N)*exp(-(t0+tau)/2/N))
  
  tau3 <- tau[(tau-3*(tau%/%3))==1]   #### plot 1 out of 3 value of tau
  
  bridge <- c(max(tau3[tau3 <= 2*N]),min(tau3[tau3 > 2*N]))
  
  
  if (il ==1)
  {
    plot(slope[1,tau3]/4/N0_sel,rel_depth[tau3], xlim=c(0,200),ylim=c(0.92,1),pch=16,cex=1,xlab=NA,ylab=NA,cex.axis=2,col="red")
  }
  
  else
  {
    points(slope[1,tau3]/4/N0_sel,rel_depth[tau3],xlim=c(0,200),ylim=c(0.92,1),pch=16,cex.axis=2,col="red")
  }
  
  lines(slope[1,bridge]/4/N0_sel,rel_depth[bridge],col="red")  #### plot bridge
  
  t1 <- floor(34/30*N)    
  t25 <- floor(72/30*N)
  
  slope1 <- c(slope1,slope[1,t1]/4/N0_sel)
  depth1 <- c(depth1,rel_depth[t1])
  slope25 <- c(slope25,slope[1,t25]/4/N0_sel)
  depth25 <- c(depth25,rel_depth[t25])
  
}

points(slope1,depth1,pch=25,cex=1.5,bg="cyan",col="black")      ###### plot quantiles
points(slope25,depth25,pch=23,cex=1.8,bg="cyan",col="black")

#t1 <- floor(34/30*N)           ##### quantiles 
#t25 <- floor(72/30*N)
#t50 <- floor(100/30*N)
#t75 <- floor(144/30*N)
#t99 <- floor(336/30*N)

####
###
##
#


#
##
###
####  selection

N <- N0_sel
tau <-  1:floor(100/30*100)
N0_surf <- tau*0 + N0_sel

slope <- matrix(nrow=length(io),ncol=length(tau))
depth <- matrix(nrow=length(io),ncol=length(tau))

for (k in 1:length(io))
{
  for (j in 1:length(tau))
  {
    t <- 1:tau[j]
    
    dem <- 2*N*io[k]/(io[k]+(2*N-io[k])*exp(-2*log(2*N)/tau[j]*(tau[j]-t)))     ## selection
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
    
    TBhom <- r
    TBhet <- r
    
    for (i in 1:length(r))
    {
      TBhom[i] <- sum(t*probaTa*exp(-2*r[i]*sum_1minusX))
      TBhet[i] <- 1-sum(exp(-2*r[i]*sum_1minusX)*probaTa)
    }
    slope[k,j] <- (TBhet[2]*(2*N0_surf[j]+tau[j])+TBhom[2]-TBhet[1]*(2*N0_surf[j]+tau[j])-TBhom[1])/(r[2]-r[1])
    depth[k,j] <- TBhet[1]*(2*N0_surf[j]+tau[j])+TBhom[1]
  }
}

rel_sel_depth <- (2*N0_sel-depth[1,])/(2*N0_sel)
tau3 <- tau[(tau-3*(tau%/%3))==1]
points(slope[1,tau3]/4/N0_sel,rel_sel_depth[tau3],pch=16,col="black")

lines(1:300, 1-(1:300)/N0_sel,col="black",lty=2,lwd=3)  ### star tree approximation


taus <- c(40,100,180,275)
points(slope[1,taus]/4/N0_sel,rel_sel_depth[taus],pch=21,bg="green",cex=2)   ### plot values of s

library(pracma)
for (i in 1:length(taus))
{
  tau <- taus[i]
  function_tau_to_s <- function(x) tau*x-2*log(4*N0_sel*x)
  s <- findzeros(function_tau_to_s,0,1)
  print(paste("for tau = ",taus[i]," we have s = ",s[2],sep=""))
}

####
###
##
#


#
##
###
####   plot N0_surf/N0_sel

for (il in 1:5)
{
  N <- 20*il
  tau <- 1:floor(100/30*N) ##### t50
  
  N0_surf <- (N0_sel-N)*exp((t0+tau)/2/N)+N
  if (il == 1)
  {
    plot(tau,N0_surf/N0_sel,xlim=c(0,350),xlab=NA,ylab=NA,cex.axis=2)
  }
  else
  {
    points(N0_surf/N0_sel,xlab=NA,ylab=NA,cex.axis=2)
  }
}

####
###
##
#