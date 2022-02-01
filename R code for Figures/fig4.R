setwd("C:/Users/Antoine/Dropbox/SurfSweepCodeBlocks/history")

pointSize = 1   ### 2 for tau = 20 and 1 otherwise

N <- 20
io <- 1
tau <- 40
t <- 1:tau

### P(Tf) from simulations
###

name <- paste("TAdistrN20tau",as.character(tau),".dat",sep="")  
TAdistr <- read.table(name)
probaTa <- TAdistr$V3
probaTa <- tail(probaTa,length(probaTa)-1)
sum(probaTa)

plot(1:tau,probaTa,col="black",type="p",xlab=NA,ylab=NA,cex.axis=2,cex=pointSize,ylim=c(0,0.03),pch=19,xlim=c(0,tau))
if(tau == 20)
{
  lines(1:tau,probaTa,col="black")
}




### P(Tf) from analytical trajectory
###

dem <- 2*N*(1-t*(1-io/2/N)/tau)   ### tau < 2Ne

if (tau > 2*N)
{
  dem <- 2*N*(0.5 + 1/2*exp(-t/N) - (1-io/N)/2*exp(-(tau-t)/N)) ##### tau > 2Ne
}

#dem[tau] <- io

probaTa <- t
probaTa[1] <- 1/dem[1]

for (i in 2:tau)
{
  probaTa[i] <- prod(1-1/head(dem,i-1))/dem[i]
}
sum(probaTa)

points(1:tau,probaTa,col="blue",type="p",xlab=NA,ylab=NA,cex.axis=2,cex=pointSize,ylim=c(0,0.053),pch=19,xlim=rev(c(0,120)))
if(tau == 20)
{
  lines(1:tau,probaTa,col="blue")
}



### P(Tf) from simulated trajectory
###


demog <- read.table(paste("demHistN20tau",as.character(tau),".dat",sep=""))
dem <- demog$V2
dem <- c(head(dem,tau-1)[(tau-1):1],io) 

probaTa <- t
probaTa[1] <- 1/dem[1]

for (i in 2:tau)
{
  probaTa[i] <- prod(1-1/head(dem,i-1))/dem[i]
}

points(1:tau,probaTa,col="red",type="p",xlab=NA,ylab=NA,cex.axis=2,cex=pointSize,ylim=c(0,0.053),pch=19,xlim=rev(c(0,120)))
if(tau == 20)
{
  lines(1:tau,probaTa,col="red")
}





