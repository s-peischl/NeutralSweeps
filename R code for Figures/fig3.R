setwd("C:/Users/Antoine/Dropbox/SurfSweepCodeBlocks/history")

library(EnvStats)

demog <- read.table("demHistN20tau80.dat")
par(mar=c(5.1,5,1,2.1))  ## bottom,left,top,right
plot(c(0,demog$V1),c(1,demog$V2),type="p",col="red",pch=19,cex.axis=2,cex.lab=2,cex=1,xlab="time",ylab="number of mutant copies",xlim=c(0,120))
lines((0:80),40*(0.5 + 1/2*exp(-(80-(0:80))/20) - (1-1/20)/2*exp(-(0:80)/20)),lty=2,lwd=3)
points(c(0,demog$V1),c(1,demog$V2),type="p",col="red",pch=19,xlab=NA,ylab=NA,cex.axis=2,cex.lab=2,cex=1)

demog <- read.table("demHistN20tau10.dat")
points(c(0,demog$V1),c(1,demog$V2),type="p",col="blue",pch=19,xlab=NA,ylab=NA,cex.axis=2,cex.lab=2,cex=1.5)
lines(0:10, (1+(0:10)*39/10),lty=2,lwd=3)
points(c(0,demog$V1),c(1,demog$V2),type="p",col="blue",pch=19,xlab=NA,ylab=NA,cex.axis=2,cex.lab=2,cex=1.5)

demog <- read.table("demHistN20tau40.dat")
points(c(0,demog$V1),c(1,demog$V2),type="p",col="green",pch=19,xlab=NA,ylab=NA,cex.axis=2,cex.lab=2,cex=1.2)
lines(0:40, (1+(0:40)*39/40),lty=2,lwd=3)
points(c(0,demog$V1),c(1,demog$V2),type="p",col="green",pch=19,xlab=NA,ylab=NA,cex.axis=2,cex.lab=2,cex=1.2)

demog <- read.table("demHistN20tau120.dat")
points(c(0,demog$V1),c(1,demog$V2),type="p",col="cyan",pch=19,xlab=NA,ylab=NA,cex.axis=2,cex.lab=2,cex=0.8)
lines((0:120),40*(0.5 + 1/2*exp(-(120-(0:120))/20) - (1-1/20)/2*exp(-(0:120)/20)),lty=2,lwd=3)
points(c(0,demog$V1),c(1,demog$V2),type="p",col="cyan",pch=19,xlab=NA,ylab=NA,cex.axis=2,cex.lab=2,cex=0.8)


