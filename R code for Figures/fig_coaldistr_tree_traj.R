## run until line line 290 then run arbre.R for tau = 20, 40 and 80 and then run this one from line 208 to the end

dev.off()
library(ggplot2)
par(mfrow=c(2,4),mar=c(6,6.5,3,1))

setwd("C:/Users/Antoine/OneDrive - Universitaet Bern/SurfSweepCodeBlocks/history")

a = 1.5

pointSize_tau20 = 1.2*a   ### 2 for tau = 20 and 1 otherwise
pointSize_tau40 = 0.8*a
pointSize_tau80 = 0.6*a
dashed_line_wd = 2*a


N <- 20
io <- 1
#tau <- 20
#t <- 1:tau

#
##
###
####
##### first row

demog <- read.table("demHistN20tau20.dat")
plot(c(0,demog$V1),c(1,demog$V2),type="p",col="red",pch=19,main=expression(t[m]==20),cex.main=2*a,xlab=NA,cex.axis=2*a,cex.lab=2*a,ylab=NA,cex=pointSize_tau20,xaxt="n")
#plot(c(0,demog$V1),c(1,demog$V2),type="p",col="red",pch=19,main=expression(t[m]==20),cex.main=2*a,xlab="time since mutation",cex.axis=2*a,cex.lab=2*a,ylab="number of mutant copies",cex=pointSize_tau20)
axis(1,at=5*(0:4),labels=c("0","5","10","15","20"),cex.axis=2*a,padj=0.7)
title(xlab="time since mutation",cex.lab=2*a,line=5)
title(ylab="number of mutant copies",cex.lab=2*a,line=4)


demog <- read.table("demHistN20tau20Distr.dat")
for (i in 1:length(demog[,1]))
{
  lines(0:20,c(1,demog[i,]),col=alpha("black",0.006),lwd=7.5)
}

demog <- read.table("demHistN20tau20.dat")
lines(0:20, (1+(0:20)*39/20),lty=2,lwd=dashed_line_wd,col="white")
points(c(0,demog$V1),c(1,demog$V2),type="p",col="red",pch=19,xlab=NA,ylab=NA,cex.axis=2*a,cex.lab=2*a,cex=pointSize_tau20)
text(3, 30, "(a)",cex = 3*a)



demog <- read.table("demHistN20tau40.dat")
plot(c(0,demog$V1),c(1,demog$V2),type="p",col="red",pch=19,xlab=NA,main=expression(t[m]==40),cex.main=2*a,ylab=NA,cex.axis=2*a,cex.lab=2*a,cex=pointSize_tau40,xaxt="n")
axis(1,at=10*(0:4),labels=c("0","10","20","30","40"),cex.axis=2*a,padj=0.7)
title(xlab="time since mutation",cex.lab=2*a,line=5)


demog <- read.table("demHistN20tau40Distr.dat")
for (i in 1:1777)
{
  lines(0:40,c(1,demog[i,]),col=alpha("black",0.006),lwd=7.5)
}

demog <- read.table("demHistN20tau40.dat")
lines(0:40, (1+(0:40)*39/40),lty=2,lwd=dashed_line_wd,col="white")
points(c(0,demog$V1),c(1,demog$V2),type="p",col="red",pch=19,xlab="time since mutation",ylab=NA,cex.axis=2,cex.lab=2,cex=pointSize_tau40)
text(6, 30, "(b)",cex = 3*a)


demog <- read.table("demHistN20tau80.dat")
plot(c(0,demog$V1),c(1,demog$V2),type="p",col="red",pch=19,cex.axis=2*a,cex.lab=2*a,main=expression(t[m]==80),cex.main=2*a,cex=pointSize_tau80,xlab=NA,ylab=NA,xaxt="n")
axis(1,at=20*(0:4),labels=c("0","20","40","60","80"),cex.axis=2*a,padj=0.7)
title(xlab="time since mutation",cex.lab=2*a,line=5)

demog <- read.table("demHistN20tau80Distr.dat")
for (i in 1:1777)
{
  lines(0:80,c(1,demog[i,]),col=alpha("black",0.007),lwd=4)
}

demog <- read.table("demHistN20tau80.dat")
lines((0:80),40*(0.5 + 1/2*exp(-(80-(0:80))/20) - (1-1/20)/2*exp(-(0:80)/20)),lty=2,lwd=dashed_line_wd,col="white")
points(c(0,demog$V1),c(1,demog$V2),type="p",col="red",pch=19,xlab=NA,ylab=NA,cex.axis=2,cex.lab=2,cex=pointSize_tau80)
text(12, 30, "(c)",cex = 3*a)

demog <- read.table("demHistN20tau120.dat")
plot(c(0,demog$V1),c(1,demog$V2),type="p",col="red",pch=19,cex.axis=2*a,cex.lab=2*a,main=expression(t[m]==120),cex.main=2*a,cex=pointSize_tau80,xlab=NA,ylab=NA,xaxt="n")
axis(1,at=30*(0:4),labels=c("0","30","60","90","120"),cex.axis=2*a,padj=0.7)
title(xlab="time since mutation",cex.lab=2*a,line=5)

demog <- read.table("demHistN20tau120Distr.dat")
for (i in 1:1777)
{
  lines(0:120,c(1,demog[i,]),col=alpha("black",0.007),lwd=3)
}

demog <- read.table("demHistN20tau120.dat")
lines((0:120),40*(0.5 + 1/2*exp(-(120-(0:120))/20) - (1-1/20)/2*exp(-(0:120)/20)),lty=2,lwd=dashed_line_wd,col="white")
points(c(0,demog$V1),c(1,demog$V2),type="p",col="red",pch=19,xlab=NA,ylab=NA,cex.axis=2,cex.lab=2,cex=pointSize_tau80)
text(12, 30, "(d)",cex = 3*a)

#
##
###
####
##### second row

#
##
### 

tau <- 20
t <- 1:tau
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


#plot(0:(tau-1),probaTa[tau:1],col="black",type="l",lty=2,lwd=dashed_line_wd,xlab=expression(t~(T[f]~"="~tau~"-"~t)),ylab=expression(distr.~of~coal.~times~P(T[f])),cex.lab=2,cex.axis=2,cex=pointSize,ylim=c(0,0.2),pch=19,xlim=c(0,tau))
plot(1:tau,probaTa,col="black",type="l",lty=2,lwd=dashed_line_wd,xlab=NA,ylab=expression(distr.~of~ coal. ~times~ P(T^(f) )),cex.lab=2*a,cex.axis=2*a,cex=pointSize_tau40,ylim=c(0,0.2),pch=19,xlim=rev(c(0,tau)),xaxt="n")
text(17, 0.15, "(e)",cex = 3*a)
axis(1,at=5*(0:4),labels=c("0","5","10","15","20"),cex.axis=2*a,padj=0.7)
title(xlab=expression(T^(f)),cex.lab=2*a,line=5)


### P(Tf) from simulations
###

name <- paste("TAdistrN20tau",as.character(tau),".dat",sep="")  
TAdistr <- read.table(name)
probaTa <- TAdistr$V3
probaTa <- tail(probaTa,length(probaTa)-1)
sum(probaTa)

points(1:tau,probaTa,col="red",type="p",xlab=NA,ylab=NA,cex.axis=2,cex=pointSize_tau20,ylim=c(0,0.03/2),pch=19,xlim=c(0,tau))


#
##
###

tau <- 40
t <- 1:tau
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


plot(1:tau,probaTa,col="black",type="l",lty=2,lwd=dashed_line_wd,xlab=NA,ylab=NA,cex.lab=2*a,cex.axis=2*a,cex=pointSize_tau40,ylim=c(0,0.2/2),pch=19,xlim=rev(c(0,tau)),xaxt="n")
text(34, 0.075, "(f)",cex = 3*a)
axis(1,at=10*(0:4),labels=c("0","10","20","30","40"),cex.axis=2*a,padj=0.7)
title(xlab=expression(T^(f)),cex.lab=2*a,line=5)


### P(Tf) from simulations
###

name <- paste("TAdistrN20tau",as.character(tau),".dat",sep="")  
TAdistr <- read.table(name)
probaTa <- TAdistr$V3
probaTa <- tail(probaTa,length(probaTa)-1)
sum(probaTa)

points(1:tau,probaTa,col="red",type="p",xlab=NA,ylab=NA,cex.axis=2,cex=pointSize_tau40,ylim=c(0,0.03),pch=19,xlim=c(0,tau))


#
##
###

tau <- 80
t <- 1:tau
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


plot(1:tau,probaTa,col="black",type="l",lty=2,lwd=dashed_line_wd,xlab=NA,ylab=NA,cex.lab=2*a,cex.axis=2*a,ylim=c(0,0.2/4),pch=19,xlim=rev(c(0,tau)),xaxt="n")
text(68, 0.0375, "(g)",cex = 3*a)
axis(1,at=20*(0:4),labels=c("0","20","40","60","80"),cex.axis=2*a,padj=0.7)
title(xlab=expression(T^(f)),cex.lab=2*a,line=5)


### P(Tf) from simulations
###

name <- paste("TAdistrN20tau",as.character(tau),".dat",sep="")  
TAdistr <- read.table(name)
probaTa <- TAdistr$V3
probaTa <- tail(probaTa,length(probaTa)-1)
sum(probaTa)

points(1:tau,probaTa,col="red",type="p",xlab=NA,ylab=NA,cex.axis=2,cex=pointSize_tau80,ylim=c(0,0.03),pch=19,xlim=c(0,tau))

#
##
###

tau <- 120
t <- 1:tau
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


plot(1:tau,probaTa,col="black",type="l",lty=2,lwd=dashed_line_wd,xlab=NA,ylab=NA,cex.lab=2*a,cex.axis=2*a,ylim=c(0,0.2/6),pch=19,xlim=rev(c(0,tau)),xaxt="n")
text(102, 0.0250, "(h)",cex = 3*a)
axis(1,at=30*(0:4),labels=c("0","30","60","90","120"),cex.axis=2*a,padj=0.7)
title(xlab=expression(T^(f)),cex.lab=2*a,line=5)


### P(Tf) from simulations
###

name <- paste("TAdistrN20tau",as.character(tau),".dat",sep="")  
TAdistr <- read.table(name)
probaTa <- TAdistr$V3
probaTa <- tail(probaTa,length(probaTa)-1)
sum(probaTa)

points(1:tau,probaTa,col="red",type="p",xlab=NA,ylab=NA,cex.axis=2,cex=pointSize_tau80,ylim=c(0,0.03),pch=19,xlim=c(0,tau))























## skip ####


# setwd("C:/Users/Antoine/Dropbox/SurfSweepCodeBlocks/history")
# 
# #
# ##
# ###
# #### creates an inset plot
# 
# par(fig = c(0.4,0.65, 0.43, 0.63), new = T) 
# 
# tau <- 40
# t <- 1:tau
# 
# dem <- 2*N*(1-t*(1-io/2/N)/tau)   ### tau < 2Ne
# 
# if (tau > 2*N)
# {
#   dem <- 2*N*(0.5 + 1/2*exp(-t/N) - (1-io/N)/2*exp(-(tau-t)/N)) ##### tau > 2Ne
# }
# 
# #dem[tau] <- io
# 
# probaTa <- t
# probaTa[1] <- 1/dem[1]
# 
# for (i in 2:tau)
# {
#   probaTa[i] <- prod(1-1/head(dem,i-1))/dem[i]
# }
# sum(probaTa)
# 
# 
# plot(0:(tau-1),probaTa[tau:1],col="black",type="l",lty=2,lwd=dashed_line_wd,xlab=NA,ylab=NA,cex.lab=2,cex.axis=2,ylim=c(0.015,0.03),pch=19,xlim=c(0,tau))
# 
# 
# 
# ### P(Tf) from simulations
# ###
# 
# name <- paste("TAdistrN20tau",as.character(tau),".dat",sep="")  
# TAdistr <- read.table(name)
# probaTa <- TAdistr$V3
# probaTa <- tail(probaTa,length(probaTa)-1)
# sum(probaTa)
# 
# points(0:(tau-1),probaTa[tau:1],col="red",type="p",xlab=NA,ylab=NA,cex.axis=2,cex=pointSize_tau40,ylim=c(0,0.03),pch=19,xlim=c(0,tau))
# 
# 
# #
# ##
# ###
# #### creates an inset plot
# 
# par(fig = c(0.73,0.98, 0.43, 0.63), new = T)
# 
# tau <- 80
# t <- 1:tau
# 
# dem <- 2*N*(1-t*(1-io/2/N)/tau)   ### tau < 2Ne
# 
# if (tau > 2*N)
# {
#   dem <- 2*N*(0.5 + 1/2*exp(-t/N) - (1-io/N)/2*exp(-(tau-t)/N)) ##### tau > 2Ne
# }
# 
# #dem[tau] <- io
# 
# probaTa <- t
# probaTa[1] <- 1/dem[1]
# 
# for (i in 2:tau)
# {
#   probaTa[i] <- prod(1-1/head(dem,i-1))/dem[i]
# }
# sum(probaTa)
# 
# 
# plot(0:(tau-1),probaTa[tau:1],col="black",type="l",lty=2,lwd=dashed_line_wd,xlab=NA,ylab=NA,cex.lab=2,cex.axis=2,ylim=c(0,0.03),pch=19,xlim=c(0,tau))
# 
# 
# 
# ### P(Tf) from simulations
# ###
# 
# name <- paste("TAdistrN20tau",as.character(tau),".dat",sep="")
# TAdistr <- read.table(name)
# probaTa <- TAdistr$V3
# probaTa <- tail(probaTa,length(probaTa)-1)
# sum(probaTa)
# 
# points(0:(tau-1),probaTa[tau:1],col="red",type="p",xlab=NA,ylab=NA,cex.axis=2,cex=pointSize_tau80,ylim=c(0,0.03),pch=19,xlim=c(0,tau))





















### P(Tf) from simulated trajectory
###


# demog <- read.table(paste("demHistN20tau",as.character(tau),".dat",sep=""))
# dem <- demog$V2
# dem <- c(head(dem,tau-1)[(tau-1):1],io) 
# 
# probaTa <- t
# probaTa[1] <- 1/dem[1]
# 
# for (i in 2:tau)
# {
#   probaTa[i] <- prod(1-1/head(dem,i-1))/dem[i]
# }
# 
# points(1:tau,probaTa,col="red",type="p",xlab=NA,ylab=NA,cex.axis=2,cex=pointSize,ylim=c(0,0.053),pch=19,xlim=rev(c(0,120)))





