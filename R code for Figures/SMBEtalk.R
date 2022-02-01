

#
##
###
####
#####  compare neutral VS selection troughs


dev.off()
par(mfrow=c(1,1),mar=c(5,5,0.5,0.5))


#### surf

#index <- c(1:5,9)
index <- c(1,2,3,6)
mycolor=viridis(6)
for (il in index)
{
  N <- 20*il
  
  
  io <- 1
  t0 <- 0
  r <- seq(0,0.5,0.0001)
  d <- -log(1-2*r)/2
  
  
  tau <- 150
  
  N0 <- (1500-N)*exp(tau/2/N)+N
  BG <- 3000
  
  t <- 1:tau
  dem <- 2*N*(1-t*(1-io/2/N)/tau)   ### tau < 2Ne
  
  if (tau > 2*N)
  {
    dem <- 2*N*(0.5 + 1/2*exp(-t/N) - (1-io/N)/2*exp(-(tau-t)/N))  ##### tau > 2Ne
  }
  
  dem[tau] <- 1
  
  probaTa <- t
  sum_1minusX <- t
  probaTa[1] <- 1/dem[1]
  sum_1minusX[1] <- 1-dem[1]/2/N
  
  if (tau > 1)
  {
    for (i in 2:tau)
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
  
  c <- 1 - sum(t*probaTa)/BG
  slope <- sum(2*(tau+2*N0-t)*sum_1minusX*probaTa)
  a <- slope/BG/c
  
  exp_fit <-  BG*( 1 - c*exp(-a*r))
  if (il == 1)
  {
    #plot(100*d,exp_fit,lty=1,lwd=2,col="red",xlim=c(4.5,5.5),ylim=c(2400,2850),type="l",cex.axis=2,cex.lab=2,ylab="average coalescence time",xlab="distance (cM)")
    plot(100*d,exp_fit,lty=1,lwd=2,col=mycolor[(1:6)[index==il]],xlim=c(-3,3),ylim=c(0,3000),type="l",cex.axis=2,cex.lab=2,ylab="average coalescence time",xlab="distance (cM)")
  }
  lines(100*d,exp_fit,lty=1,lwd=2,col="red",type="l",cex.axis=2,cex.lab=2,ylab="average coalescence time",xlab="distance (cM)")
  lines(-100*d,exp_fit,lty=1,lwd=2,col="red")
  
  #lines(100*d,TBhom+TBhet*(2*N0+tau),col="black",lwd=2)
  #lines(-100*d,TBhom+TBhet*(3000+tau),col="black",lwd=2)
  
}

### draw a box
# segments(4.5,2400,4.5,2850,lwd=2.5)
# segments(5.5,2400,5.5,2850,lwd=2.5)
# segments(4.5,2400,5.5,2400,lwd=2.5)
# segments(4.5,2850,5.5,2850,lwd=2.5)







####  selection
N1 <- 1500
N <- N1


dem <- 2*N*io/(io+(2*N-io)*exp(-2*log(2*N)/tau*(tau-t)))     ## selection
dem[tau] <- 1


probaTa <- t
sum_1minusX <- t
probaTa[1] <- 1/dem[1]
sum_1minusX[1] <- 1-dem[1]/2/N

if (tau > 1)
{
  for (i in 2:tau)
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

c <- 1 - sum(t*probaTa)/BG
slope <- sum(2*(tau+2*N1-t)*sum_1minusX*probaTa)
a <- slope/BG/c

exp_fit <-  BG*( 1 - c*exp(-a*r))

lines(-100*d,exp_fit,lty=1,lwd=2,col="black",type="l",xlim=c(-10,10))
lines(100*d,exp_fit,lty=1,lwd=2,col="black")

#lines(100*d,TBhom+TBhet*(2*N1+tau),col="black",lwd=2)
#lines(-100*d,TBhom+TBhet*(2*N1+tau),col="black",lwd=2)



#####
####
###
##
#
































#
##
###
####
##### compare sel demog - prob. distr.  VS neutrl demog - prob. distr.




dev.off()
par(mfrow=c(2,2),mar=c(5,5,0.5,0.5))

setwd("C:/Users/Antoine/Dropbox/SurfSweepCodeBlocks/history")

pointSize_tau20 = 1   ### 2 for tau = 20 and 1 otherwise
pointSize_tau40 = 1
pointSize_tau80 = 1
dashed_line_wd = 2

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
plot(c(0,demog$V1),c(1,demog$V2)/40,type="p",col="red",pch=19,xlab="time (in generations)",cex.axis=2,cex.lab=2,ylab="frequency of the allele",cex=pointSize_tau20)
lines(0:20, (1+(0:20)*39/20)/40,lty=2,lwd=dashed_line_wd)
points(c(0,demog$V1),c(1,demog$V2)/40,type="p",col="red",pch=19,xlab=NA,ylab=NA,cex.axis=2,cex.lab=2,cex=pointSize_tau20)

seldem <- 3000*io/(io+(3000-io)*exp(-2*log(3000)/20*(0:20)))
points(0:20,seldem/3000,pch=19)

demog <- read.table("demHistN20tau80.dat")
plot(c(0,demog$V1),c(1,demog$V2)/40,type="p",col="red",pch=19,cex.axis=2,cex.lab=2,cex=pointSize_tau80,xlab="time (in generations)",ylab=NA)
lines((0:80),(0.5 + 1/2*exp(-(80-(0:80))/20) - (1-1/20)/2*exp(-(0:80)/20)),lty=2,lwd=dashed_line_wd)
points(c(0,demog$V1),c(1,demog$V2)/40,type="p",col="red",pch=19,xlab=NA,ylab=NA,cex.axis=2,cex.lab=2,cex=pointSize_tau80)

seldem <- 3000*io/(io+(3000-io)*exp(-2*log(3000)/80*(0:80)))
points(0:80,seldem/3000,pch=19)

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
### selection
###

dem <- 3000*io/(io+(3000-io)*exp(-2*log(3000)/20*(19:0)))


dem[tau] <- io

probaTa <- t
probaTa[1] <- 1/dem[1]

for (i in 2:tau)
{
  probaTa[i] <- prod(1-1/head(dem,i-1))/dem[i]
}
sum(probaTa)


#plot(0:(tau-1),probaTa[tau:1],col="black",type="l",lty=2,lwd=dashed_line_wd,xlab=expression(t~(T[f]~"="~tau~"-"~t)),ylab=expression(distr.~of~coal.~times~P(T[f])),cex.lab=2,cex.axis=2,cex=pointSize,ylim=c(0,0.2),pch=19,xlim=c(0,tau))
plot(0:(tau-1),probaTa[tau:1],col="black",type="p",xlab=expression(t~(T[f]~"="~tau~"-"~t)),ylab="distr. of coal. times P(Tf)",cex.lab=2,cex.axis=2,cex=pointSize_tau40,pch=19,xlim=c(0,tau),ylim=c(0,0.37))



### P(Tf) from simulations
###

name <- paste("TAdistrN20tau",as.character(tau),".dat",sep="")  
TAdistr <- read.table(name)
probaTa <- TAdistr$V3
probaTa <- tail(probaTa,length(probaTa)-1)
sum(probaTa)

points(0:(tau-1),probaTa[tau:1],col="red",type="p",xlab=NA,ylab=NA,cex.axis=2,cex=1,ylim=c(0,0.03),pch=19,xlim=c(0,tau))


#
##
###

tau <- 80
t <- 1:tau
### P(Tf) from analytical trajectory
###

dem <- 3000*io/(io+(3000-io)*exp(-2*log(3000)/80*(79:0)))

#dem[tau] <- io

probaTa <- t
probaTa[1] <- 1/dem[1]

for (i in 2:tau)
{
  probaTa[i] <- prod(1-1/head(dem,i-1))/dem[i]
}
sum(probaTa)


plot(0:(tau-1),probaTa[tau:1],col="black",type="p",xlab=expression(t~(T[f]~"="~tau~"-"~t)),ylab=NA,cex.lab=2,cex.axis=2,ylim=c(0,0.37),pch=19,xlim=c(0,tau))



### P(Tf) from simulations
###

name <- paste("TAdistrN20tau",as.character(tau),".dat",sep="")  
TAdistr <- read.table(name)
probaTa <- TAdistr$V3
probaTa <- tail(probaTa,length(probaTa)-1)
sum(probaTa)

points(0:(tau-1),probaTa[tau:1],col="red",type="p",xlab=NA,ylab=NA,cex.axis=2,cex=1,ylim=c(0,0.03),pch=19,xlim=c(0,tau))















































