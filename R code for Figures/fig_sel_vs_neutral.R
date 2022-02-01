### FIGURE selection vs neutral trough + frequency + coal. time distr.


#par(mfrow=c(1,2), mar=c(5,5,0.5,0.5))
par(mar=c(5,5,0.5,0.5))
layout(matrix(c(1,2,3,2),ncol=2,byrow=TRUE))

#
##
###
#### surf vs selection trajectory ####


setwd("C:/Users/Antoine/Dropbox/SurfSweepCodeBlocks/history")

tau <- 120
N0 <- 1500
io <- 1
N <- 20

#### surf frequency(t)

#demog <- read.table("demHistN20tau120.dat")
t <- 0:tau
demog <- 2*N*(0.5 + 1/2*exp(-(tau-t)/N) - (1-io/N)/2*exp(-t/N))
if (tau <= 2*N)
{
  demog <- 2*N*t*(1-io/2/N)/tau +1
}

plot(0:tau,demog/40,type="p",col="red",pch=19,cex.axis=2,cex.lab=2,cex=1,xlab="time since mutation",ylab="mutant average frequency")
#lines(0:tau,c(1,demog$V2)/40,type="p",col="blue",pch=19,xlab=NA,ylab=NA,cex.axis=2,cex.lab=2,cex=1)

#### selective frequency(t)


dem <- 2*N0*io/(io+(2*N0-io)*exp(-2*log(2*N0)*t/tau))
points(t,dem/2/N0,col="black",pch=19)

text(5, 0.95, "(a)",cex = 2.5)

####
###
##
#




#
##
###
####   surf vs selective trough ####

#### surf

io <- 1
t0 <- 0
r <- seq(0,0.5,0.0001)
d <- -log(1-2*r)/2


N <- 20

N0 <- (1500-20)*exp(tau/2/N)+20
#N0 <- 1500
BG <- 3000
#BG <- 2*((1500-20)*exp(-tau/2/N)+20)

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
plot(100*d,exp_fit,lty=1,lwd=2,col="red",xlim=c(-5,5),ylim=c(0,3100),type="l",cex.axis=2,cex.lab=2,ylab="average coalescence time",xlab="distance (cM)")
lines(-100*d,exp_fit,lty=1,lwd=2,col="red")

#lines(100*d,TBhom+TBhet*(2*N0+tau),col="black",lwd=2)
#lines(-100*d,TBhom+TBhet*(3000+tau),col="black",lwd=2)


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

lines(-100*d,exp_fit,lty=1,lwd=2,col="black")
lines(100*d,exp_fit,lty=1,lwd=2,col="black")

#lines(100*d,TBhom+TBhet*(2*N1+tau),col="black",lwd=2)
#lines(-100*d,TBhom+TBhet*(2*N1+tau),col="black",lwd=2)

text(-4.5, 3100, "(c)",cex = 2.5)


####
###
##
#






#
##
###
####  surf vs selection focal coalescence time distribution ####

#### selection

t <- 1:tau
N0<-1500
dem <-  2*N0*io/(io+(2*N0-io)*exp(-2*log(2*N0)/tau*(tau-t)))
dem[tau] <- 1

probaTa <- t
probaTa[1] <- 1/dem[1]

if (tau > 1)
{
  for (i in 2:tau)
  {
    probaTa[i] <- prod(1-1/head(dem,i-1))/dem[i]
  }
}

sum(probaTa)
plot(1:tau,probaTa,col="black",type="p",cex.axis=2,cex.lab=2,cex=1,pch=19,xlim=rev(c(0,120)),ylim=c(0,0.055),xlab=expression(coalescence~time~T^(f)),ylab="probability density")
lines(1:tau,probaTa,col="black")

#### surfing
N <- 20
dem <- 2*N*(0.5 + 1/2*exp(-t/N) - (1-io/N)/2*exp(-(tau-t)/N))
dem[tau] <- 1

probaTa[1] <- 1/dem[1]

if (tau > 1)
{
  for (i in 2:tau)
  {
    probaTa[i] <- prod(1-1/head(dem,i-1))/dem[i]
  }
}

sum(probaTa)
points(1:tau,probaTa,col="red",cex=1,pch=19)

text(115, 0.053, "(b)",cex = 2.5)

####
###
##
#


library(pracma)
N1<-1500
tau <- 120
function_tau_to_s <- function(x) tau*x-2*log(4*N1*x)
s <- findzeros(function_tau_to_s,0,1)
s[2]



