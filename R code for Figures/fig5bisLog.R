#setwd("C:/Users/Antoine/Dropbox/SurfSweepCodeBlocks/coalNiter109")
setwd("C:/Users/Antoine/Dropbox/SurfSweepCodeBlocks/coal_high_rec")

par(mfrow=c(1,3))

#rec<-c(0,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001)
#0.002,0.004,0.006,0.008,0.02,0.03,0.04
rec <- c(0,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.002,0.004,0.006,0.008,0.01,0.02,0.03,0.04,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45)

#rec <- (1-exp(-2*d))/2
d <- -log(1-2*rec)/2 ### Haldane's mapping and the centiMorgans log10(log10(100*d))

pointsize <- 1.2

io <- 1
N <- 20


#
##
###
####
#####
######
#######

j <- 15
couleur <- "black"

rel_coal <- 0*rec
sigma <- 0*rec


for (i in 1:length(rec))
{
  name <- paste("coalN20rec","0",substr(formatC(rec[i],format="f"),3,nchar(formatC(rec[i],format="g"))),".dat",sep="")
  if (i==1)
  {
    name <- "coalTimeN20rec0.dat"
  }
  coal <- read.table(name)
  
  k <- 1
  while (coal$V1[k] != j)
  {
    k <- k+1
  }
  print(coal$V1[k])
   
  if (i==1)
  {
    rel_coal[i] <- coal$V3[k]
    sigma[i] <- coal$V4[k]-coal$V3[k]^2  
    sigma[i] <- sqrt(sigma[i])
    sigma[i] <- sigma[i]/(coal$V2[k])^0.5 
    
  }
  
  else
  {
    rel_coal[i] <- coal$V2[k]*(coal$V1[k]+3000)+coal$V4[k]
    sigma[i] <- coal$V5[k]*(3000+coal$V1[k])^2+2*coal$V7[k]*(3000+coal$V1[k])+coal$V6[k]-(coal$V2[k]*(3000+coal$V1[k])+coal$V4[k])^2
    sigma[i] <- sqrt(sigma[i])
    sigma[i] <- sigma[i]/(coal$V3[k])^0.5 
  }
}

plot(log10(100*d),log10(rel_coal),col=couleur,xlab=NA,ylab=NA,type="p",xlim=c(-2,2),ylim=c(1,3.5),pch=19,cex.axis=2,cex.lab=2,cex=pointsize)  ## high rec
#plot(log10(100*d),rel_coal,col=couleur,xlab=NA,ylab=NA,type="p",xlim=c(-2,2),ylim=c(0,500),pch=19,cex.axis=2,cex.lab=2,cex=pointsize)  ## intermediate rec
#plot(log10(100*d),rel_coal,col=couleur,xlab=NA,ylab=NA,type="p",xlim=c(-0.1,0.1),ylim=c(0,60),pch=19,cex.axis=2,cex.lab=2,cex=pointsize) # small rec

#points(-log10(100*d),rel_coal,col=couleur,xlab=NA,ylab=NA,type="p",xlim=c(-0.001,0.001),ylim=c(0,60),pch=19,cex.axis=2,cex.lab=2,cex=pointsize)

#arrows(log10(100*d), rel_coal-2*sigma, log10(100*d), rel_coal + 2*sigma, length=0, angle=90, code=3,col=couleur,lwd=2)
#arrows(-log10(100*d), rel_coal-2*sigma, -log10(100*d), rel_coal + 2*sigma, length=0, angle=90, code=3,col=couleur,lwd=2)


#### analytical with simulated traj and focal coalescence distribution

tau <- j

TAdistr <- read.table(paste("TAdistrN20tau",as.character(j),".dat",sep=""))

probaTa <- TAdistr$V3
probaTa <- tail(probaTa,length(probaTa)-1)

demog <- read.table(paste("demHistN20tau",as.character(j),".dat",sep=""))
dem <- demog$V2
dem <- c(head(dem,tau-1)[(tau-1):1],io) 

t <- 1:tau

sum_1minusX <- t
sum_1minusX[1] <- 1-dem[1]/2/N

for (i in 2:tau)
{
  sum_1minusX[i] <- sum(head(1-dem/2/N,i))
}

r <- seq(0,0.5,0.0001)

TBhom <- r
TBhet <- r

for (i in 1:length(r))
{
  TBhom[i] <- sum(t*probaTa*exp(-2*r[i]*sum_1minusX))
  TBhet[i] <- 1-sum(exp(-2*r[i]*sum_1minusX)*probaTa)
}

points(log10(-100*log(1-2*r)/2),log10(TBhet*(3000+tau)+TBhom),type="l",col=couleur,xlim=c(-0.001,0.001),ylim=c(0,1412),lwd=2)   #### ylim=c(0,3000),xlim=c(-0.5,0.5)
#points(log10(-100*log(1-2*r)/2),TBhet*(3000+tau)+TBhom,type="l",col=couleur,xlim=c(-0.001,0.001),ylim=c(0,1412),lwd=2)


#### analytical

# dem <- 2*N*(1-t*(1-io/2/N)/tau)   
# 
# sum_1minusX <- t
# sum_1minusX[1] <- 1-dem[1]/2/N
# 
# probaTa <- t
# probaTa[1] <- 1/dem[1]
# 
# for (i in 2:tau)
# {
#   probaTa[i] <- prod(1-1/head(dem,i-1))/dem[i]
#   sum_1minusX[i] <- sum(head(1-dem/2/N,i))
# }
# 
# r <- rec
# 
# TBhom <- r
# TBhet <- r
# 
# for (i in 1:length(r))
# {
#   TBhom[i] <- sum(t*probaTa*exp(-2*r[i]*sum_1minusX))
#   TBhet[i] <- 1-sum(exp(-2*r[i]*sum_1minusX)*probaTa)
# }
# 
# points(log10(100*d),TBhet*(3000+tau)+TBhom,type="l",col=couleur,xlim=c(-0.001,0.001),ylim=c(0,1412),lwd=2,lty=2)   #### ylim=c(0,3000),xlim=c(-0.5,0.5)
# points(-log10(100*d),TBhet*(3000+tau)+TBhom,type="l",col=couleur,xlim=c(-0.001,0.001),ylim=c(0,1412),lwd=2,lty=2)

#######
######
#####
####
###
##
#


#
##
###
####
#####
######
#######

j <- 20
couleur <- "red"

rel_coal <- 0*rec
sigma <- 0*rec


for (i in 1:length(rec))
{
  name <- paste("coalN20rec","0",substr(formatC(rec[i],format="f"),3,nchar(formatC(rec[i],format="g"))),".dat",sep="")
  if (i==1)
  {
    name <- "coalTimeN20rec0.dat"
  }
  coal <- read.table(name)
  
  k <- 1
  while (coal$V1[k] != j)
  {
    k <- k+1
  }
  print(coal$V1[k])
  
  if (i==1)
  {
    rel_coal[i] <- coal$V3[k]
    sigma[i] <- coal$V4[k]-coal$V3[k]^2  
    sigma[i] <- sqrt(sigma[i])
    sigma[i] <- sigma[i]/(coal$V2[k])^0.5 
    
  }
  
  else
  {
    rel_coal[i] <- coal$V2[k]*(coal$V1[k]+3000)+coal$V4[k]
    sigma[i] <- coal$V5[k]*(3000+coal$V1[k])^2+2*coal$V7[k]*(3000+coal$V1[k])+coal$V6[k]-(coal$V2[k]*(3000+coal$V1[k])+coal$V4[k])^2
    sigma[i] <- sqrt(sigma[i])
    sigma[i] <- sigma[i]/(coal$V3[k])^0.5 
    
  }
}


plot(log10(100*d),log10(rel_coal),col=couleur,xlab=NA,ylab=NA,type="p",xlim=c(-2,2),ylim=c(1,3.5),pch=19,cex.axis=2,cex.lab=2,cex=pointsize)
#points(-log10(100*d),rel_coal,col=couleur,xlab=NA,ylab=NA,type="p",xlim=c(-0.001,0.001),ylim=c(0,60),pch=19,cex.axis=2,cex.lab=2,cex=pointsize)

#arrows(log10(100*d), log10(rel_coal-2*sigma), log10(100*d), log10(rel_coal + 2*sigma), length=0, angle=90, code=3,col=couleur,lwd=2)
#arrows(-log10(100*d), rel_coal-2*sigma, -log10(100*d), rel_coal + 2*sigma, length=0, angle=90, code=3,col=couleur,lwd=2)


#### analytical with simulated traj and focal coalescence distribution

tau <- j

TAdistr <- read.table(paste("TAdistrN20tau",as.character(j),".dat",sep=""))

probaTa <- TAdistr$V3
probaTa <- tail(probaTa,length(probaTa)-1)

demog <- read.table(paste("demHistN20tau",as.character(j),".dat",sep=""))
dem <- demog$V2
dem <- c(head(dem,tau-1)[(tau-1):1],io) 

t <- 1:tau

sum_1minusX <- t
sum_1minusX[1] <- 1-dem[1]/2/N

for (i in 2:tau)
{
  sum_1minusX[i] <- sum(head(1-dem/2/N,i))
}

r <- seq(0,0.5,0.0001)

TBhom <- r
TBhet <- r

for (i in 1:length(r))
{
  TBhom[i] <- sum(t*probaTa*exp(-2*r[i]*sum_1minusX))
  TBhet[i] <- 1-sum(exp(-2*r[i]*sum_1minusX)*probaTa)
}

points(log10(-100*log(1-2*r)/2),log10(TBhet*(3000+tau)+TBhom),type="l",col=couleur,xlim=c(-0.001,0.001),ylim=c(0,1412),lwd=2)   #### ylim=c(0,3000),xlim=c(-0.5,0.5)
#points(log10(-100*log(1-2*r)/2),TBhet*(3000+tau)+TBhom,type="l",col=couleur,xlim=c(-0.001,0.001),ylim=c(0,1412),lwd=2)


#### analytical

# dem <- 2*N*(1-t*(1-io/2/N)/tau)   
# 
# sum_1minusX <- t
# sum_1minusX[1] <- 1-dem[1]/2/N
# 
# probaTa <- t
# probaTa[1] <- 1/dem[1]
# 
# for (i in 2:tau)
# {
#   probaTa[i] <- prod(1-1/head(dem,i-1))/dem[i]
#   sum_1minusX[i] <- sum(head(1-dem/2/N,i))
# }
# 
# r <- rec
# 
# TBhom <- r
# TBhet <- r
# 
# for (i in 1:length(r))
# {
#   TBhom[i] <- sum(t*probaTa*exp(-2*r[i]*sum_1minusX))
#   TBhet[i] <- 1-sum(exp(-2*r[i]*sum_1minusX)*probaTa)
# }
# 
# points(log10(100*d),TBhet*(3000+tau)+TBhom,type="l",col=couleur,xlim=c(-0.001,0.001),ylim=c(0,1412),lwd=2,lty=2)   #### ylim=c(0,3000),xlim=c(-0.5,0.5)
# points(-log10(100*d),TBhet*(3000+tau)+TBhom,type="l",col=couleur,xlim=c(-0.001,0.001),ylim=c(0,1412),lwd=2,lty=2)

#######
######
#####
####
###
##
#

#
##
###
####
#####
######
#######

j <- 40
couleur <- "blue"

rel_coal <- 0*rec
sigma <- 0*rec


for (i in 1:length(rec))
{
  name <- paste("coalN20rec","0",substr(formatC(rec[i],format="f"),3,nchar(formatC(rec[i],format="g"))),".dat",sep="")
  if (i==1)
  {
    name <- "coalTimeN20rec0.dat"
  }
  coal <- read.table(name)
  
  k <- 1
  while (coal$V1[k] != j)
  {
    k <- k+1
  }
  print(coal$V1[k])
  
  if (i==1)
  {
    rel_coal[i] <- coal$V3[k]
    sigma[i] <- coal$V4[k]-coal$V3[k]^2  
    sigma[i] <- sqrt(sigma[i])
    sigma[i] <- sigma[i]/(coal$V2[k])^0.5 
    
  }
  
  else
  {
    rel_coal[i] <- coal$V2[k]*(coal$V1[k]+3000)+coal$V4[k]
    sigma[i] <- coal$V5[k]*(3000+coal$V1[k])^2+2*coal$V7[k]*(3000+coal$V1[k])+coal$V6[k]-(coal$V2[k]*(3000+coal$V1[k])+coal$V4[k])^2
    sigma[i] <- sqrt(sigma[i])
    sigma[i] <- sigma[i]/(coal$V3[k])^0.5 
    
  }
}


plot(log10(100*d),log10(rel_coal),col=couleur,xlab=NA,ylab=NA,type="p",xlim=c(-2,2),ylim=c(1,3.5),pch=19,cex.axis=2,cex.lab=2,cex=pointsize)
#points(-log10(100*d),rel_coal,col=couleur,xlab=NA,ylab=NA,type="p",xlim=c(-0.001,0.001),ylim=c(0,60),pch=19,cex.axis=2,cex.lab=2,cex=pointsize)

#arrows(log10(100*d), rel_coal-2*sigma, log10(100*d), rel_coal + 2*sigma, length=0, angle=90, code=3,col=couleur,lwd=2)
#arrows(-log10(100*d), rel_coal-2*sigma, -log10(100*d), rel_coal + 2*sigma, length=0, angle=90, code=3,col=couleur,lwd=2)


#### analytical with simulated traj and focal coalescence distribution

tau <- j

TAdistr <- read.table(paste("TAdistrN20tau",as.character(j),".dat",sep=""))

probaTa <- TAdistr$V3
probaTa <- tail(probaTa,length(probaTa)-1)

demog <- read.table(paste("demHistN20tau",as.character(j),".dat",sep=""))
dem <- demog$V2
dem <- c(head(dem,tau-1)[(tau-1):1],io) 

t <- 1:tau

sum_1minusX <- t
sum_1minusX[1] <- 1-dem[1]/2/N

for (i in 2:tau)
{
  sum_1minusX[i] <- sum(head(1-dem/2/N,i))
}

r <- seq(0,0.5,0.0001)

TBhom <- r
TBhet <- r

for (i in 1:length(r))
{
  TBhom[i] <- sum(t*probaTa*exp(-2*r[i]*sum_1minusX))
  TBhet[i] <- 1-sum(exp(-2*r[i]*sum_1minusX)*probaTa)
}

points(log10(-100*log(1-2*r)/2),log10(TBhet*(3000+tau)+TBhom),type="l",col=couleur,xlim=c(-0.001,0.001),ylim=c(0,1412),lwd=2)   #### ylim=c(0,3000),xlim=c(-0.5,0.5)
#points(log10(-100*log(1-2*r)/2),TBhet*(3000+tau)+TBhom,type="l",col=couleur,xlim=c(-0.001,0.001),ylim=c(0,1412),lwd=2)


#### analytical

# dem <- 2*N*(1-t*(1-io/2/N)/tau)   
# 
# sum_1minusX <- t
# sum_1minusX[1] <- 1-dem[1]/2/N
# 
# probaTa <- t
# probaTa[1] <- 1/dem[1]
# 
# for (i in 2:tau)
# {
#   probaTa[i] <- prod(1-1/head(dem,i-1))/dem[i]
#   sum_1minusX[i] <- sum(head(1-dem/2/N,i))
# }
# 
# r <- rec
# 
# TBhom <- r
# TBhet <- r
# 
# for (i in 1:length(r))
# {
#   TBhom[i] <- sum(t*probaTa*exp(-2*r[i]*sum_1minusX))
#   TBhet[i] <- 1-sum(exp(-2*r[i]*sum_1minusX)*probaTa)
# }
# 
# points(log10(100*d),TBhet*(3000+tau)+TBhom,type="l",col=couleur,xlim=c(-0.001,0.001),ylim=c(0,1412),lwd=2,lty=2)   #### ylim=c(0,3000),xlim=c(-0.5,0.5)
# points(-log10(100*d),TBhet*(3000+tau)+TBhom,type="l",col=couleur,xlim=c(-0.001,0.001),ylim=c(0,1412),lwd=2,lty=2)

#######
######
#####
####
###
##
#
