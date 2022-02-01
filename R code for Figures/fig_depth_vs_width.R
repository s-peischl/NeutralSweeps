N1 <- 1500
BG <- 2*N1
io <- 1
t0 <- 0
r <- seq(0,0.5,0.00001)
dev.off()

### read timeDistr30 from SurfSweepCodeBlocks/coalescent 

setwd("C:/Users/Antoine/OneDrive - Universitaet Bern/SurfSweepCodeBlocks/coalescent")
probb <- read.table("timeDistr30.dat")
mysize <- probb$V2
mysize=0.99/300+probb$V2*0 



#par(mfrow=c(1,2),mar=c(5,5,0.5,0.5))
par(mfrow=c(1,1),mar=c(5,5,0.5,0.5))

slope120 <- vector()
depth120 <- vector()

library(viridis)
mycolor = viridis(floor(100/30*100))
mycolor <- mycolor[floor(100/30*100):1]

index <- c(1:5,9)
#index <- c(1,2,3,6)
for (il in index)
{
  N <- 20*il
  #tau <- 1:floor(180/30*N) ##### t50
  tau <- 1:floor(100/30*100) ##### t50
  
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
  
  # rel_depth <- (2*N+(2*N0-2*N)*exp(-(t0+tau)/2/N)-depth[1,])/(2*N+(2*N0-2*N)*exp(-(t0+tau)/2/N))
  # 
  # tau3 <- tau[(tau-3*(tau%/%3))==1]   #### plot 1 out of 3 value of tau
  # 
  # bridge <- c(max(tau3[tau3 <= 2*N]),min(tau3[tau3 > 2*N]))
  
  for (j in 1:length(tau))
  {
    if (il == 1 && j == 1)
    {
      plot(FWHD[1,j],rel_depth[1,j],pch=16,cex=mysize[floor(j*30/il/20)+(floor(j*30/il/20)==0)]*300+0.01,xlab="width (cM)",ylab="relative depth",cex.lab=2,cex.axis=2,col=mycolor[j],ylim=c(0.91,1),log="x",xlim=c(0.001,200),xaxt="n")
      #axis(2,cex.axis=2,cex.lab=2)
      #axis(2,at=c(0,2))
      #axis(4,at=c(0,2))
      axis(1,at=c(0.01,1,100),cex.axis=2)
      axis(3,at=c(0.0001,1000))
      axis(1,at=c(0.0001,0.001,0.1,10,1000),labels=NA)
      axis(1,at=0.01*c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),tcl=-0.3,labels=NA)
      axis(1,at=0.1*c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),tcl=-0.3,labels=NA)
      axis(1,at=c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),tcl=-0.3,labels=NA)
      axis(1,at=10*c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),tcl=-0.3,labels=NA)
      axis(1,at=100*c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),tcl=-0.3,labels=NA)
      axis(1,at=1000*c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),tcl=-0.3,labels=NA)
      axis(1,at=10000*c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),tcl=-0.3,labels=NA)
      axis(1,at=10000*c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),tcl=-0.3,labels=NA)
      #plot(FWHD[1,j],rel_depth[1,j],pch=16,cex=1,xlab=NA,ylab=NA,cex.axis=2,col=mycolor[j],xlim=c(0,150),ylim=c(0.9,1))
      #plot(j,log10(FWHD[1,j]),pch=16,cex=1,xlab=NA,ylab=NA,cex.axis=2,col=mycolor[j],xlim=c(0,600),ylim=c(-1,5))
    }
    
    else
    {
      points(FWHD[1,j],rel_depth[1,j],pch=16,cex=mysize[floor(j*30/il/20)+(floor(j*30/il/20)==0)]*300+0.01,xlab=NA,ylab=NA,cex.axis=2,col=mycolor[j])
      #points(FWHD[1,j],rel_depth[1,j],pch=16,cex=1,xlab=NA,ylab=NA,cex.axis=2,col=mycolor[j])
      #points(j,log10(FWHD[1,j]),pch=16,cex=1,xlab=NA,ylab=NA,cex.axis=2,col=mycolor[j],xlim=c(0,600),ylim=c(0,10))
    }
  }
  # jj <- 150
  # points(FWHD[1,jj],rel_depth[1,jj],pch=16,cex=2.5,xlab=NA,ylab=NA,cex.axis=2,col="red")
  
  
  
  # lines(slope[1,bridge]/4/N0,rel_depth[bridge],col="red")  #### plot bridge
  # 
  # slope120 <- c(slope120,slope[1,120]/4/N0)    #### tau = 120
  # depth120 <- c(depth120,rel_depth[120])
}

####
###
##
#




#
##
###
####  selection

N <- N1
tau <-  1:floor(100/30*100)

FWHD <- matrix(nrow=length(io),ncol=length(tau))
rel_depth <- matrix(nrow=length(io),ncol=length(tau))

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
    
    # TBhom <- r
    # TBhet <- r
    # 
    # for (i in 1:length(r))
    # {
    #   TBhom[i] <- sum(t*probaTa*exp(-2*r[i]*sum_1minusX))
    #   TBhet[i] <- 1-sum(exp(-2*r[i]*sum_1minusX)*probaTa)
    # }
    
    c <- 1 - sum(t*probaTa)/BG
    slope <- sum(2*(tau[j]+2*N1-t)*sum_1minusX*probaTa)
    a <- slope/BG/c
    
    exp_fit <-  BG*( 1 - c*exp(-a*r))
    
    #slope[k,j] <- (TBhet[2]*(2*N0+tau[j])+TBhom[2]-TBhet[1]*(2*N0+tau[j])-TBhom[1])/(r[2]-r[1])
    #depth[k,j] <- TBhet[1]*(2*N0+tau[j])+TBhom[1]
    FWHD[k,j] <- -100*log(1-2*r[(1:length(r))[exp_fit>=(BG+sum(t*probaTa))/2][1]])
    rel_depth[k,j] <- 1-sum(t*probaTa)/BG
  }
}


for (j in 1:length(tau))
{
  points(FWHD[1,j],rel_depth[1,j],pch=16,cex=1,xlab=NA,ylab=NA,cex.axis=2,col=mycolor[j])
}
# jj <- 150
# points(FWHD[1,jj],rel_depth[1,jj],pch=16,cex=2.5,xlab=NA,ylab=NA,cex.axis=2,col="black")

# rel_sel_depth <- (2*N0-depth[1,])/(2*N0)
# tau3 <- tau[(tau-3*(tau%/%3))==1]
# points(slope[1,tau3]/4/N0,rel_sel_depth[tau3],pch=16,col="black")
# 
# slope120 <- c(slope120,slope[1,120]/4/N0)    #### tau = 120
# depth120 <- c(depth120,rel_sel_depth[120])
# 
# 
# lines(1:300, 1-(1:300)/N0,col="black",lty=2,lwd=3)  ### star tree approximation
# 
# slope120 <- c(slope120,60)    #### tau = 120
# depth120 <- c(depth120,1-60/N0)
# 
# 
# lines(slope120,depth120,col="black",lwd=2)
# points(slope120,depth120,pch=21,cex=1.5,bg="cyan",col="black")  ### tau = 120

####
###
##
#

library(plotrix)
#color.legend(13,0.92,20,0.95,align="rb",c(" 1 (23.5)",rep(" ",48)," 50 (0.3)",rep(" ",99)," 150 (0.08)",rep(" ",99)," 250 (0.045)",rep(" ",82)," 333 (0.03)"),mycolor,gradient="y",cex=1.5)
color.legend(13,0.92,20,0.95,align="rb",rev(c(rep(" ",49)," 50",rep(" ",99)," 150",rep(" ",99)," 250",rep(" ",82)," 333")),rev(mycolor),gradient="y",cex=1.5)
color.legend(13,0.92,20,0.95,align="lt",rev(c(rep(" ",49),"0.3 ",rep(" ",99),"0.08 ",rep(" ",99),"0.045 ",rep(" ",82),"0.03 ")),rev(mycolor),gradient="y",cex=1.5)
text(17, 0.955, expression(s~~~~~italic(t)[m]),cex = 2)

text_size = 1.8
text(0.003, 0.995, expression(N[c] == 20),cex = text_size)
text(0.02, 0.9825, expression(N[c] == 40),cex = text_size)
text(0.051, 0.974, expression(N[c] == 60),cex = text_size)
text(0.071, 0.966, expression(N[c] == 80),cex = text_size)
text(0.087, 0.956, expression(N[c] == 100),cex = text_size)
text(0.12, 0.942, expression(N[c] == 180),cex = text_size)
text(0.2, 0.9185, "selection",cex = text_size)

taus <- c(1,20,80,50,150,250,333)
library(pracma)
for (i in 1:length(taus))
{
  tau <- taus[i]
  function_tau_to_s <- function(x) tau*x-2*log(4*N1*x)
  s <- findzeros(function_tau_to_s,0,1)
  print(paste("for tau = ",taus[i]," we have s = ",s[2],sep=""))
}




























#########  depth vs width D. mealnogaster
########
#######
######
#####
####
###
##
#


N1 <- 1.85*10^6
BG <- 2*N1
io <- 1
t0 <- 0
r <- seq(0,0.5,0.000001)
#dev.off()
par(mar=c(5,5,1.1,1))

slope120 <- vector()
depth120 <- vector()

library(viridis)
mycolor = viridis(30)

index <- c(800,1200,1600,2400,4400,15000,1500000)
for (il in 1:6)
{
  N <- index[il]
  #tau <- 1:floor(180/30*N) ##### t50
  tau <- 200*(1:30)
  
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
      #FWHD[k,j] <- -100*log(1-2*r[(1:length(r))[exp_fit>=(BG+sum(t*probaTa))/2][1]])
      FWHD[k,j] <- -100*log(1-2*log(2)/a)  ### r_hd = log(2)/a
      rel_depth[k,j] <- c
    }
  }
  
  # rel_depth <- (2*N+(2*N0-2*N)*exp(-(t0+tau)/2/N)-depth[1,])/(2*N+(2*N0-2*N)*exp(-(t0+tau)/2/N))
  # 
  # tau3 <- tau[(tau-3*(tau%/%3))==1]   #### plot 1 out of 3 value of tau
  # 
  # bridge <- c(max(tau3[tau3 <= 2*N]),min(tau3[tau3 > 2*N]))
  
  for (j in 1:length(tau))
  {
    if (il == 1 && j == 1)
    {
      plot(FWHD[1,j],rel_depth[1,j],pch=16,cex=1,xlab="FWHD (cM)",ylab="relative depth",cex.lab=2,cex.axis=2,col=mycolor[j],ylim=c(0.9985,1),log="x",xlim=c(0.005,1),xaxt="n")  ## ylim=(0.9985,1)
      #axis(2,cex.axis=2,cex.lab=2)
      #axis(2,at=c(0,2))
      #axis(4,at=c(0,2))
      axis(1,at=c(0.001,0.01,0.1,1,10),cex.axis=2)
      axis(3,at=c(0.001,10))
      axis(1,at=0.1*c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),tcl=-0.3,labels=NA)
      axis(1,at=c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),tcl=-0.3,labels=NA)
      axis(1,at=10*c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),tcl=-0.3,labels=NA)
      axis(1,at=100*c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),tcl=-0.3,labels=NA)
      axis(1,at=1000*c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),tcl=-0.3,labels=NA)
      #plot(FWHD[1,j],rel_depth[1,j],pch=16,cex=1,xlab=NA,ylab=NA,cex.axis=2,col=mycolor[j],xlim=c(0,150),ylim=c(0.9,1))
      #plot(j,log10(FWHD[1,j]),pch=16,cex=1,xlab=NA,ylab=NA,cex.axis=2,col=mycolor[j],xlim=c(0,600),ylim=c(-1,5))
    }
    
    else
    {
      points(FWHD[1,j],rel_depth[1,j],pch=16,cex=1,xlab=NA,ylab=NA,cex.axis=2,col=mycolor[j])
      #points(FWHD[1,j],rel_depth[1,j],pch=16,cex=1,xlab=NA,ylab=NA,cex.axis=2,col=mycolor[j])
      #points(j,log10(FWHD[1,j]),pch=16,cex=1,xlab=NA,ylab=NA,cex.axis=2,col=mycolor[j],xlim=c(0,600),ylim=c(0,10))
    }
  }
  
  
  
  # lines(slope[1,bridge]/4/N0,rel_depth[bridge],col="red")  #### plot bridge
  # 
  # slope120 <- c(slope120,slope[1,120]/4/N0)    #### tau = 120
  # depth120 <- c(depth120,rel_depth[120])
}

####
###
##
#
abline(v=0.069,lty=2,lwd=2)  ## (real value 0.06862354) and relative depth value 0.9993096

library(plotrix)
color.legend(0.2,0.9987,0.25,0.9992,align="rb",c(200,rep(" ",3),1000,rep(" ",4),2000,rep(" ",4),3000,rep(" ",4),4000,rep(" ",4),5000,rep(" ",4),6000),mycolor,gradient="y",cex=1.5)
text(0.22, 0.99925, expression(tau),cex = 2)

text_size = 1.8
text(0.007, 0.99976, "Nc = 800",cex = text_size)
text(0.008, 0.999615, "Nc = 1200",cex = text_size)
text(0.009, 0.99945, "Nc = 1600",cex = text_size)
text(0.011, 0.99929, "Nc = 2400",cex = text_size)
text(0.013, 0.99904, "Nc = 4400",cex = text_size)
text(0.011, 0.99865, "Nc = 1.5e+4",cex = text_size)

abline(h=0.99931)
abline(h=0.999525)
