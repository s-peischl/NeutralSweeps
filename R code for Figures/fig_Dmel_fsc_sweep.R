setwd("C:/Users/Antoine/OneDrive - Universitaet Bern/SurfingPaper/DmelanogasterSweep/N1_2_37e6_100")

library(ggplot2)

long_div <- vector()
where <- vector()
mean_div <- 0


for (ik in 84)
{
  gen <- read.table(paste("rangeExpansion2D-c1-1x1_1_",as.character(ik),".gen",sep=""),sep="",header=TRUE)
  
  print(ik)
  
  num_genes <- 30
  interval_size <- 10^3
  chrom_length <- 2*10^7
  window_size <- 10^4  ## sliding window size
  num_windows <- chrom_length/interval_size ## 1kb intervals
  
  num_ss <- length(gen$Pos)    ##  number of segregating sites
  avg_dist_ss <- chrom_length/num_ss  ## avg distance in the file (lines) between segregating sites
  
  diff <- rep(0,length(gen$Pos))
  cum_diff <- rep(0,length(gen$Pos))
  
  div <- rep(0,num_windows)
  
  diff <- rowSums(gen[,-(1:4)])
  diff <- diff*(num_genes-diff)
  
  cum_diff <- cumsum(diff)
  head(cum_diff)
  
  mean_div <- mean_div + cum_diff[length(cum_diff)]/num_genes/(num_genes-1)*2/chrom_length
  
  search_low <- 1
  search_high <- 1
  
  for (i in 1:num_windows)
  {
    if (i*interval_size - window_size/2 > gen$Pos[1] && i*interval_size + window_size/2 < gen$Pos[length(gen$Pos)])
    {
      while(gen$Pos[search_low] >= i*interval_size - window_size/2)
      {
        search_low <- search_low -1
      }
      while (gen$Pos[search_low] < i*interval_size - window_size/2)
      {
        search_low <- search_low + 1
      }
      lowerBound <- search_low
      
      
      while(gen$Pos[search_high] <= i*interval_size + window_size/2)
      {
        search_high <- search_high + 1
      }
      while (gen$Pos[search_high] > i*interval_size + window_size/2)
      {
        search_high <- search_high - 1
      }
      upperBound <- search_high
      
      div[i] <- (cum_diff[upperBound]-cum_diff[lowerBound-1])/(window_size+1)/num_genes/(num_genes-1)*2
    }
  }
  
  long_div <- c(long_div,div[20:(length(div)-20)])
  where <- c(where,rep(ik,num_windows-39))
}


histt <- hist(long_div,breaks = 100,add=F,col=alpha("green",0.2),freq=F,xlim=c(0,0.01),border=alpha("green",0.2))



#
##
###
####
##### FIGURE 8 D. mel sweep

dev.off()
layout(matrix(c(1,1,1,1,2), nrow=1))  ### in 100 run 84


par(mar=c(6,6.1,1.5,0.5))
plot(-100:100,div[19570:19770],ylim=c(0,0.006),xlab=NA,ylab=NA,type="l",col="blue",cex.axis=3,cex.lab=3,axes=F)  ## sim 84 in 100
axis(1,mgp=c(4,2,0),cex.axis=3)
axis(3,at=c(-1000,1000))
axis(1,at=c(-1000,1000))
axis(2,cex.axis=3)
axis(4,at=c(-1,1))
axis(2,at=c(-1,1))
title(ylab="nucleotide diversity",cex.lab=3,line=4)
title(xlab="position on the chromosome (Kb)",cex.lab=3,line=4.5)




#
##
### compute sel sweep exp_fit

N0 <- 1.85*10^6   ### ne in the mathematica file
s <- 0.0098
muChr <- 5.8*10^-9   ### substitutions per site per generation 
mu <- 5.42*10^-10  ### effective mutation rate = piChr / (4*ne)   piChr = 0.004
rec <- 3.5*10^-8  ### recombinations per site per generation

r <- seq(0,0.5,10^-4)
d <- -log(1-2*r)/2/rec/1000  ### distance in kbp  

#tau_s <- floor(2*log(4*N0*s)/s)  ## Barton correction
tau_s <- floor(2*log(2*N0)/s)     ## without correction gives a result closer to Rogers et al.

t <- 1:tau_s

dem <- 2*N0/(1+(2*N0-1)*exp(-2*log(2*N0)/tau_s*(tau_s-t)))     ## selection
dem[tau_s] <- 1

sum_1minusX <- t
sum_1minusX[1] <- 1-dem[1]/2/N0

probaTa <- t
probaTa[1] <- 1/dem[1]

for (i in 2:tau_s)
{
  probaTa[i] <- prod(1-1/head(dem,i-1))/dem[i]
  sum_1minusX[i] <- sum(head(1-dem/2/N0,i))
}

BG <- 2*N0
c <- 1 - sum(t*probaTa)/BG
slope <- sum(2*(tau_s+2*N0-t)*sum_1minusX*probaTa)
a <- slope/BG/c

exp_fit <-  2*mu*BG*( 1 - c*exp(-a*r))

###
##
#

lines(d,exp_fit,lty=1,lwd=1.5,col="red")  
lines(-d,exp_fit,lty=1,lwd=1.5,col="red")


#
##
### compute neutral sweep exp fit

N <- 4400
tau <- 2200
N1 <- 1.85*10^6
N0 <- (N1-N)*exp(tau/2/N)+N
t <- 1:tau
dem <- 2*N*(1-t*(1-1/2/N)/tau)   ### tau < 2Ne

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

c <- 1 - sum(t*probaTa)/BG
slope <- sum(2*(tau+2*N0-t)*sum_1minusX*probaTa)
a <- slope/BG/c

exp_fit <-  2*mu*BG*( 1 - c*exp(-a*r))

###
##
#

lines(d,exp_fit,lty=2,lwd=3,col="black")  
lines(-d,exp_fit,lty=2,lwd=3,col="black")


par(mar=c(6,0.5,1.5,1.5))
plot(histt$density,histt$mids,type="l",xlab = NA,ylab=NA,xaxt="n",yaxt="n",bty="n")
polygon( c(0,histt$density,0), c(histt$mids[1],histt$mids,histt$mids[length(histt$mids)]),col=gray(0.8),border=NA)
lines(histt$density,histt$mids,type="l")

#####
####
###
##
#