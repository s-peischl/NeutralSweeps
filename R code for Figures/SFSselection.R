setwd("C:/Users/Antoine/OneDrive - Universitaet Bern/SurfingPaper/DmelanogasterSweep/selectiveSweep/lastSims/vcf_files")
library(ggplot2)
library(vcfR)

par(mfcol=c(2,2),mar=c(5.5,5.5,2.5,0.5))


aa = 1.5  ## axis labels size

#
##
###
#### exp fit
N0 <- 1.85*10^6/100 ### ne in the mathematica file
s <- 0.01
muChr <- 5.8*10^-9   ### substitutions per site per generation 
mu <- 3*5.42*10^-10*100  ### effective mutation rate = piChr / (4*ne)   piChr = 0.004
rec <- 3.5*10^-8  ### recombinations per site per generation

d <- 0:140
r <- (1-exp(-d*1000*rec*2))/2

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
####
###
##
#

lower_limit <- c(2/3*exp_fit[97:0],rep(0,7),2/3*exp_fit[0:97])


long_div <- vector()
where <- vector()
SFS_vector1 <- vector()
SFS_vector2 <- vector()
SFS_vector3 <- vector()
SFS_vector4 <- vector()
SFS_vector5 <- vector()
avg_div <- rep(0,300)
avg_fifty <- rep(0,300)



filenamess <- list.files(path = ".", pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
num_troughs <- length(filenamess)


first_one = 0
for (runss in 1:length(filenamess))
{
  num_genes <- 20
  
  vcf <- read.vcfR(filenamess[runss], verbose = FALSE )
  num_ss <- length(vcf@gt[,1])
  
  
  gen <- matrix(nrow = num_ss,ncol = num_genes)
  
  for (j in 1:length(vcf@gt[,1]))
  {
    locvec <- vector()
    for (i in 2:11)
    {
      locvec <- c(locvec,as.integer(unlist(strsplit(vcf@gt[j,i], split="|"))[c(1,3)]))
    }
    gen[j,] <- locvec
  }
  
  print(runss)
  
  
  interval_size <- 10^3
  chrom_length <- 300*10^3
  window_size <- 10*10^3  ## sliding window size
  sfs_window_size <- 10*10^3
  num_windows <- chrom_length/interval_size ## 1kb intervals
  
  avg_dist_ss <- chrom_length/num_ss  ## avg distance in the file (lines) between segregating sites
  
  diff <- rep(0,num_ss)
  cum_diff <- rep(0,num_ss)
  
  div <- rep(0,num_windows)
  fifty <- rep(0,num_windows)
  
  diff <- rowSums(gen)
  diff <- diff*(num_genes-diff)
  
  cum_diff <- cumsum(diff)
  
  search_low <- 1
  search_high <- 1
  
  
  
  for (i in 1:num_windows)
  {
    if (i*interval_size - window_size/2 > as.integer(vcf@fix[1,2]) & i*interval_size + window_size/2 < as.integer(vcf@fix[num_ss,2]))
    {
      while(as.integer(vcf@fix[search_low,2]) >= i*interval_size - window_size/2)
      {
        search_low <- search_low -1
      }
      while (as.integer(vcf@fix[search_low,2])  < i*interval_size - window_size/2)
      {
        search_low <- search_low + 1
      }
      lowerBound <- search_low


      while(as.integer(vcf@fix[search_high,2])  <= i*interval_size + window_size/2)
      {
        search_high <- search_high + 1
      }
      while (as.integer(vcf@fix[search_high,2]) > i*interval_size + window_size/2)
      {
        search_high <- search_high - 1
      }
      upperBound <- search_high

      div[i] <- (cum_diff[upperBound]-cum_diff[lowerBound-1])/(window_size+1)/num_genes/(num_genes-1)*2

      fifty[i] <- length(which(rowSums(gen)[lowerBound:upperBound]==10))/length(rowSums(gen)[lowerBound:upperBound])
    }
  }
  
  
  
  avg_div <- avg_div + div
  #avg_fifty <- avg_fifty + fifty
  
  if (first_one == 0)
  {
    plot(-140:139,div[11:290],type="l",ylim=c(0,0.018),xlim=c(-150,150),xlab="position (Kb)",ylab="nucleotide diversity",col=alpha("black",0.006),lwd=7.5,cex.axis=aa,cex.lab=aa,main="Selective sweeps",cex.main=aa)
    first_one = 1
  }
  lines(-140:139,div[11:290],col=alpha("black",0.006),lwd=7.5)

  #
  ##
  ###
  
  centre_trough <- 150
  
  search_low <- 1
  search_high <- 1

  while(as.integer(vcf@fix[search_low,2]) >= centre_trough*interval_size - sfs_window_size/2)
  {
    search_low <- search_low -1
  }
  while (as.integer(vcf@fix[search_low,2]) < centre_trough*interval_size - sfs_window_size/2)
  {
    search_low <- search_low + 1
  }
  lowerBound <- search_low


  while(as.integer(vcf@fix[search_high,2]) <= centre_trough*interval_size + sfs_window_size/2)
  {
    search_high <- search_high + 1
  }
  while (as.integer(vcf@fix[search_high,2]) > centre_trough*interval_size + sfs_window_size/2)
  {
    search_high <- search_high - 1
  }
  upperBound <- search_high
  SFS_vector1 <- c(SFS_vector1,rowSums(gen)[lowerBound:upperBound])

  
  #
  ##
  ###
  
  while(as.integer(vcf@fix[search_low,2]) >= centre_trough*interval_size - 1.5*sfs_window_size/2)
  {
    search_low <- search_low -1
  }
  while (as.integer(vcf@fix[search_low,2]) < centre_trough*interval_size - 1.5*sfs_window_size/2)
  {
    search_low <- search_low + 1
  }
  lowerBound <- search_low
  
  
  while(as.integer(vcf@fix[search_high,2]) <= centre_trough*interval_size + 1.5*sfs_window_size/2)
  {
    search_high <- search_high + 1
  }
  while (as.integer(vcf@fix[search_high,2]) > centre_trough*interval_size + 1.5*sfs_window_size/2)
  {
    search_high <- search_high - 1
  }
  upperBound <- search_high
  
  SFS_vector2 <- c(SFS_vector2,rowSums(gen)[lowerBound:upperBound])
  
  
  #
  ##
  ###
  
  while(as.integer(vcf@fix[search_low,2]) >= centre_trough*interval_size - 2*sfs_window_size/2)
  {
    search_low <- search_low -1
  }
  while (as.integer(vcf@fix[search_low,2]) < centre_trough*interval_size - 2*sfs_window_size/2)
  {
    search_low <- search_low + 1
  }
  lowerBound <- search_low
  
  
  while(as.integer(vcf@fix[search_high,2]) <= centre_trough*interval_size + 2*sfs_window_size/2)
  {
    search_high <- search_high + 1
  }
  while (as.integer(vcf@fix[search_high,2]) > centre_trough*interval_size + 2*sfs_window_size/2)
  {
    search_high <- search_high - 1
  }
  upperBound <- search_high

  SFS_vector3 <- c(SFS_vector3,rowSums(gen)[lowerBound:upperBound])

  #
  ##
  ###
  
  while(as.integer(vcf@fix[search_low,2]) >= centre_trough*interval_size - 3*sfs_window_size/2)
  {
    search_low <- search_low -1
  }
  while (as.integer(vcf@fix[search_low,2]) < centre_trough*interval_size - 3*sfs_window_size/2)
  {
    search_low <- search_low + 1
  }
  lowerBound <- search_low
  
  
  while(as.integer(vcf@fix[search_high,2]) <= centre_trough*interval_size + 3*sfs_window_size/2)
  {
    search_high <- search_high + 1
  }
  while (as.integer(vcf@fix[search_high,2]) > centre_trough*interval_size + 3*sfs_window_size/2)
  {
    search_high <- search_high - 1
  }
  upperBound <- search_high
  SFS_vector4 <- c(SFS_vector4,rowSums(gen)[lowerBound:upperBound])

  # if (first_one == 0)
  # {
  #   plot(hist(rowSums(gen[,-(1:4)])[lowerBound:upperBound],breaks=1:20-0.5,col=alpha("green",0.5),add=F,freq=T)$density,type="p",ylim=c(0,0.32))
  #   first_one = 1
  #   sfs_cloud <- vector()
  # }
  # lines(hist(rowSums(gen[,-(1:4)])[lowerBound:upperBound],breaks=1:20-0.5,plot=F)$density,type="p")
  # sfs_cloud <- rbind(sfs_cloud,hist(rowSums(gen[,-(1:4)])[lowerBound:upperBound],breaks=1:20-0.5,plot=F)$density)
  # 
  # #
  # ##
  # ###
  # 
  # 
#   while(gen$Pos[search_low] >= (centre_trough+5)*interval_size - sfs_window_size/2)
#   {
#     search_low <- search_low -1
#   }
#   while (gen$Pos[search_low] < (centre_trough+5)*interval_size - sfs_window_size/2)
#   {
#     search_low <- search_low + 1
#   }
#   lowerBound <- search_low
#   
#   
#   while(gen$Pos[search_high] <= (centre_trough+5)*interval_size + sfs_window_size/2)
#   {
#     search_high <- search_high + 1
#   }
#   while (gen$Pos[search_high] > (centre_trough+5)*interval_size + sfs_window_size/2)
#   {
#     search_high <- search_high - 1
#   }
#   upperBound <- search_high
#   SFS_vector2 <- c(SFS_vector2,rowSums(gen[,-(1:4)])[lowerBound:upperBound])
#   
#   ###
#   
#   while(gen$Pos[search_low] >= (centre_trough-5)*interval_size - sfs_window_size/2)
#   {
#     search_low <- search_low -1
#   }
#   while (gen$Pos[search_low] < (centre_trough-5)*interval_size - sfs_window_size/2)
#   {
#     search_low <- search_low + 1
#   }
#   lowerBound <- search_low
#   
#   
#   while(gen$Pos[search_high] <= (centre_trough-5)*interval_size + sfs_window_size/2)
#   {
#     search_high <- search_high + 1
#   }
#   while (gen$Pos[search_high] > (centre_trough-5)*interval_size + sfs_window_size/2)
#   {
#     search_high <- search_high - 1
#   }
#   upperBound <- search_high
#   SFS_vector2 <- c(SFS_vector2,rowSums(gen[,-(1:4)])[lowerBound:upperBound])
#   
#   #
#   ##
#   ###
#   
#   while(gen$Pos[search_low] >= (centre_trough+10)*interval_size - sfs_window_size/2)
#   {
#     search_low <- search_low -1
#   }
#   while (gen$Pos[search_low] < (centre_trough+10)*interval_size - sfs_window_size/2)
#   {
#     search_low <- search_low + 1
#   }
#   lowerBound <- search_low
#   
#   
#   while(gen$Pos[search_high] <= (centre_trough+10)*interval_size + sfs_window_size/2)
#   {
#     search_high <- search_high + 1
#   }
#   while (gen$Pos[search_high] > (centre_trough+10)*interval_size + sfs_window_size/2)
#   {
#     search_high <- search_high - 1
#   }
#   upperBound <- search_high
#   SFS_vector3 <- c(SFS_vector3,rowSums(gen[,-(1:4)])[lowerBound:upperBound])
#   
#   ###
#   
#   while(gen$Pos[search_low] >= (centre_trough-10)*interval_size - sfs_window_size/2)
#   {
#     search_low <- search_low -1
#   }
#   while (gen$Pos[search_low] < (centre_trough-10)*interval_size - sfs_window_size/2)
#   {
#     search_low <- search_low + 1
#   }
#   lowerBound <- search_low
#   
#   
#   while(gen$Pos[search_high] <= (centre_trough-10)*interval_size + sfs_window_size/2)
#   {
#     search_high <- search_high + 1
#   }
#   while (gen$Pos[search_high] > (centre_trough-10)*interval_size + sfs_window_size/2)
#   {
#     search_high <- search_high - 1
#   }
#   upperBound <- search_high
#   SFS_vector3 <- c(SFS_vector3,rowSums(gen[,-(1:4)])[lowerBound:upperBound])
  
  #
  ##
  ###
  
  # while(gen$Pos[search_low] >= (centre_trough+15)*interval_size - sfs_window_size/2)
  # {
  #   search_low <- search_low -1
  # }
  # while (gen$Pos[search_low] < (centre_trough+15)*interval_size - sfs_window_size/2)
  # {
  #   search_low <- search_low + 1
  # }
  # lowerBound <- search_low
  # 
  # 
  # while(gen$Pos[search_high] <= (centre_trough+15)*interval_size + sfs_window_size/2)
  # {
  #   search_high <- search_high + 1
  # }
  # while (gen$Pos[search_high] > (centre_trough+15)*interval_size + sfs_window_size/2)
  # {
  #   search_high <- search_high - 1
  # }
  # upperBound <- search_high
  # SFS_vector4 <- c(SFS_vector4,rowSums(gen[,-(1:4)])[lowerBound:upperBound])
  # 
  # ###
  # 
  # while(gen$Pos[search_low] >= (centre_trough-15)*interval_size - sfs_window_size/2)
  # {
  #   search_low <- search_low -1
  # }
  # while (gen$Pos[search_low] < (centre_trough-15)*interval_size - sfs_window_size/2)
  # {
  #   search_low <- search_low + 1
  # }
  # lowerBound <- search_low
  # 
  # 
  # while(gen$Pos[search_high] <= (centre_trough-15)*interval_size + sfs_window_size/2)
  # {
  #   search_high <- search_high + 1
  # }
  # while (gen$Pos[search_high] > (centre_trough-15)*interval_size + sfs_window_size/2)
  # {
  #   search_high <- search_high - 1
  # }
  # upperBound <- search_high
  # SFS_vector4 <- c(SFS_vector4,rowSums(gen[,-(1:4)])[lowerBound:upperBound])
  
  #
  ##
  ###
  
#   while(gen$Pos[search_low] >= (centre_trough+80)*interval_size - sfs_window_size/2)
#   {
#     search_low <- search_low -1
#   }
#   while (gen$Pos[search_low] < (centre_trough+80)*interval_size - sfs_window_size/2)
#   {
#     search_low <- search_low + 1
#   }
#   lowerBound <- search_low
#   
#   
#   while(gen$Pos[search_high] <= (centre_trough+80)*interval_size + sfs_window_size/2)
#   {
#     search_high <- search_high + 1
#   }
#   while (gen$Pos[search_high] > (centre_trough+80)*interval_size + sfs_window_size/2)
#   {
#     search_high <- search_high - 1
#   }
#   upperBound <- search_high
#   SFS_vector5 <- c(SFS_vector5,rowSums(gen[,-(1:4)])[lowerBound:upperBound])
#   
#   ###
#   
#   while(gen$Pos[search_low] >= (centre_trough-80)*interval_size - sfs_window_size/2)
#   {
#     search_low <- search_low -1
#   }
#   while (gen$Pos[search_low] < (centre_trough-80)*interval_size - sfs_window_size/2)
#   {
#     search_low <- search_low + 1
#   }
#   lowerBound <- search_low
#   
#   
#   while(gen$Pos[search_high] <= (centre_trough-80)*interval_size + sfs_window_size/2)
#   {
#     search_high <- search_high + 1
#   }
#   while (gen$Pos[search_high] > (centre_trough-80)*interval_size + sfs_window_size/2)
#   {
#     search_high <- search_high - 1
#   }
#   upperBound <- search_high
#   SFS_vector5 <- c(SFS_vector5,rowSums(gen[,-(1:4)])[lowerBound:upperBound])
# 


}

lines(0:140,exp_fit,lwd=2.5,col="red")
lines(0:(-140),exp_fit,lwd=2.5,col="red")
#lines(-100:100,lower_limit,lwd=2,lty=2,col="red")
lines(-140:139,avg_div[11:290]/num_troughs,lwd=2,col="green")
#lines(-150:49,avg_fifty/num_troughs,lwd=2,col="green")


text(-120, 0.0175, "(a)",cex = 2)

hist(SFS_vector1[which(SFS_vector1!=20)],breaks=1:20-0.5,col=alpha("green",0.5),add=F,freq=F,xlab="site frequency",ylab="probability",main=NA,cex.axis=aa,cex.lab=aa)
length(SFS_vector1)
points(1:19,1/(1:19)/sum(1/(1:19)),pch=19)
lines(1:19,1/(1:19)/sum(1/(1:19)))

text(3, 0.45, "(b)",cex = 2)

#hist(SFS_vector4,breaks=1:20-0.5,col=alpha("blue",0.5),add=T,freq=F)
#length(SFS_vector5)

#histt <- hist(long_div,breaks = 100,add=F,col=alpha("green",0.2),freq=F,xlim=c(0,0.01),border=alpha("green",0.2)) ## hist of genetic diversity

# gen <- read.table(paste("sim",as.character(centerss$V2[ik]),"_1_1.gen",sep=""),sep="",header=TRUE)
# vectt <-rowSums(gen[,-(1:4)])
# hist(vectt,breaks=1:20-0.5,col=alpha("red",0.5),add=T,freq=F)  ## whole genome SFS
# length(vectt)
# 
# ciii <- vector()
# for (i in 1:19)
# {
#   ciii[i]<-quantile(sfs_cloud[,i],probs=0.95)
# }
# 
# lines(ciii,col="blue",lwd=3)
# lines(colMeans(sfs_cloud),lwd=3,col="red")
