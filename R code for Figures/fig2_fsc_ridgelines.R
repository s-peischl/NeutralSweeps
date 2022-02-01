library(ggplot2)
library(ggridges)
library(ggforce)
library(ggpubr)
library(ggprism)

#### panel a ####

setwd("C:/Users/Antoine/OneDrive - Universitaet Bern/SurfingPaper/DmelanogasterSweep")

etvla <- read.csv("fsc_ridgelines.csv")

ggridge_yo <- ggplot(etvla, aes(time,linear_coal, width = popo, group=time))+
  #geom_density_ridges(stat="identity",size=0.1,scale=2,aes(fill=time))+ ## different colours
  geom_vridgeline(stat="identity",size=0.1,scale=4,fill="lightblue",alpha=0.8)+
  theme_minimal()+
  xlab(expression(time~since~contraction~("in"~2*N[c]~generations)))+
  ylab(expression(nucleotide~ diversity~pi))+
  #scale_fill_viridis(alpha=0.9)+  ## different colours
  #scale_x_continuous(trans = scales::reverse_trans(),labels = trans_format("log10", math_format(10^.x)),position="top") +
  scale_y_continuous(breaks=c(0,0.002,0.004,0.00514726,0.006),labels=c(expression(4*N[c]~mu),"0.002","0.004",expression(4*N[0]~mu),"0.006"))+ 
  #scale_x_log10(limits = c(1, 1e4),labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous(breaks=c(0,2200,4400,6600,8800,11100),labels=c("0","0.25","0.5","0.75","1","1.25"))+
  #annotation_logticks(sides="t",mid = unit(0.3, "cm"),short = unit(0.3, "cm"),trans="reverse")+
  theme(axis.line = element_line(color="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  #theme(plot.margin = unit(c(0,0,0.5,0.5), "cm"))+
  theme(axis.text=element_text(size=20),axis.ticks = element_line(size=0.5),axis.ticks.length = unit(-8, "pt"),axis.text.x.top = element_text(margin=unit(c(0,0,0.8,0.8), "cm"),angle=90,vjust=0.5,hjust=0.5), axis.text.y.right = element_text(margin=unit(c(10,0,0.8,0.8), "cm"),angle=90,hjust=0.5))+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text.x = element_text(margin=unit(c(0.3,0.8,0.25,0.12), "cm")), axis.text.y = element_text(margin=unit(c(0.8,0.2,0.12,0.12), "cm")))+
  annotate("text", x=1100, y=0.007, label= "(a)",size=10) 
  #annotate("text", x=7500, y=0.0027, label= "(c)",size=10) +  
  #annotate("text", x = 3200, y=0.004, label = "(b)",size=10)

#
##
###
#### panel b ####
##### 

setwd("C:/Users/Antoine/OneDrive - Universitaet Bern/SurfingPaper/DmelanogasterSweep/N1_2_37e6_t_0")


long_div <- vector()
where <- vector()
mean_div <- 0
for (ik in 1:1)
{
  gen <- read.table(paste("rangeExpansion2D-c1-1x1_1_",as.character(ik),".gen",sep=""),sep="",header=TRUE)
  #gen <- read.table("N0_1_85e6_no_bneck.gen",sep="",header=TRUE)
  
  
  print(ik)
  
  num_genes <- 30
  interval_size <- 10^3
  chrom_length <- 2*10^7
  #chrom_length <- 10^6
  window_size <- 10^4  ## sliding window size
  num_windows <- chrom_length/interval_size ## 1kb intervals
  
  num_ss <- length(gen$Pos)    ##  number of segregating sites
  avg_dist_ss <- chrom_length/num_ss  ## avg distance in the file (lines) between segregating sites
  
  diff <- rep(0,length(gen$Pos))
  cum_diff <- rep(0,length(gen$Pos))
  
  div <- rep(0,num_windows)
  
  diff <- rowSums(gen[,-(1:4)])
  #diff <- rowSums(gen[,-(1:34)]) # if both populations are sampled
  diff <- diff*(num_genes-diff)
  
  #long_div <- c(long_div,sum(diff)/10000/num_genes/(num_genes-1)*2)  ### for 10kbp genome
  
  
  cum_diff <- cumsum(diff)
  head(cum_diff)
  
  mean_div <- mean_div + cum_diff[length(cum_diff)]/num_genes/(num_genes-1)*2/chrom_length
  
  # hist(num_of_1s[174923:175032],breaks=100)
  # gen$Pos[174923]
  # gen$Pos[175032]
  
  
  #length(num_of_1s[num_of_1s==0 | num_of_1s == num_genes])
  
  search_low <- 1
  search_high <- 1
  
  for (i in 1:num_windows)
    #for (i in 20) ## for 10kbp genome
  {
    #i <- num_windows
    if (i*interval_size - window_size/2 > gen$Pos[1] && i*interval_size + window_size/2 < gen$Pos[length(gen$Pos)])
    {
      #print(paste(as.character(ik)," ok"))  ## 10kbp genome
      
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
      
      # if (i %% 200 == 0)
      # {
      #   print(i)
      #   #print(floor((i*interval_size - window_size/2)/avg_dist_ss) - lowerBound)
      #   #print(lowerBound-upperBound)
      #   #print(gen$Pos[upperBound]-gen$Pos[lowerBound])
      # }
      
      div[i] <- (cum_diff[upperBound]-cum_diff[lowerBound-1])/(window_size+1)/num_genes/(num_genes-1)*2
    }
  }
  
  long_div <- c(long_div,div[20:(length(div)-20)])
  #long_div <- c(long_div,div[20])   ## 10kbp genome
  where <- c(where,rep(ik,num_windows-39))
}

t <- 0
mu=5.42*10^-10

liste <- long_div
pos <- 20:(length(div)-20)

opla <- data.frame(liste,pos/1000)

require(scales) ## for trans_format


ggDiv0 <- ggplot(opla, aes(pos/1000,liste))+
  geom_line(size=0.1)+
  scale_y_continuous(limits=c(0,0.007),breaks=c(0,0.002,0.004,0.00514726,0.006),labels=c(expression(4*N[c]~mu),"0.002","0.004",expression(4*N[0]~mu),"0.006"))+ 
  #annotation_logticks(sides="l",mid = unit(0.3, "cm"),short = unit(0.3, "cm"))+
  theme_bw()+
  annotate("text", x=3, y=0.007, label= "(b)",size=10) +
  annotate("text", x=4, y=0.002, label= expression(italic(t)[c]==0),size=7) +
  scale_x_continuous(limits=c(0,20),breaks=c(0,5,10,15,20),minor_breaks=c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19),labels = c("0","5","10","15","20"))+
  annotation_ticks(sides="b",type="both",outside=F,minor.length = unit(7,"pt"),tick.length = unit(13,"pt"),size=0.5)+
  #theme(axis.line = element_line(color="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab(expression(nucleotide~ diversity~pi))+
  xlab("pos. on the chromosome (in Mb)")+
  geom_hline(size=1,linetype="dashed",col="red",yintercept=2*mu*((4748395-8800)*exp(-t/8800) +8800))+
  #geom_hline(size=2,linetype="dashed",col="cyan",yintercept=40)+
  #geom_hline(size=2,linetype="dashed",col="green",yintercept=3000)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20),axis.title.y=element_text(angle=90,vjust=0.5),axis.ticks = element_line(size=0.5),axis.ticks.length.y = unit(-13, "pt"),axis.ticks.length.x = unit(0, "pt"),axis.text.x = element_text(margin=unit(c(0.3,0.8,0.25,0.12), "cm")), axis.text.y = element_text(margin=unit(c(0.8,0.2,0.12,0.12), "cm")))+
  
  theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm"))


#
##
###
#### panel c ####
#####

setwd("C:/Users/Antoine/OneDrive - Universitaet Bern/SurfingPaper/DmelanogasterSweep/N1_2_37e6_100")


long_div <- vector()
where <- vector()
mean_div <- 0
for (ik in 84:84)
{
  gen <- read.table(paste("rangeExpansion2D-c1-1x1_1_",as.character(ik),".gen",sep=""),sep="",header=TRUE)
  #gen <- read.table("N0_1_85e6_no_bneck.gen",sep="",header=TRUE)
  
  
  print(ik)
  
  num_genes <- 30
  interval_size <- 10^3
  chrom_length <- 2*10^7
  #chrom_length <- 10^6
  window_size <- 10^4  ## sliding window size
  num_windows <- chrom_length/interval_size ## 1kb intervals
  
  num_ss <- length(gen$Pos)    ##  number of segregating sites
  avg_dist_ss <- chrom_length/num_ss  ## avg distance in the file (lines) between segregating sites
  
  diff <- rep(0,length(gen$Pos))
  cum_diff <- rep(0,length(gen$Pos))
  
  div <- rep(0,num_windows)
  
  diff <- rowSums(gen[,-(1:4)])
  #diff <- rowSums(gen[,-(1:34)]) # if both populations are sampled
  diff <- diff*(num_genes-diff)
  
  #long_div <- c(long_div,sum(diff)/10000/num_genes/(num_genes-1)*2)  ### for 10kbp genome
  
  
  cum_diff <- cumsum(diff)
  head(cum_diff)
  
  mean_div <- mean_div + cum_diff[length(cum_diff)]/num_genes/(num_genes-1)*2/chrom_length
  
  # hist(num_of_1s[174923:175032],breaks=100)
  # gen$Pos[174923]
  # gen$Pos[175032]
  
  
  #length(num_of_1s[num_of_1s==0 | num_of_1s == num_genes])
  
  search_low <- 1
  search_high <- 1
  
  for (i in 1:num_windows)
    #for (i in 20) ## for 10kbp genome
  {
    #i <- num_windows
    if (i*interval_size - window_size/2 > gen$Pos[1] && i*interval_size + window_size/2 < gen$Pos[length(gen$Pos)])
    {
      #print(paste(as.character(ik)," ok"))  ## 10kbp genome
      
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
      
      # if (i %% 200 == 0)
      # {
      #   print(i)
      #   #print(floor((i*interval_size - window_size/2)/avg_dist_ss) - lowerBound)
      #   #print(lowerBound-upperBound)
      #   #print(gen$Pos[upperBound]-gen$Pos[lowerBound])
      # }
      
      div[i] <- (cum_diff[upperBound]-cum_diff[lowerBound-1])/(window_size+1)/num_genes/(num_genes-1)*2
    }
  }
  
  long_div <- c(long_div,div[20:(length(div)-20)])
  #long_div <- c(long_div,div[20])   ## 10kbp genome
  where <- c(where,rep(ik,num_windows-39))
}

t <- 2200
mu=5.42*10^-10

liste <- long_div
pos <- 20:(length(div)-20)

opla <- data.frame(liste,pos/1000)

require(scales) ## for trans_format


ggDiv2200 <-  ggplot(opla, aes(pos/1000,liste))+
  geom_line(size=0.1)+
  scale_y_continuous(limits=c(0,0.007),breaks=c(0,0.002,0.004,0.00514726,0.006),labels=c(expression(4*N[c]~mu),"0.002","0.004",expression(4*N[0]~mu),"0.006"))+ 
  #annotation_logticks(sides="l",mid = unit(0.3, "cm"),short = unit(0.3, "cm"))+
  theme_bw()+
  scale_x_continuous(limits=c(0,20),breaks=c(0,5,10,15,20),minor_breaks=c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19),labels = c("0","5","10","15","20"))+
  annotation_ticks(sides="b",type="both",outside=F,minor.length = unit(7,"pt"),tick.length = unit(13,"pt"),size=0.5)+
  #theme(axis.line = element_line(color="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab(expression(nucleotide~ diversity~pi))+
  xlab("pos. on the chromosome (in Mb)")+
  geom_hline(size=1,linetype="dashed",col="red",yintercept=2*mu*((4748395-8800)*exp(-t/8800) +8800))+
  #geom_hline(size=2,linetype="dashed",col="cyan",yintercept=40)+
  #geom_hline(size=2,linetype="dashed",col="green",yintercept=3000)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20),axis.title.y=element_text(angle=90,vjust=0.5),axis.ticks = element_line(size=0.5),axis.ticks.length.y = unit(-13, "pt"),axis.ticks.length.x = unit(0, "pt"),axis.text.x = element_text(margin=unit(c(0.3,0.8,0.25,0.12), "cm")), axis.text.y = element_text(margin=unit(c(0.8,0.2,0.12,0.12), "cm")))+
  annotate("text", x=3, y=0.0065, label= "(c)",size=10) +
  annotate("text", x=18, y=0.0065, label= expression(italic(t)[c]==0.25),size=7) +
  theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm"))

  ## axis.ticks.length.x.bottom = unit(c(rep(c(1,0.5,0.5,0.5,0.5),4),1)*15,"pt")
  ## axis.ticks.length = unit(-15, "pt")

#
##
###
#### panel d ####
#####

setwd("C:/Users/Antoine/OneDrive - Universitaet Bern/SurfingPaper/DmelanogasterSweep/N1_2_37e6_t_6600")


long_div <- vector()
where <- vector()
mean_div <- 0
for (ik in 1:1)
{
  gen <- read.table(paste("rangeExpansion2D-c1-1x1_1_",as.character(ik),".gen",sep=""),sep="",header=TRUE)
  #gen <- read.table("N0_1_85e6_no_bneck.gen",sep="",header=TRUE)
  
  
  print(ik)
  
  num_genes <- 30
  interval_size <- 10^3
  chrom_length <- 2*10^7
  #chrom_length <- 10^6
  window_size <- 10^4  ## sliding window size
  num_windows <- chrom_length/interval_size ## 1kb intervals
  
  num_ss <- length(gen$Pos)    ##  number of segregating sites
  avg_dist_ss <- chrom_length/num_ss  ## avg distance in the file (lines) between segregating sites
  
  diff <- rep(0,length(gen$Pos))
  cum_diff <- rep(0,length(gen$Pos))
  
  div <- rep(0,num_windows)
  
  diff <- rowSums(gen[,-(1:4)])
  #diff <- rowSums(gen[,-(1:34)]) # if both populations are sampled
  diff <- diff*(num_genes-diff)
  
  #long_div <- c(long_div,sum(diff)/10000/num_genes/(num_genes-1)*2)  ### for 10kbp genome
  
  
  cum_diff <- cumsum(diff)
  head(cum_diff)
  
  mean_div <- mean_div + cum_diff[length(cum_diff)]/num_genes/(num_genes-1)*2/chrom_length
  
  # hist(num_of_1s[174923:175032],breaks=100)
  # gen$Pos[174923]
  # gen$Pos[175032]
  
  
  #length(num_of_1s[num_of_1s==0 | num_of_1s == num_genes])
  
  search_low <- 1
  search_high <- 1
  
  for (i in 1:num_windows)
    #for (i in 20) ## for 10kbp genome
  {
    #i <- num_windows
    if (i*interval_size - window_size/2 > gen$Pos[1] && i*interval_size + window_size/2 < gen$Pos[length(gen$Pos)])
    {
      #print(paste(as.character(ik)," ok"))  ## 10kbp genome
      
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
      
      # if (i %% 200 == 0)
      # {
      #   print(i)
      #   #print(floor((i*interval_size - window_size/2)/avg_dist_ss) - lowerBound)
      #   #print(lowerBound-upperBound)
      #   #print(gen$Pos[upperBound]-gen$Pos[lowerBound])
      # }
      
      div[i] <- (cum_diff[upperBound]-cum_diff[lowerBound-1])/(window_size+1)/num_genes/(num_genes-1)*2
    }
  }
  
  long_div <- c(long_div,div[20:(length(div)-20)])
  #long_div <- c(long_div,div[20])   ## 10kbp genome
  where <- c(where,rep(ik,num_windows-39))
}

t <- 6600
mu=5.42*10^-10

liste <- long_div
pos <- 20:(length(div)-20)

opla <- data.frame(liste,pos/1000)

require(scales) ## for trans_format


ggDiv6600 <- ggplot(opla, aes(pos/1000,liste))+
  geom_line(size=0.1)+
  scale_y_continuous(limits=c(0,0.007),breaks=c(0,0.002,0.004,0.00514726,0.006),labels=c(expression(4*N[c]~mu),"0.002","0.004",expression(4*N[0]~mu),"0.006"))+ 
  #annotation_logticks(sides="l",mid = unit(0.3, "cm"),short = unit(0.3, "cm"))+
  theme_bw()+
  scale_x_continuous(limits=c(0,20),breaks=c(0,5,10,15,20),minor_breaks=c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19),labels = c("0","5","10","15","20"))+
  annotation_ticks(sides="b",type="both",outside=F,minor.length = unit(7,"pt"),tick.length = unit(13,"pt"),size=0.5)+
  #theme(axis.line = element_line(color="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab(expression(nucleotide~ diversity~pi))+
  xlab("pos. on the chromosome (in Mb)")+
  geom_hline(size=1,linetype="dashed",col="red",yintercept=2*mu*((4748395-8800)*exp(-t/8800) +8800))+
  #geom_hline(size=2,linetype="dashed",col="cyan",yintercept=40)+
  #geom_hline(size=2,linetype="dashed",col="green",yintercept=3000)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20),axis.title.y=element_text(angle=90,vjust=0.5),axis.ticks = element_line(size=0.5),axis.ticks.length.y = unit(-13, "pt"),axis.ticks.length.x = unit(0, "pt"),axis.text.x = element_text(margin=unit(c(0.3,0.8,0.25,0.12), "cm")), axis.text.y = element_text(margin=unit(c(0.1,0.2,0.12,0.12), "cm")))+
  annotate("text", x=3, y=0.0065, label= "(d)",size=10) +
  annotate("text", x=18, y=0.0065, label= expression(italic(t)[c]==0.75),size=7) +
  theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm"))


ggarrange(ggridge_yo, ggDiv0, ggDiv2200, ggDiv6600, heights = c(2, 2),ncol = 2, nrow = 2)



