##Ridge lines figure


setwd("C:/Users/Antoine/Dropbox/SurfSweepCodeBlocks/coalNiter108")
setwd("C:/Users/Antoine/OneDrive - Universitaet Bern/SurfSweepCodeBlocks/coalNiter108")



pr <- read.table("NucDivDistrT2200mu542e-8.dat")
pr <- read.table("pairProbaDistrT22ago0.dat")
lines(pr$V1[-length(pr$V2)]/15/29/10000,pr$V2[-length(pr$V2)]*5,type="l",col="black",xlim=c(0,0.006))
lines(pr$V1*2*5.42*10^-8,pr$V2,type="l",col="red",xlim=c(0,0.02))

sum(pr$V2*pr$V1)/sum(pr$V2)

sum(pr$V2[pr$V1<152200])/sum(pr$V2)
tail(pr$V2)
(1:length(pr$V2))[pr$V2!=0]
abline(v=15)

sum(pr$V2[-length(pr$V2)]*pr$V1[-length(pr$V2)])/sum(pr$V2[-length(pr$V2)])/15/29/10000  ### divide by width of window and numbers of pairs of haplotypes

#
##
###
####
#####

time <- vector()
popo <- vector()
log_coal <- rep(log(1:10000),20)
linear_coal <- rep(1:10000,20)



### t = 5

pr <- read.table("pairProbaDistrT5.dat")

pr$V2 <- pr$V2/sum(pr$V2)

y0 <- ksmooth(pr$V1,pr$V2,kernel="normal",bandwidth=100,range.x = c(1,length(pr$V1)), n.points=length(pr$V1)-1+1)
x5 <- pr$V1
y5 <- y0$y

plot(x5,y5,type="l")
#plot(pr$V1,pr$V2,col="red")

popo <- c(popo,head(y5,10000))
time <- c(time,rep(1,10000)*5*1)




### t = 10

pr <- read.table("pairProbaDistrT10.dat")

pr$V2 <- pr$V2/sum(pr$V2)

y2 <- ksmooth(pr$V1,pr$V2,kernel="normal",bandwidth=10,range.x = c(10,100), n.points=100-10+1)
y3 <- ksmooth(pr$V1,pr$V2,kernel="normal",bandwidth=100,range.x = c(60,length(pr$V1)),n.points=length(pr$V1)-60+1)

x10 <- pr$V1
y10 <- c(pr$V2[1:10],pr$V2[10]+(1:5)*(y2$y[1]-pr$V2[10])/5,y2$y[y2$x>15],y3$y[y3$x > 100])

plot(log(x10),y10,type="l")
lines(log(pr$V1),pr$V2,col="red")

popo <- c(popo,head(y10,10000))
time <- c(time,rep(1,10000)*5*2)




### t = 15

pr <- read.table("pairProbaDistrT15.dat")

pr$V2 <- pr$V2/sum(pr$V2)

y2 <- ksmooth(pr$V1,pr$V2,kernel="normal",bandwidth=10,range.x = c(15,100), n.points = 100-15+1)
y3 <- ksmooth(pr$V1,pr$V2,kernel="normal",bandwidth=100,range.x = c(60,length(pr$V1)),n.points=length(pr$V1)-60+1)

x15 <- pr$V1
y15 <- c(pr$V2[1:13],pr$V2[13]+(1:8)*(y2$y[8]-pr$V2[13])/9,y2$y[8:length(y2$y)],y3$y[y3$x > 100])

y4 <- ksmooth(x15,y15,kernel="normal",bandwidth=20,range.x = c(50,100), n.points = 100-50+1)

y15 <- c(y15[1:49],y4$y[1:51],y15[101:length(y15)])

#plot(log(x15),y15,type="l",col="blue")
#lines(log(pr$V1),pr$V2,col="red",xlim=c(0,200))

popo <- c(popo,head(y15,10000))
time <- c(time,rep(1,10000)*5*3)



library(ggplot2)
library(ggridges)
library(ggforce)

for (i in 4:20)
{
  name <- paste("pairProbaDistrT",as.character(5*i),".dat",sep="")
  pr <- read.table(name)
  pr$V2 <- pr$V2/sum(pr$V2)
  
  y1 <- ksmooth(pr$V1,pr$V2,kernel="normal",bandwidth=100,range.x = c(100+3*i,length(pr$V1)),n.points=length(pr$V1)-100-3*i+1)
  
  y5i <- c(pr$V2[1:(5*i)],pr$V2[5*i]+(1:(100+3*i-5*i-1))*(y1$y[1]-pr$V2[5*i])/(100+3*i-5*i),y1$y[1:length(y1$y)])
  
  popo <- c(popo,head(y5i,10000))
  time <- c(time,rep(1,10000)*5*i)
}

popo <- log(1000000*popo)

for (i in 1:length(popo))
{
  if (popo[i] < 0)
  {
    popo[i] = 0
  }
}

###### underdominance
time <- vector()
popo <- vector()
linear_coal <- vector()
for (i in 5*(38:80))
{
  time <- c(time,rep(i*100,300))
  #popo <- c(popo,1-w[i*100,])
  popo <- c(popo,aB[i*100,])
  linear_coal <- c(linear_coal,1:300)
}
######


etvla <- data.frame(time,popo,linear_coal)

ggplot(etvla, aes(linear_coal,time, height = popo, group = time))+
  geom_density_ridges(stat = "identity",fill="red",size=0.1,scale=2.5,alpha=0.9)+
  theme_minimal()+
  xlab("space")+
  ylab("time")+
  #scale_x_continuous(trans = trans_reverser('log10'),labels = trans_format("log10", math_format(10^.x)),position="top") +
  #scale_x_log10(limits = c(1, 1e4),labels = trans_format("log10", math_format(10^.x)))+
  #scale_y_continuous(position="right")+
  #annotation_logticks(sides="t",mid = unit(0.3, "cm"),short = unit(0.3, "cm"),trans="reverse")+
  theme(axis.line = element_line(color="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.margin = unit(c(0,0,0.5,0.5), "cm"))
  #theme(axis.text=element_text(size=40),axis.ticks = element_line(size=1),axis.ticks.length = unit(-15, "pt"),axis.text.x.top = element_text(margin=unit(c(0,0,0.8,0.8), "cm"),angle=90,vjust=0.5,hjust=0.5), axis.text.y.right = element_text(margin=unit(c(10,0,0.8,0.8), "cm"),angle=90,hjust=0.5))









 
##### plot diversity along genome
#####
#####


t <- 200

pr <- read.table("pairProbaListT200.dat")

liste <- head(pr$V1,3000)
pos <- 1:3000

opla <- data.frame(liste,pos)

require(scales) ## for trans_format

ggplot(opla, aes(pos,liste))+
  geom_line(size=0.1)+
  scale_y_log10(limits = c(1, 1e4),labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides="l",mid = unit(0.3, "cm"),short = unit(0.3, "cm"))+
  theme_bw()+
  #theme(axis.line = element_line(color="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab(expression(bar(T)))+
  xlab("pos. on the chromosome")+
  geom_hline(size=2,linetype="dashed",col="red",yintercept=(3000-40)*exp(-t/40) +40)+
  geom_hline(size=2,linetype="dashed",col="cyan",yintercept=40)+
  geom_hline(size=2,linetype="dashed",col="green",yintercept=3000)+
  theme(axis.text=element_text(size=40),axis.title=element_text(size=40),axis.title.y=element_text(angle=0,vjust=0.5),axis.ticks = element_line(size=1),axis.ticks.length = unit(-15, "pt"),axis.text.x = element_text(margin=unit(c(0.8,0.8,0.25,0.25), "cm")), axis.text.y = element_text(margin=unit(c(0.8,0.8,0.12,0.12), "cm")))+
  theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm"))

  




