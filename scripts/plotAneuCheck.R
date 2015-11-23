#!/usr/bin/Rscript

#Holmes Rscript to plot aneuploidy check results

#Positional argument (1): full path to Holmes aneuploidyCheck.fractions.txt file (input)
#Positional argument (2): full path to Holmes aneuploidyCheck.counts.txt file (input)
#Positional argument (3): full path to fractions plot pdf (graphical output)
#Positional argument (4): full path to counts plot pdf (graphical output)

#Read command-line arguments
args <- commandArgs(TRUE)

##Ensures Standard Decimal Notation##
options(scipen=1000)

#Load files
frac <- read.table(args[1],header=T)
cop <- read.table(args[2],header=T)

#Extract expected fractions & reformat df for plotting
exp_frac_1 <- frac[,2]/2
exp_frac_2 <- frac[,2]
exp_frac_3 <- (3*frac[,2])/2
frac <- t(frac[,-c(1,2)])

#Plot fractions
pdf(args[3],height=7,width=12)
plot(c(1:24),exp_frac_3,xaxt="n",yaxt="n",ylim=c(0,max(exp_frac_3)),type="n",
     xlab="Chromosome",ylab="Pct. of Reads in Library")
rect(xleft=par("usr")[1],xright=par("usr")[2],
     ybottom=par("usr")[3],ytop=par("usr")[4],
     col="gray90")
for(i in c(1:24)){
  if((i %% 2) > 0){
    rect(xleft=i-0.5,xright=i+0.5,
         ybottom=par("usr")[3],ytop=par("usr")[4],
         col="gray95",border="gray80")
  }
  points(x=jitter(rep(i,nrow(frac)),amount=0.25),y=frac[,i],col=rep(c("orange","forestgreen"),12)[i],cex=0.5)
}
points(c(1:24),exp_frac_3,pch="-",col="blue",cex=2.5)
points(c(1:24),exp_frac_2,pch="-",cex=2.5)
points(c(1:24),exp_frac_1,pch="-",col="red",cex=2.5)
axis(1,at=c(1:24),labels=c(1:22,"X","Y"))
axis(2,at=seq(0,0.12,by=0.02),labels=paste(seq(0,12,by=2),"%"),las=2)
dev.off()

#Reformat copies df for plotting
cop <- t(cop[,-c(1,2)])

#Plot copies
pdf(args[4],height=7,width=12)
plot(c(1:24),xaxt="n",yaxt="n",ylim=c(0,3),type="n",
     xlab="Chromosome",ylab="Projected Copy State")
rect(xleft=par("usr")[1],xright=par("usr")[2],
     ybottom=par("usr")[3],ytop=par("usr")[4],
     col="gray90")
for(i in c(1:24)){
  if((i %% 2) > 0){
    rect(xleft=i-0.5,xright=i+0.5,
         ybottom=par("usr")[3],ytop=par("usr")[4],
         col="gray95",border="gray80")
  }
}
abline(h=c(0,1,2,3),lwd=c(2,1,1,1),col=c("red","red","black","blue"),lty=c(2,1,1,1))
for(i in c(1:24)){
  points(x=jitter(rep(i,nrow(cop)),amount=0.25),y=cop[,i],col=rep(c("orange","forestgreen"),12)[i],cex=0.5)
}
axis(1,at=c(1:24),labels=c(1:22,"X","Y"))
axis(2,at=seq(0,3,by=0.5),las=2)
dev.off()