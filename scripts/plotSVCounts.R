#!/usr/bin/Rscript

#Positional argument (1): full path to varCounts_cohort_toPlot.txt (usually in ${WRKDIR})
#Positional argument (2): full path to varCounts_sample_toPlot.txt (usually in ${WRKDIR})
#Positional argument (3): full path to output pdf for grouped  (graphics out)

#Reads args
args <- commandArgs(TRUE)

#Sets params
options(scipen=1000)

#Load interleave helper function
interleave <- function(v1,v2){
  ord1 <- 2*(1:length(v1))-1
  ord2 <- 2*(1:length(v2))
  c(v1,v2)[order(c(ord1,ord2))]
}

#Reads data
cohort <- t(read.table(args[1],header=F))
sample <- t(read.table(args[2],header=F))
#Scales data per sample to same axis as cohort
scaler <- max(apply(cohort,2,sum))/max(apply(sample,2,sum))
max_cohort <- ceiling(max(apply(cohort,2,sum))/100)*100
max_sample <- ceiling(max(apply(sample,2,sum))/10)*10

#Instantiates fake categories in matrix for correct color plotting
x <- matrix(c(cohort[,1],rep(0,28),sample[,1]*scaler,rep(0,28),
              cohort[,2],rep(0,28),sample[,2]*scaler,rep(0,28),
              cohort[,3],rep(0,28),sample[,3]*scaler,rep(0,28),
              cohort[,4],rep(0,28),sample[,4]*scaler,rep(0,28),
              cohort[,5],rep(0,28),sample[,5]*scaler,rep(0,28),
              cohort[,6],rep(0,28),sample[,6]*scaler,rep(0,28),
              cohort[,7],rep(0,28),sample[,7]*scaler),
             byrow=F,nrow=28)
xsum <- matrix(c(sum(cohort[,1]),rep(0,14),sum(sample[,1])*scaler,rep(0,14),
                 sum(cohort[,2]),rep(0,14),sum(sample[,2])*scaler,rep(0,14),
                 sum(cohort[,3]),rep(0,14),sum(sample[,3])*scaler,rep(0,14),
                 sum(cohort[,4]),rep(0,14),sum(sample[,4])*scaler,rep(0,14),
                 sum(cohort[,5]),rep(0,14),sum(sample[,5])*scaler,rep(0,14),
                 sum(cohort[,6]),rep(0,14),sum(sample[,6])*scaler,rep(0,14),
                 sum(cohort[,7]),rep(0,14),sum(sample[,7])*scaler),
               byrow=F,nrow=14)

#Plots barplot to args[3]
pdf(args[3],height=6,width=9)
par(bty="n",mar=c(3.1,4.1,3.1,4.1))
bpt <- barplot(x,space=c(1,0.2),plot=F)
#Create plot area
plot(x=c(0,23),y=c(0,max_cohort),type="n",
     main="Structural Variant Counts",
     ylab="",yaxt="n",xlab="",xaxt="n")
#Background grid lines
segments(x0=par("usr")[1],
         x1=bpt[2],
         y0=seq(0,max_cohort,by=ceiling(max_cohort/600)*100),
         y1=seq(0,max_cohort,by=ceiling(max_cohort/600)*100),
         col="gray20",lty=3)
segments(x0=bpt[13],
         x1=par("usr")[2],
         y0=seq(0,max_sample*scaler,by=scaler*ceiling(max_sample/60)*10),
         y1=seq(0,max_sample*scaler,by=scaler*ceiling(max_sample/60)*10),
         col="gray20",lty=3)
#Add shaded rectangles to distinguish per-cohort and per-sample counts
rect(xleft=interleave(bpt[c(1,3,5,7,9,11,13)]-0.7,bpt[c(1,3,5,7,9,11,13)]+0.6),
     xright=interleave(bpt[c(2,4,6,8,10,12,14)]+0.7,bpt[c(2,4,6,8,10,12,14)]+0.7),
     ybottom=par("usr")[3],ytop=par("usr")[4],
     col=c("gray90","gray95"),border="gray70")
#Hashed axis gridlines within shaded rectanges
segments(x0=rep(bpt[c(1,3,5,7,9,11,13)]-0.7,6),
         x1=rep(bpt[c(1,3,5,7,9,11,13)]+0.6,6),
         y0=seq(0,max_cohort,by=ceiling(max_cohort/600)*100),
         y1=seq(0,max_cohort,by=ceiling(max_cohort/600)*100),
         col="gray60")
segments(x0=rep(bpt[c(1,3,5,7,9,11,13)]+0.6,6),
         x1=bpt[c(2,4,6,8,10,12,14)]+0.7,
         y0=seq(0,max_sample*scaler,by=scaler*ceiling(max_sample/60)*10),
         y1=seq(0,max_sample*scaler,by=scaler*ceiling(max_sample/60)*10),
         col="gray70")
barplot(xsum,space=c(1,0.2),col="white",yaxt="n",add=T)
barplot(x,space=c(1,0.2),
        col=c(rep("firebrick3",2),rep("firebrick1",2),
              rep("dodgerblue3",2),rep("dodgerblue1",2),
              rep("darkorchid4",2),rep("darkorchid2",2),
              rep("darkorange3",2),rep("darkorange1",2),
              rep("forestgreen",2),rep(adjustcolor("forestgreen",alpha=0.4),2),
              rep("aquamarine3",2),rep("aquamarine1",2),
              rep("gray40",2),rep("gray75",2)),
        density=rep(c(NA,20),28),yaxt="n",add=T)
#Y-axes
axis(2,at=seq(0,max_cohort,by=ceiling(max_cohort/600)*100),las=2,cex.axis=0.7)
rightTicks <- seq(0,max_sample*scaler,by=round(scaler*ceiling(max_sample/60)*10,0))
axis(4,at=seq(0,max_sample*scaler,by=round(scaler*ceiling(max_sample/60)*10,0)),
     labels=seq(0,max_sample,by=round(ceiling(max_sample/60)*10,0))[1:length(rightTicks)],
     las=2,cex.axis=0.7)
#Axis labels
mtext(text=c("Deletion","Tandem\nDuplication","Insertion","Simple\nInversion",
             "Translocation","Complex","Incompletely\nResolved"),
      side=1,line=1,cex=0.7,font=3,
      at=c(mean(bpt[1:2]),mean(bpt[3:4]),
           mean(bpt[5:6]),mean(bpt[7:8]),
           mean(bpt[9:10]),mean(bpt[11:12]),
           mean(bpt[13:14])))
mtext("Total Variants in Cohort",side=2,line=2.5)
mtext("Mean Variants per Sample",side=4,line=2.5,srt=180)
legend("topright",legend=c("Total Variants in Cohort","Mean Variants per Sample","VAF > 50%"),
       fill=c("gray20","gray70","black"),density=c(NA,NA,20),
       cex=0.7,bg="white")
dev.off()