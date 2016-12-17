#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#Plot variant frequencies from Holmes

#Positional argument (1): full path to frequencies directory (subdir of ${WRKDIR})
#Positional argument (2): full path to output pdf (graphics out)

#Reads args
args <- commandArgs(TRUE)

#Sets params
options(scipen=1000)

#Reads data
del <- as.vector(read.table(paste(args[1],"deletion.list",sep=""),header=F)[,1])
dup <- as.vector(read.table(paste(args[1],"duplication.list",sep=""),header=F)[,1])
inv <- as.vector(read.table(paste(args[1],"inversion.list",sep=""),header=F)[,1])
ins <- as.vector(read.table(paste(args[1],"insertion.list",sep=""),header=F)[,1])
cpx <- as.vector(read.table(paste(args[1],"complex.list",sep=""),header=F)[,1])
unres <- as.vector(read.table(paste(args[1],"unresolved.list",sep=""),header=F)[,1])

#Gather CDF at 0.5% intervals
cdfx <- t(sapply(seq(0,1,by=0.005),function(i){
  return(unlist(lapply(list(del,dup,inv,ins,cpx,unres),function(vals){
    length(vals[which(vals<=i)])/length(vals)
  })))
}))

#Plot CDF
colors <- c("firebrick","dodgerblue","darkorange","darkorchid4","aquamarine","gray25")
pdf(args[2],height=5,width=8)
plot(x=c(0,200),y=c(0,1),type="n",
     main="Structural Variant Frequencies by Class",
     xlab="Maximum Variant Allele Frequency",xaxt="n",
     ylab="Percent of Variants",yaxt="n",xaxs="i",yaxs="i")
rect(xleft=100,xright=par("usr")[2],
     ybottom=par("usr")[3],ytop=par("usr")[4],
     col="gray80")
rect(xleft=c(0,10,20,50),xright=c(10,20,50,100),
     ybottom=par("usr")[3],ytop=par("usr")[4],
     col=c("gray100","gray97","gray94","gray91"),
     border="gray60")
rect(xleft=100,xright=par("usr")[2],
     ybottom=par("usr")[3],ytop=par("usr")[4],
     col="gray60",density=15)
abline(h=seq(0,1,by=0.25),lty=2,col="gray40")
text(x=150,y=0,pos=3,
     labels="VAF > 50%")
for(i in 6:1){
  points(cdfx[,i],type="l",lwd=3,col=colors[i])
}
legend("right",
       legend=c("Deletion",
                "Tandem Duplication",
                "Simple Inversion",
                "Insertion",
                "Complex",
                "Incompletely Resolved Site"),
       col=colors,lwd=3,cex=0.7,
       bg="white")
axis(1,at=c(0,10,20,50,100,200),labels=paste(c(0,5,10,25,50,100),"%",sep=""),cex.axis=0.75)
axis(2,at=seq(0,1,by=0.25),labels=paste(seq(0,100,by=25),"%",sep=""),las=2,cex.axis=0.75)
dev.off()
