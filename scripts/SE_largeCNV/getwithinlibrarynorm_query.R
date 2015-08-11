#!/usr/bin/Rscript

#### This script is written by Serkan Erdin to within-library query samples in July, 2015 #####
#Modified by RC on 8/11/15

args<-commandArgs(TRUE)
fileinput <- args[1]  #insert coverage file for query samples 
outputfilehandler <- args[2]    
outdir <- args[3]

data <- read.table(file=fileinput,head=T,sep="\t")
counts <- data[,-c(1,2,3)]
colMedians <- apply(counts,2,median)
cpm <- t(t(counts)/colMedians)
newdata <- data.frame(data[,c(1:3)],cpm)
colnames(newdata) <- c("data.Chr","data.Start","data.Stop",colnames(newdata)[-c(1:3)])
write.table(newdata,file=paste(outdir,"/",outputfilehandler,".query.bindata.txt",sep=""),sep="\t",row.names=F,quote=F)

		
