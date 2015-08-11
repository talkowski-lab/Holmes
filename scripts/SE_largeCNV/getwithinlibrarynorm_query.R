#!/usr/bin/Rscript

#### This script is written by Serkan Erdin to within-library query samples in July, 2015 #####


args<-commandArgs(TRUE)
fileinput <- args[1]  #insert coverage file for query samples 
outputfilehandler <- args[2]    
outdir <- args[3]
option <- args[4] #it has two options: nozerorows and keep. The former removes bins with zero insert coverage from analysis, while the latter keeps all bins


if(option == "nonzerorows"){
	data <- read.table(file=fileinput,head=T,sep="\t")
	data$rowsums <- rowSums(data[,-c(1,2,3)])
	data <- data[which(data$rowsums > 0),]
	dim <- dim(data)
	counts <- data[,-c(1,2,3,dim[2])]
	colMedians <- apply(counts,2,median)
	cpm <- t(t(counts)/colMedians)
	newdata <- data.frame(data$Chr,data$Start,data$End,cpm)
	write.table(newdata,file=paste(outdir,"/",outputfilehandler,".query.bindata.txt",sep=""),sep="\t",row.names=F,quote=F)
}else if(option == "keep"){
	 data <- read.table(file=fileinput,head=T,sep="\t")
	counts <- data[,-c(1,2,3)]
        colMedians <- apply(counts,2,median)
        cpm <- t(t(counts)/colMedians)
	newdata <- data.frame(data$Chr,data$Start,data$End,cpm)
	write.table(newdata,file=paste(outdir,"/",outputfilehandler,".query.bindata.txt",sep=""),sep="\t",row.names=F,quote=F)
}	
		
