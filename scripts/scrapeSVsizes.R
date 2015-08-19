#!/usr/bin/Rscript

#Positional argument (1): full path to sizes directory (subdir of ${WRKDIR})
#Positional argument (2): full path to run_summary.txt (text out)
#Positional argument (3): full path to output pdf (graphics out)

#Reads args
args <- commandArgs(TRUE)

#Sets params
options(scipen=1000)

#Reads data
del <- as.vector(read.table(paste(args[1],"deletion.size",sep=""),header=F)[,1])
dup <- as.vector(read.table(paste(args[1],"duplication.size",sep=""),header=F)[,1])
inv <- as.vector(read.table(paste(args[1],"inversion.size",sep=""),header=F)[,1])
ins_src <- as.vector(read.table(paste(args[1],"insertion_source.size",sep=""),header=F)[,1])
ins_snk <- as.vector(read.table(paste(args[1],"insertion_sink.size",sep=""),header=F)[,1])
cpx <- as.vector(read.table(paste(args[1],"complex.size",sep=""),header=F)[,1])

#Calculate summary
res <- as.data.frame(t(matrix(unlist(lapply(list(del,dup,inv,ins_src,ins_snk,cpx),summary)),nrow=6)))
rownames(res) <- c("Deletion","Duplication","Inversion","Insertion (Source)","Insertion (Sink)","Complex (Approx.)")
colnames(res) <- c("Min","Q1","Med","Mean","Q3","Max")
write.table(res,args[2],sep="\t",append=T,row.names=F,col.names=F,quote=F)

#Plot