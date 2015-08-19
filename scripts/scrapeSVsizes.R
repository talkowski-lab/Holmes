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
cpx <- as.vector(read.table(paste(args[1],"complex.size",sep=""),header=F)[,1])
densities <- lapply(list(del,dup,inv,ins_src,cpx),function(vals){return(density(log10(vals)))})

#Calculate summary
res <- data.frame(c("Deletion","Duplication","Inversion","Insertion (Source)","Complex (Approx.)"),
                  t(matrix(unlist(lapply(list(del,dup,inv,ins_src,cpx),summary)),nrow=6)))
colnames(res) <- c("Class","Min","Q1","Med","Mean","Q3","Max")
write.table(res,args[2],sep="\t",append=T,row.names=F,col.names=T,quote=F)

#Plot
colors <- c("firebrick","dodgerblue","darkorange","darkorchid4","aquamarine")
pdf(args[3],height=6,width=8)
plot(0,0,type="n",
     xlim=c(3,7),ylim=c(0,1.05*max(unlist(lapply(densities,function(list){return(max(list$y))})))),
     lwd=2,col="firebrick",
     xaxs="i",yaxs="i",xaxt="n",
     ylab="Density",xlab="Variant Size",main="Variant Sizes by Class")
grid(nx=NULL,ny=NA,col="black",lty=2)
abline(v=log10(c(seq(2000,9000,by=1000),
                 seq(20000,90000,by=10000),
                 seq(200000,900000,by=100000),
                 seq(2000000,9000000,by=1000000))),
       col="gray20")
for(i in 1:5){
  polygon(x=c(densities[[i]]$x,rev(densities[[i]]$x)),
          y=c(densities[[i]]$y,rep(0,length(densities[[i]]$y))),
          border=colors[i],lwd=3,col=adjustcolor(colors[i],alpha=0.1))
}
axis(1,at=c(3:7),labels=c("1kb","10kb","100kb","1Mb","10Mb"))
legend("topright",
       legend=c("Deletion",
                "Tandem Duplication",
                "Simple Inversion",
                "Insertion",
                "Complex (Approx.)"),
       col=colors,lwd=3,
       bg="white")
dev.off()
