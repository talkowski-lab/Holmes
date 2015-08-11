#!/usr/bin/Rscript



args<-commandArgs(TRUE)
fileinput <- args[1]
header <- args[2]
outdir <- args[3]
#lengthThreshold <- args[4]

chromlength=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566,16569)

data <- read.table(file=fileinput,head=T,sep="\t")
data2 <- data[which(data$Log2 != "NA"),]
data2 <- data2[which(data2$Chr != "Y"),]
levels(data2$Chr)[levels(data2$Chr)=="X"] <- "23"
data2$scaled = data2$Log2 - median(data2$Log2)

movingav=filter(data2$scaled,rep(1/201,201,sides=2))
#movingmed=runmed(data2$scaled,201)
gfcoeff <- function(n,s){
t <- seq(-n,n,1)
return(exp(-(t^2/(2*s^2)))/sqrt(2*pi*s^2))
}
gaussianfilter=filter(data2$scaled,gfcoeff((2*100+1)/6,100))

for(i in 1:23){
chr=data2[which(data2$Chr == i),]
movingmed=runmed(chr$scaled,201)
png(paste(outdir,"/",header,"cnv.chr",i,".png",sep=""),width=800,height=600)
plot(chr$Start,chr$scaled,col="gray90",xlim=c(1,chromlength[i]),pch=16,cex=0.75,ylab="Log 2 ratio",xlab="Bp",main=paste(header,":Chr",i,sep=""),cex.axis=1.2,cex.lab=1.2)
lines(chr$Start,movingmed,col="black",lwd=3)
abline(h=0)
dev.off()
rm(chr)
}
