#!/usr/bin/Rscript



args<-commandArgs(TRUE)
fileinput <- args[1]
header <- args[2]
outdir <- args[3]


data <- read.table(file=fileinput,head=T,sep="\t")
data2 <- data[which(data$Log2 != "NA"),]
data2 <- data2[which(data2$Chr != "Y"),]
levels(data2$Chr)[levels(data2$Chr)=="X"] <- "23"

data2[,8]=1:nrow(data2)
data2$scaled = data2$Log2 - median(data2$Log2)
chrOrder <- c(1:23)
aggr<-aggregate(V8~Chr,data2,max)
aggr<-aggr[match(chrOrder,aggr$Chr),]
chr1=0.5*aggr$V8[1]
chr2=aggr$V8[1]+(aggr$V8[2]-aggr$V8[1])/2
chr3=aggr$V8[2]+(aggr$V8[3]-aggr$V8[2])/2
chr4=aggr$V8[3]+(aggr$V8[4]-aggr$V8[3])/2
chr5=aggr$V8[4]+(aggr$V8[5]-aggr$V8[4])/2
chr6=aggr$V8[5]+(aggr$V8[6]-aggr$V8[5])/2
chr7=aggr$V8[6]+(aggr$V8[7]-aggr$V8[6])/2
chr8=aggr$V8[7]+(aggr$V8[8]-aggr$V8[7])/2
chr9=aggr$V8[8]+(aggr$V8[9]-aggr$V8[8])/2
chr10=aggr$V8[9]+(aggr$V8[10]-aggr$V8[9])/2
chr11=aggr$V8[10]+(aggr$V8[11]-aggr$V8[10])/2
chr12=aggr$V8[11]+(aggr$V8[12]-aggr$V8[11])/2
chr13=aggr$V8[12]+(aggr$V8[13]-aggr$V8[12])/2
chr14=aggr$V8[13]+(aggr$V8[14]-aggr$V8[13])/2
chr15=aggr$V8[14]+(aggr$V8[15]-aggr$V8[14])/2
chr16=aggr$V8[15]+(aggr$V8[16]-aggr$V8[15])/2
chr17=aggr$V8[16]+(aggr$V8[17]-aggr$V8[16])/2
chr18=aggr$V8[17]+(aggr$V8[18]-aggr$V8[17])/2
chr19=aggr$V8[18]+(aggr$V8[19]-aggr$V8[18])/2
chr20=aggr$V8[19]+(aggr$V8[20]-aggr$V8[19])/2
chr21=aggr$V8[20]+(aggr$V8[21]-aggr$V8[20])/2
chr22=aggr$V8[21]+(aggr$V8[22]-aggr$V8[21])/2
chrX=aggr$V8[22]+(aggr$V8[23]-aggr$V8[22])/2
#chrY=aggr$V8[23]+(aggr$V8[24]-aggr$V8[23])/2
#y <- lowess(data$V8,data$V8)

curcol=rep("blue4",dim(data2)[1])
curcol[which(data2$Chr==2)]="orange3"
curcol[which(data2$Chr==4)]="orange3"
curcol[which(data2$Chr==6)]="orange3"
curcol[which(data2$Chr==8)]="orange3"
curcol[which(data2$Chr==10)]="orange3"
curcol[which(data2$Chr==12)]="orange3"
curcol[which(data2$Chr==14)]="orange3"
curcol[which(data2$Chr==16)]="orange3"
curcol[which(data2$Chr==18)]="orange3"
curcol[which(data2$Chr==20)]="orange3"
curcol[which(data2$Chr==22)]="orange3"
#curcol[which(data2$Chr==24)]="orange3"

y_sym=filter(data2$scaled,rep(1/201,201,sides=2))


png(paste(outdir,"/",header,"cnv.png",sep=""),width=1200,height=800)

plot(1:nrow(data2),data2$scaled,xaxt="n",bty="n",pch=20,ylab="log2 ratios",xlab="",ylim=c(-4,4),col=curcol,cex=0.08)
axis(side=1,at=c(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX),labels=c(1:22,"X"),tick=F)
abline(h=0,col="black")
abline(h=0.3,col="black")
abline(h=-0.3,col="black")
lines(1:nrow(data2),y_sym,col="red")
dev.off()

fn <- paste(header,".SD.profile",sep="")
if (file.exists(fn)) file.remove(fn)

write.table(paste("Category","SD",sep="\t"),file=paste(header,".SD.profile",sep=""),quote=F,row.names=F,sep="\t",col.names=F)
for(i in 1:23){
data3 <- data2[which(data2$Chr==i),]
sd=sd(data3$scaled)

write.table(data.frame(i,sd),file=paste(header,".SD.profile",sep=""),quote=F,append=T,row.names=F,sep="\t",col.names=F)
rm(data3)
}

sdall=sd(data2$scaled)
write.table(data.frame("ALL",sdall),file=paste(header,".SD.profile",sep=""),quote=F,append=T,row.names=F,sep="\t",col.names=F)

library(DNAcopy)
#data2$scaled = data2$Log2 - median(data2$Log2)
CNA.object <- CNA(genomdat=data2$scaled,chrom=data2$Chr,maploc=data2$Start,data.type='logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segment <- segment(CNA.smoothed,verbose=0,alpha=0.001,nperm=10000,p.method=c("perm"),min.width=5,undo.splits=c("sdundo"),undo.SD=3)

#y_sym=filter(data2$scaled,rep(1/301,301,sides=2))
p.segment <- segments.p(segment)
write.table(p.segment,file=paste(header,".log2.segments.p_value",sep=""),sep="\t",quote=F,col.names=T)

