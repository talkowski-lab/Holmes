#!/usr/bin/Rscript

#---------------------------
# Batch CNV Interval Genotyping
# Talkowski Laboratory
#
# Harrison Brand, Ryan Collins, and Joseph Glessner
# Update October 2015 (load sample set prior)
# Update August 2015 (for incorporation into Holmes liWGS-SV 1.0)
#---------------------------

##POSITIONAL ARGUMENTS##
# [1] Bed file of CNVs to check. No header. ID as fourth column. SampleIDs of interest comma delimited as fifth column.
# [2] Full path to 1kb binned coverage matrix for entire cohort. Generate with binCov.sh if necessary.
# [3] Output file name
# [4] Optional: Verbose TRUE. Default:FALSE
# [5] Optional: outType K. Default:Z
# [6] Optional: prior Default:FALSE
# [7] Optional: delPrior FALSE for Dups Default:TRUE
# [8] Optional: groupC Default:FALSE

#Read command-line arguments
args <- commandArgs(TRUE)

##Ensures Standard Decimal Notation##
options(scipen=1000)

##Loads required packages; installs if necessary##
if("plyr" %in% rownames(installed.packages()) == FALSE)
{install.packages("plyr",repos="http://cran.rstudio.com")}
suppressMessages(library(plyr))
if("MASS" %in% rownames(installed.packages()) == FALSE)
{install.packages("MASS",repos="http://cran.rstudio.com")}
suppressMessages(library(MASS))
if("zoo" %in% rownames(installed.packages()) == FALSE)
{install.packages("zoo",repos="http://cran.rstudio.com")}
suppressMessages(library(zoo))
if("e1071" %in% rownames(installed.packages()) == FALSE)
{install.packages("e1071",repos="http://cran.rstudio.com")}
suppressMessages(library(e1071))
if("mclust" %in% rownames(installed.packages()) == FALSE)
{install.packages("mclust",repos="http://cran.rstudio.com")}
suppressMessages(library(mclust))
if("fpc" %in% rownames(installed.packages()) == FALSE)
{install.packages("fpc",repos="http://cran.rstudio.com")}
suppressMessages(library(fpc))
if("DAAG" %in% rownames(installed.packages()) == FALSE)
{install.packages("DAAG",repos="http://cran.rstudio.com")}
suppressMessages(library(DAAG))
if("reshape" %in% rownames(installed.packages()) == FALSE)
{install.packages("reshape",repos="http://cran.rstudio.com")}
suppressMessages(library(reshape))
if("perm" %in% rownames(installed.packages()) == FALSE)
{install.packages("perm",repos="http://cran.rstudio.com")}
suppressMessages(library(perm))


##Loads regions to plot
intervals <- read.table(args[1],sep="\t",header=F)[,c(1:5)]
intervals[,c(2:3)] <- apply(intervals[,c(2:3)],2,function(vals){return(as.numeric(as.character(vals)))})

#Load plotting values
rawcov <- read.table(args[2],header=T,sep="\t")

#Calculate median bin value per sample
allnorm<-apply(rawcov[,-c(1:3)],2,median)

#Rebinning helper function
rebin <- function(df,compression){
  Chr <- df[1,1]
  Start <- df[1,2]
  End <- df[compression,3]
  for(i in 2:(floor(nrow(df)/compression))){
    Chr <- c(Chr,as.character(df[((i-1)*compression)+1,1]))
    Start <- c(Start,as.integer(df[((i-1)*compression)+1,2]))
    End <- c(End,as.integer(df[i*compression,3]))
  }
  newvals <- apply(df[,4:ncol(df)],2,
                   function(vals,compression){
                     newcol <- sum(vals[1:compression])
                     for(i in 2:(floor(length(vals)/compression))) {
                       newcol <- c(newcol,as.integer(sum(vals[(((i-1)*compression)+1):(i*compression)])))
                     }
                     return(newcol)
                   },compression)
  return(as.data.frame(cbind(Chr,Start,End,newvals)))
}

timeLimit=60              #Time limit for loci since some run long

#Main genotyping function
genotypeCov <- function(chr,start,end,            #region to be genotyped
                        cov=rawcov,               #path to coverage matrix
                        norm=allnorm,             #vector of normalized values per sample
                        ID="Default",      	  #CNV ID
                        type="both",		  #CNV type both or (del,dup)
                        bins=10,		  #How many bins			
                        normDist=30000000,        #distance outside region to normalize
                        tscore=TRUE,              #Include a t-score for calls T of F
                        z.pvaluecutoff=1e-6,      #z test pvalue cutoff (was 1e-9)
			verbose="FALSE",	  #Extra outputs
			outType="Z",		  #"Z"test or "K"means genotyping output type
			IDs="",			  #Sample IDs of interest comma delimited
			prior="FALSE",		  #Prior samples used to be tested by t-score
			delPrior="TRUE",	  #If prior samples are deletions or not (duplications)
			groupC="FALSE"		  #If group C call
){       
  #Checks for appropriate input
  if(!(is.numeric(c(start,end)) & end > start)){
    stop("INPUT ERROR: Improper input coordinates")}

 if(end-start<1000)
 {
	cat("Region too small so inflating 1kb\n")
	start<-start-500
	end<-end+500
 }

 if(is.na(verbose)){verbose="FALSE"} 
 if(is.na(outType)){outType="Z"}
 if(is.na(prior)){prior="FALSE"}
 if(is.na(delPrior)){delPrior="TRUE"}
 if(is.na(groupC)){groupC="FALSE"}

  #Round bins to ensure even counts#
  compression_size=round_any(((round_any(end,1000,floor)-(round_any(start,1000,ceiling)))/bins),1000,floor)
  Rstart<-round_any(start,1000,ceiling)+round_any(((round_any(end,1000,floor))-(round_any(start,1000,ceiling))-(compression_size*bins))/2,1000,floor)
  Rend<-round_any(end,1000,floor)-round_any(((round_any(end,1000,floor))-(round_any(start,1000,ceiling))-(compression_size*bins))/2,1000,ceiling)
  compression=compression_size/1000
  
  #Set compression
  if(Rend-Rstart<=5000){
    Rstart<-start #round_any(start,1000,ceiling)
    Rend<-end #round_any(end,1000,floor)
    compression=1
  }
 
   
  #Approximates rebinned per-sample averages (approximated for speed & memory)
  norm[which(norm==0)]<-1
  norm <- compression*norm
  
  #Subset region of interest from coverage matrix
  cov1 <- cov[which(cov$Start<Rend & cov$Stop>(Rstart) & cov$Chr==chr),]
  #cat(dim(cov1))
  if(dim(cov1)<1){cat(chr,start,end,ID,"No Data Found For Queried Region.\n");return(c(ID,chr,start,end,rep(-9,(dim(cov1)[2]-3))))}
  if(dim(cov1)==1){cat(chr,start,end,ID,"Not Enough (1) Data Found For Queried Region.\n");return(c(ID,chr,start,end,rep(-9,(dim(cov1)[2]-3))))}

CentromereStartsHg19<-c(121500000,90500000,87900000,48200000,46100000,58700000,58000000,43100000,47300000,38000000,51600000,33300000,16300000,16100000,15800000,34600000,22200000,15400000,24400000,25600000,10900000,12200000,58100000,11600000)
CentromereEndsHg19<-c(128900000,96800000,93900000,52700000,50700000,63300000,61700000,48100000,50700000,42300000,55700000,38200000,19500000,19100000,20700000,38600000,25800000,19000000,28600000,29400000,14300000,17900000,63000000,13400000)
chrIndex<-ifelse(chr=="X",23,ifelse(chr=="Y",24,chr))
if(is.na(CentromereEndsHg19[as.numeric(chrIndex)])){cat("ERROR: Unrecognized Chromosome (not 1-22,X,Y)\n");return(c(ID,chr,start,end,rep(-9,(dim(cov1)[2]-3))))}
if((start < CentromereStartsHg19[as.numeric(chrIndex)]) & (end > CentromereEndsHg19[as.numeric(chrIndex)]))
{cat(chr,start,end,ID,"WARNING: Spans Centromere!\n");Centromere=TRUE}
else
{Centromere=FALSE}
 
  #Replace zero values with 1 for handling normalization
  cov1[cov1==0]<-1
#print(cov1) 
#print(compression)
#print(as.numeric(as.matrix(rebin(cov1,compression))[,4]))
  #Rebins values
  if(compression>1 && !(anyNA(as.numeric(as.matrix(rebin(cov1,compression))[,4])))){
    res <- rebin(cov1,compression)
#print(res) 
    res<-apply(res[,4:ncol(res)],2,function(val){as.numeric(as.matrix(val))})  
  }
    else {
    res <- cov1[,4:ncol(cov1)]
  }
#print(res) 

  #Adds sample averages to df
  res1<-rbind(res,t(data.frame(norm)))
  #print(res1)
  #Scale each col within that sample
  res1<- apply(res1,2,
               function(vals){
                   return(as.numeric(vals[1:(nrow(res1)-1)])/as.numeric(vals[nrow(res1)]))
                 })
#print(res1)
  samplenames<-colnames(res1) 
  
  #K means; 0=homozygous deletion, 0.5=het deletion, 1=Normal, 1.5=3 copy, 2=4 copy state
  kmean_matrix<-t(res1)
  #print(cor(kmean_matrix))
#kmean_matrix[is.na(kmean_matrix)]<-1
#kmean_matrix[kmean_matrix==0]<-1
#print(kmean_matrix)
#cat(IDs)
colSds<-c()
for(j in c(1:ncol(kmean_matrix)))
{
	colSds[j]<-sd(kmean_matrix[,j])
}
cat(chr,start,end,ID," median sd all samples and sd:",median(colSds)," ",sd(colSds),"\n")
cat(chr,start,end,ID," [")
HomDelIDs=""
  #HomDel Check
  low=0
  for(i in c(1:nrow(kmean_matrix)))
  {
    if(ncol(kmean_matrix)==1)
    {a<-1;b<-1}
    else if(ncol(kmean_matrix)==2)
    {a<-1;b<-2}
    else
    {a<-2;b<-ncol(kmean_matrix-1)}
	for(j in a:b)
	{
		if(kmean_matrix[i,j]<0.1)
		{
			#cat(kmean_matrix[j,i])
			#cat(samplenames[i])
			low=low+1
		}
	}
	if(low>(ncol(kmean_matrix)-2)/2)
	{
		if(samplenames[i] %in% unlist(strsplit(IDs,split=",")))
		{
			cat(samplenames[i],",",sep="")
			HomDelIDs<-paste(HomDelIDs,samplenames[i],",")
		}
		else
		{
			cat(samplenames[i],"*,",sep="")
		}
	}
	low=0
  }
 cat("]\n")
#print(kmean_matrix)
  kwithinss=c(0,0,0,0,0)
  ks=c(0,0,0,0,0)
  ##Typically 0.5 intervals work best but for higher CN complexity observations, smaller intervals are needed
  for (i in c(1:5)){  
    ##Run Kmeans##
    if(nrow(kmean_matrix)>1){
      k<-kmeans(kmean_matrix,## cmeans
                matrix(rep(seq(0,i/10*19,length.out=20), 
                           ncol(kmean_matrix)),
                       ncol=ncol(kmean_matrix)),
                algorithm="Forgy")
      for (j in 0:20) {k$cluster[k$cluster==(j+1)]<-j}
    }else{k<-kmeans(kmean_matrix,1)
    }
    
    a<-(count(k$centers>0))
    b<-(count(k$withinss>0))
    kwithinss[i]=sum(k$withinss)+(a$freq[1]*1)+(b$freq[2]*3) #Special weighting scheme to minimize within sum of squares but penalize number of clusters to prevent overfitting by using too many clusters
    ks[i]=a$freq[1]
  }
  
  k_1<-kmeans(kmean_matrix,1) #Try 1 cluster all diploid possibility for noisy data
  i=6
  k_1$cluster=k_1$cluster+1
  ks[i]=1
  i=max(which(kwithinss == min(kwithinss))) #pick the highest if tie
  
  if(nrow(kmean_matrix)>1){  
    k<-kmeans(kmean_matrix,matrix(rep(seq(0,i/10*19,length.out=20),ncol(kmean_matrix)),ncol=ncol(kmean_matrix)),algorithm="Forgy")
    for (j in 0:20) {k$cluster[k$cluster==(j+1)]<-j}
  }else{k<-kmeans(kmean_matrix,1)}
  
  ##Calculate cluster metrics to further assess cluster solution from those already provided by kmeans
  d <- dist(kmean_matrix, method="euclidean") #distance matrix
  a<-(count(k$centers>0))
  if(a$freq[1]>1)
  {temp=1;temp$cluster=k$cluster #FIX TO MAINTAIN VALUES In cluster.stats(d, k$cluster) :clustering renumbered because maximum != number of clusters
   ClustSolStats=cluster.stats(d, k$cluster)
   k$cluster=temp$cluster
  }else{ClustSolStats="NA";ClustSolStats$n.between="NA"} #only 1 cluster so cluster.stats does not work
  
  
  ###START Shift scenario checking
  shiftDownOne=0
  ux <- unique(k$cluster)
  if(ux[which.max(tabulate(match(k$cluster, ux)))] == 3)  #MODE
  {maxGroupSize=max(tabulate(match(k$cluster, ux)))
   x<-tabulate(match(k$cluster, ux))
   n <- length(x)
   if(n>1){secondMaxGroupSize=sort(x,partial=n-1)[n-1]}else{secondMaxGroupSize=0}
   ratio=secondMaxGroupSize/maxGroupSize
   if(ratio<0.8) #Significantly more samples in another category suggesting it should be diploid
   {for (j in 0:20) {k$cluster[k$cluster==(j+1)]<-j}; #shiftDownOne=1
   }}
  
  ux <- unique(k$cluster)
  if(min(ux)>2)
  {for(zzz in 1:(min(ux)-2))
  {for (j in 0:20) {k$cluster[k$cluster==(j+1)]<-j}}
  }
  
  zzzz=1
  while((is.na(k$centers[zzzz,1]) | mean(k$centers[zzzz,])==0) && nrow(k$centers)>zzzz)  #Skip over empty clusters
  {zzzz=zzzz+1}
  
  if((is.na(k$centers[zzzz,1]) | mean(k$centers[zzzz,])==0) && nrow(k$centers)>zzzz)
  {if(min(k$centers[zzzz+1,]) < 0.25){if(min(ux)==1){k$cluster[k$cluster==1]<-0}} #Near 0 but called CN1 by clustering fix to CN0
   else{if(min(ux)==0){k$cluster[k$cluster==0]<-1}}  #Near 0.5 but called CN0 by clustering fix to CN1
  }else
  {if(min(k$centers[zzzz,]) < 0.25){if(min(ux)==1){k$cluster[k$cluster==1]<-0}}
   else{if(min(ux)==0){k$cluster[k$cluster==0]<-1}}
  }
  
  k$centers[k$centers=="NaN"]<-0 #cmeans prior cannot have NaNs
  z=1;zz=-9
  for(z in 1:nrow(k$centers))
  {if(k$centers[z]!=0)
  {if(zz==-9)
  {zz=z} #The first row of clusters with non-zero/non-NA values
  }
  }
  k$centers <- unique( k$centers[ , 1:ncol(k$centers) ] ) #cmeans prior cannot have duplicate cluster values
  
  OneNoisyClusterPlusOne=0
  cat(chr,start,end,ID," median separation:",median(ClustSolStats$separation),"\n")
  #if(table(as.numeric(t(genotype[1,c(5:ncol(genotype))])))[1]<ncol(genotype)-5)
  #ClustSolStats=cluster.stats(d, as.numeric(t(genotype[1,c(5:ncol(genotype))])))
  
  if(ClustSolStats$n.between!="NA") #Already 1 cluster solution
  {if(ClustSolStats$n.between/a$freq[1]>3000) #high betweeness suggesting overfitting of clusters
  {a<-(count(k$centers>0))
   if(a$freq[2] == 2)    ### TO DO: control noisy diploid cluster with real cn state cluster/s
   {k<-kmeans(kmean_matrix,1)
    k$cluster=k$cluster+1
    k$centers=rbind(k$centers,-1)  ### cmeans needs at least 2 centers
    i=6
    OneNoisyClusterPlusOne=1
   }
  }
  }
  
  uxMinK <- min(unique(k$cluster))
  #END Shift scenario checking
  if(verbose=="TRUE")
  {
  jpeg(paste(chr,"_",start,"_",end,"_kmean.jpg",sep=""),res=300, width=1800, height=1800)
  plot(as.zoo(res1), plot.type="single",col=k$cluster+1,main=paste("chr",chr,":",start,"-",end," (hg19) k-means",sep=""),ylab="Intensity")
  dev.off()
  }

  ###convert into standard genotype notation## 
  CNV<-matrix(c(ID,chr,start,end),nrow=1)
  colnames(CNV)<-c("ID","Chr","Start","End")
  genotype<-cbind(CNV,t(k$cluster))

  kclusteroutput<-genotype

##Take SD of each cluster
clusterSds<-c()
for(i in 0:20)
{
	for(j in c(1:ncol(kmean_matrix)))
	{
	colSds[j]<-sd(kmean_matrix[which(k$cluster==i),j])
	}
#cat(i,":",median(colSds)," ",sd(colSds),"\n")
clusterSds[i]<-median(colSds)
}
  cat(chr,start,end,ID," median separation k-means solution:",median(ClustSolStats$separation))#print(ClustSolStats)
  cat(" n=",ClustSolStats$n," cluster.number=",ClustSolStats$cluster.number," cluster.size=",ClustSolStats$cluster.size," min.cluster.size=",ClustSolStats$min.cluster.size," noisen=",ClustSolStats$noisen," diameter=",ClustSolStats$diameter," average.distance=",ClustSolStats$average.distance," median.distance=",ClustSolStats$median.distance," separation=",ClustSolStats$separation," average.toother=",ClustSolStats$average.toother," separation.matrix=",ClustSolStats$separation.matrix," ave.between.matrix=",ClustSolStats$ave.between.matrix," average.between=",ClustSolStats$average.between," average.within=",ClustSolStats$average.within," n.between=",ClustSolStats$n.between," n.within=",ClustSolStats$n.within," max.diameter=",ClustSolStats$max.diameter," min.separation=",ClustSolStats$min.separation," within.cluster.ss=",ClustSolStats$within.cluster.ss," clus.avg.silwidths=",ClustSolStats$clus.avg.silwidths," avg.silwidth=",ClustSolStats$avg.silwidth," g2=",ClustSolStats$g2," g3=",ClustSolStats$g3," pearsongamma=",ClustSolStats$pearsongamma," dunn=",ClustSolStats$dunn," dunn2=",ClustSolStats$dunn2," entropy=",ClustSolStats$entropy," wb.ratio=",ClustSolStats$wb.ratio," ch=",ClustSolStats$ch," cwidegap=",ClustSolStats$cwidegap," widestgap=",ClustSolStats$widestgap," sindex=",ClustSolStats$sindex," corrected.rand=",ClustSolStats$corrected.rand," vi=",ClustSolStats$vi," median cluster sds and sd:",median(na.omit(clusterSds)),sd(na.omit(clusterSds)),"\n",sep="_")
cat(chr,start,end,ID," median cluster sds and sd:",median(na.omit(clusterSds)),sd(na.omit(clusterSds)),"\n")
if( exists("ClustSolStats$separation") && median(ClustSolStats$separation)>0.1 && ClustSolStats$avg.silwidth>0.2 && ClustSolStats$pearsongamma>0.3 && ClustSolStats$dunn2>0.5 && ClustSolStats$wb.ratio<0.6 && median(na.omit(clusterSds))<0.2 )
{
	cat(chr,start,end,ID," PASS k-means metrics\n")
}
else
{
	cat(chr,start,end,ID," FAIL k-means metrics\n")
}
  ##cmeans using k-means optimal centers as determined above with little room for adjustment by default, mainly to get the "membership" probability of cluster assignment confidence
  if(nrow(kmean_matrix)>1)
  {b<-cmeans(kmean_matrix,k$centers,method="ufcl")# ,m=2.2,rate.par=1)
  }else{b<-cmeans(rbind(kmean_matrix,1),rbind(k$centers,1),method="ufcl")}
  
  ###START Shift scenario checking
  if(OneNoisyClusterPlusOne==1){
    b$cluster=b$cluster+1
  }else{
    if(zz==2 | zz==1){
      for (j in 0:20) {b$cluster[b$cluster==(j+1)]<-j}}
    if(shiftDownOne==1){  #MODE
      for (j in 0:20) {b$cluster[b$cluster==(j+1)]<-j} }
  }
  
  ux <- unique(b$cluster)
  if(ux[which.max(tabulate(match(b$cluster, ux)))] == 3) #MODE
  {maxGroupSize=max(tabulate(match(b$cluster, ux)))
   x<-tabulate(match(b$cluster, ux))
   n <- length(x)
   if(n>1){secondMaxGroupSize=sort(x,partial=n-1)[n-1]}else{secondMaxGroupSize=0}
   ratio=secondMaxGroupSize/maxGroupSize
   if(ratio<0.8) #Significantly more samples in another category suggesting it should be diploid
   {for (j in 0:20) {b$cluster[b$cluster==(j+1)]<-j}
   }}
  
  ux <- unique(b$cluster)
  
  zzzz=1
  while((is.na(b$centers[zzzz,1]) | mean(b$centers[zzzz,])<0.01) && nrow(b$centers)>zzzz)
  {zzzz=zzzz+1}
  
  if((is.na(b$centers[zzzz,1]) | mean(b$centers[zzzz,])<0.01) && nrow(b$centers)>zzzz){
    if(min(b$centers[zzzz+1,]) < 0.25){if(min(ux)==1){b$cluster[b$cluster==1]<-0}}
    else{if(min(ux)==0){b$cluster[b$cluster==0]<-1}}
  }else{
    if(min(b$centers[zzzz,]) < 0.25){if(min(ux)==1){b$cluster[b$cluster==1]<-0}}
    else{if(min(ux)==0){b$cluster[b$cluster==0]<-1}}
  }
  
  uxMinC <- min(unique(b$cluster))
  if(uxMinK != uxMinC){
    if(uxMinK>uxMinC){for(zzzzz in 1:(uxMinK-uxMinC)){b$cluster=b$cluster+1}}
    if(uxMinK<uxMinC){for(zzzzz in 1:(uxMinC-uxMinK)){b$cluster=b$cluster-1}}
  }
  #END Shift scenario checking
  
  if(verbose=="TRUE")
  {
  write.table(apply(b$membership, 1, max),paste(chr,"_",start,"_",end,"_Probabilities.txt",sep=""),quote=FALSE,col.names=FALSE)
  }

  ###convert into standard genotype notation## 
  CNV<-matrix(c(ID,chr,start,end),nrow=1)
  colnames(CNV)<-c("ID","Chr","Start","End")
  genotype<-cbind(CNV,t(k$cluster))
  
  #write.table(genotype,paste(chr,"_",start,"_",end,"_kclusterTable.txt",sep=""))
  #kclusteroutput<-genotype
  if(prior=="TRUE")
  {
  #cat("Make genotype prior\n")
  for(i in c(5:ncol(genotype)))
  {
	        if(samplenames[i-4] %in% unlist(strsplit(IDs,split=",")))
                {
			if(delPrior=="TRUE")
			{
                        genotype[1,i]=1;
			}
			else
			{
			genotype[1,i]=3;
			}
                }
                else
                {
                        genotype[1,i]=2;
                }
  }
  #print(t(genotype[1,c(5:ncol(genotype))]))
  if(table(as.numeric(t(genotype[1,c(5:ncol(genotype))])))[1]<ncol(genotype)-5)
  {ClustSolStats=cluster.stats(d, as.numeric(t(genotype[1,c(5:ncol(genotype))])))
  cat(chr,start,end,ID," median separation prior solution:",median(ClustSolStats$separation))#print(ClustSolStats)
  cat(" n=",ClustSolStats$n," cluster.number=",ClustSolStats$cluster.number," cluster.size=",ClustSolStats$cluster.size," min.cluster.size=",ClustSolStats$min.cluster.size," noisen=",ClustSolStats$noisen," diameter=",ClustSolStats$diameter," average.distance=",ClustSolStats$average.distance," median.distance=",ClustSolStats$median.distance," separation=",ClustSolStats$separation," average.toother=",ClustSolStats$average.toother," separation.matrix=",ClustSolStats$separation.matrix," ave.between.matrix=",ClustSolStats$ave.between.matrix," average.between=",ClustSolStats$average.between," average.within=",ClustSolStats$average.within," n.between=",ClustSolStats$n.between," n.within=",ClustSolStats$n.within," max.diameter=",ClustSolStats$max.diameter," min.separation=",ClustSolStats$min.separation," within.cluster.ss=",ClustSolStats$within.cluster.ss," clus.avg.silwidths=",ClustSolStats$clus.avg.silwidths," avg.silwidth=",ClustSolStats$avg.silwidth," g2=",ClustSolStats$g2," g3=",ClustSolStats$g3," pearsongamma=",ClustSolStats$pearsongamma," dunn=",ClustSolStats$dunn," dunn2=",ClustSolStats$dunn2," entropy=",ClustSolStats$entropy," wb.ratio=",ClustSolStats$wb.ratio," ch=",ClustSolStats$ch," cwidegap=",ClustSolStats$cwidegap," widestgap=",ClustSolStats$widestgap," sindex=",ClustSolStats$sindex," corrected.rand=",ClustSolStats$corrected.rand," vi=",ClustSolStats$vi,"\n",sep="_")
  }

  ##Take SD of each cluster
  clusterSds<-c()
  for(i in 0:20)
  {
        for(j in c(1:ncol(kmean_matrix)))
        {
        colSds[j]<-sd(kmean_matrix[which(genotype[1,5:ncol(genotype)]==i),j])
        }
  #cat(i,":",median(colSds)," ",sd(colSds),"\n")
  clusterSds[i]<-median(colSds)
  if(i==2) 
  {
	CnTwoMedianSd=median(colSds)
  }
  }
  cat(chr,start,end,ID," median prior cluster sds and sd:",CnTwoMedianSd,sd(na.omit(clusterSds)),"\n")
  if(median(na.omit(clusterSds))>0.1)
  {
	z.pvaluecutoff<-1e-4
  }

  }
  #write.table(genotype,paste(chr,"_",start,"_",end,"_kclusterTablePrior.txt",sep=""))

  ###T-scrore###
  if(tscore==TRUE){
    #Find max geno value
    genotable<-table(genotype[,5:ncol(genotype)])
    highestgeno<-max(genotable)
    genodata<-data.frame(genotable)
    if (length(which(genotable==highestgeno))==1){
      genoMode<-names(which(genotable==highestgeno))
    }else{   ###if more than one max###
      genoMode<-names(which(genotable==highestgeno))
      genotable1<-data.frame(genotable)
      if (length(which(c(genotable1[genoMode,1])==2))>0) { genoMode<-2 } 
      else if (length(which(c(genotable1[genoMode,1])==3))>0) { genoMode<-3 } 
      else if (length(which(c(genotable1[genoMode,1])==1))>0){ genoMode<-1 }
      else { genoMode<-4 }
    }
  }
  
  if(nrow(kmean_matrix)>1){ 
    names<-colnames(genotype)[5:ncol(genotype)]
    genoids<-names[which(genotype[,5:ncol(genotype)]==genoMode)]
    genoids0<-names[which(genotype[,5:ncol(genotype)]==0)]
    genoids1<-names[which(genotype[,5:ncol(genotype)]==1)]
    genoids2<-names[which(genotype[,5:ncol(genotype)]==2)]
    genoids3<-names[which(genotype[,5:ncol(genotype)]==3)]
    genoids4<-names[which(genotype[,5:ncol(genotype)]==4)]
    genoids5<-names[which(genotype[,5:ncol(genotype)]==5)]
    genoids6<-names[which(genotype[,5:ncol(genotype)]==6)]
    genoids7<-names[which(genotype[,5:ncol(genotype)]==7)]
    genoids8<-names[which(genotype[,5:ncol(genotype)]==8)]
    genoids9<-names[which(genotype[,5:ncol(genotype)]==9)]
    genoids10<-names[which(genotype[,5:ncol(genotype)]==10)]
    genoids11<-names[which(genotype[,5:ncol(genotype)]==11)]
    genoids12<-names[which(genotype[,5:ncol(genotype)]==12)]
    genoids13<-names[which(genotype[,5:ncol(genotype)]==13)]
    genoids14<-names[which(genotype[,5:ncol(genotype)]==14)]
    genoids15<-names[which(genotype[,5:ncol(genotype)]==15)]
    genoids16<-names[which(genotype[,5:ncol(genotype)]==16)]
    genoids17<-names[which(genotype[,5:ncol(genotype)]==17)]
    genoids18<-names[which(genotype[,5:ncol(genotype)]==18)]
    genoids19<-names[which(genotype[,5:ncol(genotype)]==19)]
    genoids20<-names[which(genotype[,5:ncol(genotype)]==20)]
    if(length(genoids)>1){
      kmeanwgeno<-apply(kmean_matrix[genoids,],1,mean)
      }else{
        kmeanwgeno<-mean(kmean_matrix[genoids,])
      }
    cat(chr,start,end,ID," ")
    if(prior==TRUE && length(kmean_matrix[which(genotype[,5:ncol(genotype)]==(ifelse(delPrior=="TRUE",1,3)))])>1 && length(kmean_matrix[which(genotype[,5:ncol(genotype)]==2)])>1 && !anyNA(cor(kmean_matrix)))
    {
	oneSample=FALSE
	t.power1 = function(nsamp=c(10,10),nsim=1000,means=c(1,0.5),sds=c(1,1),var.equal=TRUE)
	{
   	tps = replicate(nsim,
       	t.test(rnorm(nsamp[1],mean=means[1],sd=sds[1]),
              rnorm(nsamp[2],mean=means[2],sd=sds[2]))$p.value)

    	sum(tps < .025 | tps > .975) / nsim
	}
    Ttests<-c()
    power<-c()
    PermTTests<-c()
    cat("ncol",ncol(kmean_matrix)," ")
    if(ncol(kmean_matrix)==2)
    {a<-1;b<-2}
    else
    {a<-2;b<-ncol(kmean_matrix-1)}
    for(i in a:b)
    { 
     Control<-kmean_matrix[which(genotype[,5:ncol(genotype)]==2),i]
     #cat("mean:");cat(mean(Control));cat("sd:");cat(sd(Control));cat("\n");cat(Control);cat("\n")
     Treat<-kmean_matrix[which(genotype[,5:ncol(genotype)]==(ifelse(delPrior=="TRUE",1,3))),i]
     if(sd(Treat)==0)
     {
        Treat[1]=Treat[1]+.00001
     }
     #cat("mean:");cat(mean(Treat));cat("sd:");cat(sd(Treat));cat("\n");cat(Treat);cat("\n")
     if(delPrior=="TRUE"){Ttest<-t.test(Control,Treat,alternative="greater")}
     else{Ttest<-t.test(Control,Treat,alternative="less")}
     #cat(Ttest$p.value)
     Ttests[i-1]<-Ttest$p.value
     #cat("\nlens:",length(Control)," ",length(Treat),"\n")
     #cat("\npower:")
     power[i-1]<-t.power1(sds=c(sd(Control),sd(Treat)),nsamp=c(length(Control),length(Treat)),var.equal=FALSE)
     #cat(power[i-1])
     #cat("\n")
     if(delPrior=="TRUE"){myPermObject<-permTS(Control,Treat,alternative="greater",method='pclt');PermTTests[i-1]<-capture.output(cat(as.numeric(myPermObject$p.value)))}
     else{myPermObject<-permTS(Control,Treat,alternative="less",method='pclt');PermTTests[i-1]<-capture.output(cat(as.numeric(myPermObject$p.value)))}
    #cat("c=",c,"\n")
     #PermTTests[i-1]<-as.numeric(capture.output(cat(as.character(colsplit(c,split=" ",names="d")[2]))))
     #cat(PermTTests[i-1])
   } 
    cat("median p:",median(Ttests)," median power:",median(power)," median perm:",median(as.numeric(PermTTests[!is.na(PermTTests)]))," ")
    cat("separation=",median(ClustSolStats$separation)," avg.silwidth=",ClustSolStats$avg.silwidth," pearsongamma=",ClustSolStats$pearsongamma," dunn=",ClustSolStats$dunn," dunn2=",ClustSolStats$dunn2," wb.ratio=",ClustSolStats$wb.ratio," ")
    if(median(power)>=0.8 && groupC==FALSE)
    {
	cat("\n")
    }
     
    }
    else{ 
     if(anyNA(cor(kmean_matrix)))
     {
        cat("Correlation in Matrix so cannot Perform t-test. "); oneSample=FALSE; power=-1 }
     else{ 
	cat("Cannot run t-test on one sample "); oneSample=TRUE }}
      ###Z-Test to compare between individual samples and maximum sample group### 
      for (i in 0:20){ 
        temp<-paste("genoids",i,sep="")
        if (length(eval(parse(text = temp)))>0 && genoMode>i ){
          ztest<-function(dataset,mu) {1-pnorm((mean(dataset)-mu)/sd(dataset))} 
        }else if (length(eval(parse(text = temp)))>0 && genoMode<i){
          ztest<-function(dataset,mu) {pnorm((mean(dataset)-mu)/sd(dataset))}
        }
        if (length(eval(parse(text = temp)))>0 && genoMode!=i){
          kmeannogeno<-apply(kmean_matrix[setdiff(eval(parse(text = temp)),genoids),,drop=FALSE],1,mean)
          tscore<-matrix(apply(as.matrix(kmeannogeno),2,function(x){ztest(kmeanwgeno,x)}))
	  tscore1<-matrix(paste(ID,tscore))
          row.names(tscore)<-names(kmeannogeno)
          row.names(tscore1)<-names(kmeannogeno)
          fail<-names(tscore[which(tscore[,1]>z.pvaluecutoff),])
          pass<-names(tscore[which(tscore[,1]<=z.pvaluecutoff),])
	  genotype1<-t(genotype)
          genotype1[fail,]<-genoMode
          genotype<-t(genotype1)
      }
    }
          if(oneSample==TRUE || median(power)<0.8 || groupC==TRUE)
          {
		j=1;
		#cat("IN THE PRINT IDS BLOCK! ")
                for(i in c(5:ncol(genotype)))
                {
		#cat(samplenames[i-4])
                if(samplenames[i-4] %in% unlist(strsplit(IDs,split=",")))
                {
                        cat(genotype[1,i],":",samplenames[i-4],"#",tscore[j],"/",sep="")
			j=j+1
                }

                }
                cat("\n")
          }

  }

  #Start state shifting scenarios
  z$cluster<-genotype[1,5:ncol(genotype)]
  ux <- unique(z$cluster)
  if(ux[which.max(tabulate(match(z$cluster, ux)))] == 3)  #MODE
  {maxGroupSize=max(tabulate(match(z$cluster, ux)))
   x<-tabulate(match(z$cluster, ux))
   n <- length(x)
   if(n>1){secondMaxGroupSize=sort(x,partial=n-1)[n-1]}else{secondMaxGroupSize=0}
   ratio=secondMaxGroupSize/maxGroupSize
   if(ratio<0.8) #Significantly more samples in another category suggesting it should be diploid
   {for (j in 0:20) {z$cluster[z$cluster==(j+1)]<-j}; #shiftDownOne=1
   }}
  else if(ux[which.max(tabulate(match(z$cluster, ux)))] == 4)  #MODE
  {maxGroupSize=max(tabulate(match(z$cluster, ux)))
   x<-tabulate(match(z$cluster, ux))
   n <- length(x)
   if(n>1){secondMaxGroupSize=sort(x,partial=n-1)[n-1]}else{secondMaxGroupSize=0}
   ratio=secondMaxGroupSize/maxGroupSize
   if(ratio<0.8) #Significantly more samples in another category suggesting it should be diploid
   {for (j in 0:20) {z$cluster[z$cluster==(j+2)]<-j}; #shiftDownOne=1
   }}
  else if(ux[which.max(tabulate(match(z$cluster, ux)))] == 5)  #MODE
  {maxGroupSize=max(tabulate(match(z$cluster, ux)))
   x<-tabulate(match(z$cluster, ux))
   n <- length(x)
   if(n>1){secondMaxGroupSize=sort(x,partial=n-1)[n-1]}else{secondMaxGroupSize=0}
   ratio=secondMaxGroupSize/maxGroupSize
   if(ratio<0.8) #Significantly more samples in another category suggesting it should be diploid
   {for (j in 0:20) {z$cluster[z$cluster==(j+3)]<-j}; #shiftDownOne=1
   }}

  ###convert into standard genotype notation## 
  CNV<-matrix(c(ID,chr,start,end),nrow=1)
  colnames(CNV)<-c("ID","Chr","Start","End")
  genotype<-cbind(CNV,t(z$cluster))

  #End state shifting scenarios
  
  if(verbose=="TRUE")
  {
  jpeg(paste(chr,"_",start,"_",end,"_ztest.jpg",sep=""),res=300, width=1800, height=1800)
  plot(as.zoo(res1), plot.type="single",col=(as.numeric(genotype[1,5:ncol(genotype)])+1),main=paste("chr",chr,":",start,"-",end," (hg19) genotype",sep=""),ylab="Intensity")
  
  ###Mark samples of interest
  if(nrow(kmean_matrix)>1)
  {
  samplesPrior <- unlist(strsplit(IDs,","))
  for(i in 1:nrow(res1))
  {
        for(j in 1:ncol(res1))
    {
      for(l in 1:length(samplesPrior))
      {
                if(samplesPrior[l] == samplenames[j])
                {text(i, res1[i,j], "X",cex=0.5) }###samplenames[j],cex=0.5) }
          }
        }
  }
  }

  dev.off()
  }
	if(prior==TRUE && groupC==FALSE)
	{
	if(Centromere==TRUE)
        {
                for(i in c(5:ncol(genotype)))
                {
                        genotype[1,i]=2;
                }
        }
	if(exists("PermTTests")) ###"ClustSolStats")) ###"Ttests"))
	{
        if(median(as.numeric(PermTTests[!is.na(PermTTests)]))>=0.01 && median(power)>=0.8) ###ClustSolStats$dunn2<=1.5) ###median(Ttests)>=0.05 && median(power)>=0.8) minimum average dissimilarity between two cluster / maximum average within cluster dissimilarity, another version of the family of Dunn indexes See Halkidi et al. Clustering Validity Checking Methods: Part II
        {
                for(i in c(5:ncol(genotype)))
                {
                        genotype[1,i]=2;
                }
        }
	else if(median(as.numeric(PermTTests[!is.na(PermTTests)]))<0.01 && median(power)>=0.8) ###ClustSolStats$dunn2>1.5) ###median(Ttests)<0.05 && median(power)>=0.8)
	{
	  for(i in c(5:ncol(genotype)))
  	  {
                if(samplenames[i-4] %in% unlist(strsplit(IDs,split=",")))
                {
                        if(delPrior=="TRUE")
                        {
                        if(samplenames[i-4] %in% unlist(strsplit(HomDelIDs,split=",")))
                        {
                        genotype[1,i]=0;
                        }
			else
			{
                        genotype[1,i]=1;
                        }
			}
                        else
                        {
                        genotype[1,i]=3;
                        }
                }
                else
                {
                        genotype[1,i]=2;
                }
  	  }
	}
	}
	}
  #Returns Genotype
  if(outType=="Z")
  {
  return(genotype)
  }
  else
  {
  return(kclusteroutput)
  }
  
}

#Genotypes intervals
#cat(dim(rawcov)[2]-3)
timeOut <- function (expr, cpu = Inf, elapsed = Inf)
{
setTimeLimit(cpu = cpu, elapsed = elapsed, transient = TRUE)
on.exit(setTimeLimit()) # should not be needed, but averts future error message
expr
}
results_matrix <- as.data.frame(t(apply(intervals,1,function(row){ tryCatch(timeOut({return(suppressWarnings(genotypeCov(chr=row[1],
                                                                                                      start=as.numeric(as.character(row[2])),
                                                                                                      end=as.numeric(as.character(row[3])),
                                                                                                      ID=row[4], verbose=args[4], outType=args[5], IDs=row[5], prior=args[6], delPrior=args[7], groupC=args[8])))}, elapsed=timeLimit), error = function(e) {cat("Error: Took Too Long!\n");c(row[4],row[1],as.numeric(as.character(row[2])),end=as.numeric(as.character(row[3])),rep(-9,(dim(rawcov)[2]-3)))})
})))

colnames(results_matrix) <- c("CNV_ID","chr","start","end",colnames(rawcov[,-c(1:3)]))
results_matrix[,-c(1:4)] <- t(apply(results_matrix[,-c(1:4)],1,function(vals){vals <- as.numeric(as.character(vals)); return(vals)})) # return(vals-median(vals)+2)}))

#Writes genotypes to file
write.table(results_matrix,args[3],sep="\t",col.names=T,row.names=F,quote=F)

