#!/usr/bin/Rscript

#---------------------------
# Batch CNV Interval Genotyping
# Talkowski Laboratory
#
# Harrison Brand, Ryan Collins, and Joseph Glessner
# Update August 2015 (for incorporation into Holmes liWGS-SV 1.0)
#---------------------------

##POSITIONAL ARGUMENTS##
# [1] Bed file of CNVs to check. No header. ID as fourth column.
# [2] Full path to 1kb binned coverage matrix for entire cohort. Generate with binCov.sh if necessary.
# [3] Output file name

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
{install.packages("fpc",repos="http://cran.rstudio.com")}
suppressMessages(library(mclust))
if("fpc" %in% rownames(installed.packages()) == FALSE)
{install.packages("fpc",repos="http://cran.rstudio.com")}
suppressMessages(library(fpc))

##Loads regions to plot
intervals <- read.table(args[1],sep="\t",header=F)[,c(1:4)]
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

#Main genotyping function
genotypeCov <- function(chr,start,end,            #region to be genotyped
                        cov=rawcov,                  #path to coverage matrix
                        norm=allnorm,                #vector of normalized values per sample
                        ID="Default",      	      #CNV ID
                        type="both",				    	#CNV type both or (del,dup)
                        bins=10,			            #How many bins			
                        normDist=30000000,        #distance outside region to normalize
                        tscore=TRUE,              #Include a t-score for calls T of F
                        z.pvaluecutoff=1e-9      #z test pvalue cutoff
){       
  #Checks for appropriate input
  if(!(is.numeric(c(start,end)) & end > start)){
    stop("INPUT ERROR: Improper input coordinates")}
  
  #Round bins to ensure even counts#
  compression_size=round_any(((round_any(end,1000,floor)-(round_any(start,1000,ceiling)))/bins),1000,floor)
  Rstart<-round_any(start,1000,ceiling)+round_any(((round_any(end,1000,floor))-(round_any(start,1000,ceiling))-(compression_size*bins))/2,1000,floor)
  Rend<-round_any(end,1000,floor)-round_any(((round_any(end,1000,floor))-(round_any(start,1000,ceiling))-(compression_size*bins))/2,1000,ceiling)
  compression=compression_size/1000
  
  #Set compression
  if(Rend-Rstart<5000){
    Rstart<-start
    Rend<-end
    compression=1
  }
  
  #Approximates rebinned per-sample averages (approximated for speed & memory)
  norm <- compression*norm
  
  #Subset region of interest from coverage matrix
  cov1 <- cov[which(cov$Start<=Rend & cov$Stop>=(Rstart-1) & cov$Chr==chr),]
  
  #Replace zero values with 1 for handling normalization
  cov1[cov1==0]<-1
  
  #Rebins values
  if(compression>1){
    res <- rebin(cov1,compression) 
    res<-apply(res[,4:ncol(res)],2,function(val){as.numeric(as.matrix(val))})  
  } else {
    res <- cov1[,4:ncol(cov1)]
  }
  
  #Adds sample averages to df
  res1<-rbind(res,t(data.frame(norm)))
  
  #Scale each col within that sample
  res1<- apply(res1,2,
               function(vals){
                   return(as.numeric(vals[1:(nrow(res1)-1)])/as.numeric(vals[nrow(res1)]))
                 })
  samplenames<-colnames(res1) 
  
  #K means; 0=homozygous deletion, 0.5=het deletion, 1=Normal, 1.5=3 copy, 2=4 copy state
  kmean_matrix<-t(res1)
  
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
  
  ###convert into standard genotype notation## 
  CNV<-matrix(c(ID,chr,start,end),nrow=1)
  colnames(CNV)<-c("ID","Chr","Start","End")
  genotype<-cbind(CNV,t(k$cluster))
  
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
    kmeanwgeno<-apply(kmean_matrix[genoids,],1,mean)
    
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
        genotype1<-t(genotype)
        genotype1[fail,]<-genoMode
        genotype<-t(genotype1)
      }
    }
  }
  
  #Returns Genotype
  return(genotype)
  
}

#Genotypes intervals
results_matrix <- as.data.frame(t(apply(intervals,1,function(row){return(suppressWarnings(genotypeCov(chr=row[1],
                                                                                                      start=as.numeric(as.character(row[2])),
                                                                                                      end=as.numeric(as.character(row[3])),
                                                                                                      ID=row[4])))})))
colnames(results_matrix) <- c("CNV_ID","chr","start","end",colnames(rawcov[,-c(1:3)]))

#Writes genotypes to file
write.table(args[3],sep="\t",col.names=T,row.names=F,quote=F)