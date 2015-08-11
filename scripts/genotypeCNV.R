#---------------------------
# Genotype CNV interval by insert depth
# Talkowski Laboratory
#
# Harrison Brand, Ryan Collins, and Joseph Glessner
# Update May 15 2015
#---------------------------

#First Argument should be list (txt file) of CNVs to check in the following format chr,bp1,bp2,CNVID ##
###Second Argument vector specifying cohorts or samples for normalization ###
###Third Argument scaling of each column to mean or median
##Fourth Argument output type, either genotype or kmeans
##FIFTH Argument CNV type, recommend "both" which looks for both del and dup, but can also specify del or dup###
##Sixth Argument, How many bins to normalize across in the deletion### 
##Seventh File name for output###

#export PATH=/PHShome/jtg24/R-3.1.2/bin:$PATH
#1_54797114_54797344_1

chr<-1
start<-54797114#177264817#15793314#50776000#24280142
end<-54797344#177272879#15817097#50928000#24286682
batch<-1
normMid="median"
ID="Default"
output="genotype"
type="both"
bins=10   
normDist=30000000
tscore="T"
tscoreoutput="tscore.txt"
commonallele=NULL
z.pvaluecutoff=1e-9


args <- commandArgs(TRUE)

samplesBlacklisted <- c("SFARI_p6D12","SFARI_p6F05","SFARI_p4E11","SFARI_p4D06","SFARI_15939815","SFARI_p4F01","SFARI_15935029","SFARI_p6D04","SFARI_p4E09","SFARI_p6F04","SFARI_p6G02","SFARI_p6E05","SFARI_p4F02","SFARI_p6E03","SFARI_p6E08","SFARI_p4D12","SFARI_p4D07","SFARI_p4F10","SFARI_15937460","SFARI_p6C12","SFARI_15936142","SFARI_16042092","SFARI_p6C04","SFARI_15939811","SFARI_12591p1","SFARI_12591fa","SFARI_12591mo","SFARI_p6F02","SFARI_p6E07","SFARI_p6E06","SFARI_p6E04","SFARI_p4H01","SFARI_p4F05","SFARI_16042095","SFARI_15937371","SFARI_15936975","SFARI_15936895","SFARI_15936763","SFARI_15936760","SFARI_15936751","SFARI_15936745","SFARI_15936181","SFARI_15936159","SFARI_15936152","SFARI_15936746","SFARI_p4F08","SFARI_p4F09","SFARI_p6D01","SFARI_JSB02","SFARI_JSB04","SFARI_JSB10","SFARI_JSpD06","SFARI_JSpD10","SFARI_JSpE02","SFARI_JSp7A12")

genotypeCov <- function(chr,start,end,            #region to be genotyped
                        batch,                    #SFARI batch (1,2=Eco, 3=JS)
                        normMid="median",         #scaling of each column to mean or median
                        ID="Default",    		      #CNV ID
                        output="genotype",			  #output type, either genotype or kmeans, kmeans will not work on smaller samples 
                        type="both",				    	#CNV type both or (del,dup)
                        bins=10,			              #How many bins			
                        normDist=30000000,        #distance outside region to normalize
                        tscore="T",               #Include a t-score for calls T of F
                        tscoreoutput="tscore.txt",#tscore ouput file to use
                        commonallele=NULL,      #List providing assigned common alleles (variant name, kmean allele(0,1,2,3,4)) No header!
                        z.pvaluecutoff=1e-9      #z test pvalue cutoff
){       
  
  #dbDisconnect(TGDB)
  ##Checks for appropriate input##
  if(!(is.numeric(c(start,end)) & end > start)){
    stop("INPUT ERROR: Improper input coordinates")}
  if(!(normMid %in% c("mean","median"))){
    stop('INPUT ERROR: Specify "mean" or "median" for normMetric')}
  
  cat(paste(chr,"_",start,"_",end,"_",batch,"\n",sep=""))
  write.table(paste(chr,"_",start,"_",end,"_",batch,".txt",sep=""),paste(chr,"_",start,"_",end,"_",batch,"_RangesRan.txt",sep="")) ## To keep track of Ranges Ran since output not always generated
  
  ##Ensures Standard Decimal Notation##
  options(scipen=400)
  
  ##Loads required packages; installs if necessary##
  if("plyr" %in% rownames(installed.packages()) == FALSE)
  {install.packages("plyr",repos="http://cran.rstudio.com")}
  library(plyr)
  if("RMySQL" %in% rownames(installed.packages()) == FALSE)
  {install.packages("RMySQL",repos="http://cran.rstudio.com/")}
  library(RMySQL)
  if("MASS" %in% rownames(installed.packages()) == FALSE)
  {install.packages("MASS",repos="http://cran.rstudio.com")}
  library(MASS)
  if("zoo" %in% rownames(installed.packages()) == FALSE)
  {install.packages("zoo",repos="http://cran.rstudio.com")}
  library(zoo)
  if("HardyWeinberg" %in% rownames(installed.packages()) == FALSE)
  {install.packages("HardyWeinberg",repos="http://cran.rstudio.com")}
  library(HardyWeinberg)
  if("e1071" %in% rownames(installed.packages()) == FALSE)
  {install.packages("e1071",repos="http://cran.rstudio.com")}
  library(e1071)
  #library(mclust)
  if("fpc" %in% rownames(installed.packages()) == FALSE)
  {install.packages("fpc",repos="http://cran.rstudio.com")}
  library(fpc)
  ###Round bins to ensure even counts###
  compression_size=round_any(((round_any(end,1000,floor)-(round_any(start,1000,ceiling)))/bins),1000,floor)
  Rstart<-round_any(start,1000,ceiling)+round_any(((round_any(end,1000,floor))-(round_any(start,1000,ceiling))-(compression_size*bins))/2,1000,floor)
  Rend<-round_any(end,1000,floor)-round_any(((round_any(end,1000,floor))-(round_any(start,1000,ceiling))-(compression_size*bins))/2,1000,ceiling)
  compression=compression_size/1000
  
  ##Set compression##
  if(Rend-Rstart<5000){
    Rstart<-start
    Rend<-end
    compression=1
  }
  
  ##Connects to TGDB##
  cat("Attempting to connect to TGDB...")
  TGDB <- dbConnect(MySQL(),
                    user="talkLab",
                    password="Talkow_1",
                    dbname='TalkowskiGenomicsDB',
                    host='mysql2.dipr.partners.org')
  cat(" Success!\n")
  
  ##Table assignment##
  if(batch==1){
    table <- "SFARI_B1_FragCount1kbBins"
  }else if(batch==2){
    table <- "SFARI_B2_FragCount1kbBins"
  }else{
    table <- "SFARI_JS_FragCount1kbBins"
  }
  
  ##Define Columns to Load##
  samples <- c(unlist(sapply(c("SFARI"),function(cohort){
    return(dbGetQuery(TGDB,paste("DESCRIBE ",table,sep=""))$Field[grep("SFARI",
                                                                       dbGetQuery(TGDB,
                                                                                  paste("DESCRIBE SFARI_B",batch,"_FragCount1kbBins",sep=""))$Field,
                                                                       ignore.case=T)])}),use.names=F))
  
  samples <- sapply(samples,function(val){
  if(!val %in% samplesBlacklisted)
  {
    return(paste("`",val,"`",sep=""))
  }
  })
  samples <- samples[ ! sapply(samples, is.null) ]
  colsToLoad <- paste(c("`Chr`","`Start`","`End`",samples),collapse=", ")
  
  ##Load Plotting Values##
  cat("Querying coverage database... ")
  cov <- dbGetQuery(TGDB,
                    paste("SELECT ",colsToLoad," ",
                          "FROM ",table," ",
                          "WHERE `Chr` = '",chr,"' ",
                          "AND `Start` <= ",end+normDist," ",
                          "AND `End` >= ",start-normDist,
                          sep=""))
  cat("Complete!\n")
  
  ###Replace zero values with 1 for handling normalization###
  cov[cov==0]<-1
  
  ##Rebin FX##
  rebin <- function(df,compression){
    Chr <- df[1,1]
    Start <- df[1,2]
    End <- df[compression,3]
    for(i in 2:(floor(nrow(df)/compression))    ) {
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
  
  ##Rebins values##
  if(compression>1){
    res <- rebin(cov,compression) 
    res<-apply(res[,4:ncol(res)],2,function(val){as.numeric(as.matrix(val))})	
  } else {
    res <- cov[,4:ncol(cov)]
  }
  ###Take Average###
  if(normMid=="mean"){
    norm<-apply(res,2,mean)
  } else if(normMid=="median") {
    norm<-apply(res,2,median)
  }
  
  cov1 <- dbGetQuery(TGDB,
                     paste("SELECT ",colsToLoad," ",
                           "FROM ",table," ",
                           "WHERE `Chr` = '",chr,"' ",
                           "AND `Start` <= ",Rend-1," ",
                           "AND `End` >= ",Rstart,
                           sep=""))	
    cat("Complete2!\n")
  ###replace zero values with 1 for handling normalization###	
  cov1[cov1==0]<-1
  
  if(compression>1){
    res<-rebin(cov1,compression)    
    res<-apply(res[,4:ncol(res)],2,function(val){as.numeric(as.matrix(val))})	
  } else {
    res <- cov1[,4:ncol(cov1)]
  }
  res1<-rbind(res,t(data.frame(norm)))
  
  ##Scale each col within that sample##
  res1<- apply(res1,2,
               function(vals){
                 if(normMid=="mean"){
                   return(as.numeric(vals[1:(nrow(res1)-1)])/as.numeric(vals[nrow(res1)]))
                 } else if(normMid=="median") {
                   return(as.numeric(vals[1:(nrow(res1)-1)])/as.numeric(vals[nrow(res1)]))
                 }})
  samplenames<-colnames(res1) 
  
  ###K means###
  ## 0<-homozygous deletion, .5<-heterozygous deletion, 1<-Normal, 1.5<-3 copy, 2<-4 copy state ##
  kmean_matrix<-t(res1)
 
 kwithinss=c(0,0,0,0,0)#0,0,0,0,0)
 ks=c(0,0,0,0,0)#0,0,0,0,0)
##Typically 0.5 intervals work best but for higher CN complexity observations, smaller intervals are needed
  for (i in c(1:5)) 
  {  
  ##Run Kmeans##
  if (type=="both") {
   if(nrow(kmean_matrix)>1)
   {
    k<-kmeans(kmean_matrix,## cmeans
			  matrix(rep(seq(0,i/10*19,length.out=20), 
			  #matrix(rep(seq(0.5,4.4,0.2), ### chr15:21800000-22200000 sum(k$withinss)  5.759784 6.028274 - 0.0000000 0.1357046 0.7911720 1.2588087 1.1911400 1.2637990 0.4699501 0.9176996
			  #matrix(rep(seq(0,5.7,0.3),   # sum(k$withinss 18.98771
			  #matrix(rep(seq(0,7.6,0.4),    # 23.92638
			  #matrix(rep(seq(0,9.5,0.5), ### chr15:22200000-22650000 sum(k$withinss) 28.89642 - 0.6147357 11.9911207  9.4842575  6.8063038
			  ncol(kmean_matrix)),
              ncol=ncol(kmean_matrix)),
               algorithm="Forgy") ## method="ufcl")
	for (j in 0:20) {k$cluster[k$cluster==(j+1)]<-j}
   }else{k<-kmeans(kmean_matrix,1)}
  }else if(type=="del"){
    k<-kmeans(kmean_matrix,
              matrix(rep(c(0,.5,1),ncol(kmean_matrix)),
                     ncol=ncol(kmean_matrix)),
              algorithm="Forgy")
    k$cluster[k$cluster==1]<-0
    k$cluster[k$cluster==2]<-1
    k$cluster[k$cluster==3]<-2
  }else if(type=="dup"){
    k<-kmeans(kmean_matrix,
              matrix(rep(c(1,1.5,2),ncol(kmean_matrix)),
                     ncol=ncol(kmean_matrix)),
              algorithm="Forgy")
    k$cluster[k$cluster==3]<-4
    k$cluster[k$cluster==2]<-3
    k$cluster[k$cluster==1]<-2
  }   
  
   cat(sum(k$withinss))##cat(k$withinerror) ##
   cat("\n")
   a<-(count(k$centers>0)) ## a<-(count(k$size>0)) ## 
   cat(a$freq[1])  ## cat(a$freq[1])  ###     Number of clusters
   cat("\n") 
   b<-(count(k$withinss>0))
   cat(b$freq[2])  ## cat(a$freq[1])  ###    Number of non-one off clusters
   cat("\n") 
  kwithinss[i]=sum(k$withinss)+(a$freq[1]*1)+(b$freq[2]*3) ###(a$freq[1]*1)+(b$freq[2]*5)  ### Special weighting scheme to minimize within sum of squares but penalize number of clusters to prevent overfitting by using too many clusters
  ks[i]=a$freq[1]
 }
 k_1<-kmeans(kmean_matrix,1)   ### Try 1 cluster all diploid possibility for noisy data
 i=6
 k_1$cluster=k_1$cluster+1
 #kwithinss[i]=(sum(k_1$withinss))+1+3
 ks[i]=1
 cat(sum(k_1$withinss));cat("\n1\n1\n")
 
 cat(kwithinss)
  cat("\n")
  
  i=max(which(kwithinss == min(kwithinss))) ##pick the highest if tie ####which.min(kwithinss) 
 if(nrow(kmean_matrix)>1){  
  k<-kmeans(kmean_matrix,matrix(rep(seq(0,i/10*19,length.out=20),ncol(kmean_matrix)),ncol=ncol(kmean_matrix)),algorithm="Forgy") ## method="ufcl")##
  for (j in 0:20) {k$cluster[k$cluster==(j+1)]<-j}
 }else{k<-kmeans(kmean_matrix,1)}
  ##Calculate cluster metrics to further assess cluster solution from those already provided by kmeans
  d <- dist(kmean_matrix, method="euclidean") # distance matrix
  a<-(count(k$centers>0))
  if(a$freq[1]>1)
  {temp=1;temp$cluster=k$cluster ###FIX TO MAINTAIN VALUES In cluster.stats(d, k$cluster) :clustering renumbered because maximum != number of clusters
  #subtractOne=0;if(a$freq[1] < max(k$cluster)){subtractOne=1} 
  ClustSolStats=cluster.stats(d, k$cluster)
  cat("NBetween="); cat(ClustSolStats$n.between); cat("\n")  ### Very high values for degree of overlapping clusters
  write.table(ClustSolStats$n.between,paste("ClustSolStats_NBetween.txt",sep=""),quote=FALSE,append=TRUE,col.names=FALSE,row.names=FALSE)
  k$cluster=temp$cluster
  #if(subtractOne==1){k$cluster=k$cluster-1}
  }
  else{ClustSolStats="NA";ClustSolStats$n.between="NA"  ##only 1 cluster so cluster.stats does not work
  write.table("NA",paste("ClustSolStats_NBetween.txt",sep=""),quote=FALSE,append=TRUE,col.names=FALSE,row.names=FALSE)}  
  
  write.table(k$cluster,paste(chr,"_",start,"_",end,"_",batch,".txt",sep=""))
  
  
  jpeg(paste(chr,"_",start,"_",end,"_",batch,".jpg",sep=""),res=300, width=1800, height=1800)
  par(mfrow=c(3,2))
  #par(mar=c(4,4,4,4))
  #par(mfrow=c(2,1))
  
  ###START Shift scenario checking
  shiftDownOne=0
  ux <- unique(k$cluster)
  if(ux[which.max(tabulate(match(k$cluster, ux)))] == 3)  ### MODE
  {
  maxGroupSize=max(tabulate(match(k$cluster, ux)))
  x<-tabulate(match(k$cluster, ux))
  n <- length(x)
  if(n>1){secondMaxGroupSize=sort(x,partial=n-1)[n-1]}else{secondMaxGroupSize=0}
  ratio=secondMaxGroupSize/maxGroupSize
  if(ratio<0.8)   ### Significantly more samples in another category suggesting it should be diploid
  {
  for (j in 0:20) {k$cluster[k$cluster==(j+1)]<-j};#shiftDownOne=1
  }}
  
  ux <- unique(k$cluster)
  #if(ux[which.max(tabulate(match(k$cluster, ux)))] > 3)  ### MODE really high, no diploids problem
  if(min(ux)>2)
  {
	for(zzz in 1:(min(ux)-2))
	{for (j in 0:20) {k$cluster[k$cluster==(j+1)]<-j}}
  }
  
  zzzz=1
  while((is.na(k$centers[zzzz,1]) | mean(k$centers[zzzz,])==0) && nrow(k$centers)>zzzz)  ##Skip over empty clusters
  {zzzz=zzzz+1}
  
  if((is.na(k$centers[zzzz,1]) | mean(k$centers[zzzz,])==0) && nrow(k$centers)>zzzz)
  {
	if(min(k$centers[zzzz+1,]) < 0.25){if(min(ux)==1){k$cluster[k$cluster==1]<-0;cat("1->0")}}   #Near 0 but called CN1 by clustering fix to CN0
	else{if(min(ux)==0){k$cluster[k$cluster==0]<-1;cat("0->1")}}     #Near 0.5 but called CN0 by clustering fix to CN1
  }else
  {
	if(min(k$centers[zzzz,]) < 0.25){if(min(ux)==1){k$cluster[k$cluster==1]<-0;cat("1->0")}}
	else{if(min(ux)==0){k$cluster[k$cluster==0]<-1;cat("0->1")}}
  }
  
  k$centers[k$centers=="NaN"]<-0   ##cmeans prior cannot have NaNs
  z=1;zz=-9
  for(z in 1:nrow(k$centers))
  {
	if(k$centers[z]!=0)
	{
		if(zz==-9)
		{zz=z}   ### The first row of clusters with non-zero/non-NA values
	}
  }
  cat("zz=");cat(zz);cat("\n")
  k$centers <- unique( k$centers[ , 1:ncol(k$centers) ] ) ##cmeans prior cannot have duplicate cluster values
  
  #row_sub = apply(k$centers, 1, function(row) all(row !=0 ))  ## Remove 0 rows
  #k$centers<-k$centers[row_sub,]
  #row.names(k$centers)=c(1,2,3,4,5)
  
  OneNoisyClusterPlusOne=0
  #if(mean(k$centers[nrow(k$centers),]-k$centers[nrow(k$centers)-1,])<0.3) 
  if(ClustSolStats$n.between!="NA")###Already 1 cluster solution
  {
  if(ClustSolStats$n.between/a$freq[1]>3000) ###high betweeness suggesting overfitting of clusters
  {
	a<-(count(k$centers>0))
	if(a$freq[2] == 2)    ### TO DO: control noisy diploid cluster with real cn state cluster/s
	{
		k<-kmeans(kmean_matrix,1)
		k$cluster=k$cluster+1
		k$centers=rbind(k$centers,-1)  ### cmeans needs at least 2 centers
		i=6
		OneNoisyClusterPlusOne=1
	}
  }
  }
  
  uxMinK <- min(unique(k$cluster))
  ###END Shift scenario checking
  
  
    plot(as.zoo(res1), plot.type="single",col=k$cluster+1,main=paste("chr",chr,":",start,"-",end," (hg19) k-means",sep=""),ylab="Intensity") ######1:ncol(res2)


  ###Model based clustering did not work well and had error when co-linearity in data so removed 
  #######k$centers=k$centers[!(apply(k$centers, 1, function(y) any(y == 0))),]
  #a<-Mclust(kmean_matrix,prior=priorControl(mean=k$centers),G=ks[i])
  #a$class=a$class+1
  #plot(as.zoo(res1), plot.type="single",col=a$class+1,main=paste("Mclust class",sep=" "),ylab="Intensity") ######1:ncol(res2)
  #######plot(as.zoo(res1), plot.type="single",col=a$uncertainty*10+1,main=paste("Mclust uncertainty",sep=" "),ylab="Intensity") ######1:ncol(res2)
  #write.table(a$class,paste(chr,"_",start,"_",end,"_",batch,"_Class.txt",sep=""),quote=FALSE,col.names=FALSE)
  #write.table(a$uncertainty,paste(chr,"_",start,"_",end,"_",batch,"_Uncertainty.txt",sep=""),quote=FALSE,col.names=FALSE)

  ##cmeans using k-means optimal centers as determined above with little room for adjustment by default, mainly to get the "membership" probability of cluster assignment confidence
 if(nrow(kmean_matrix)>1)
 {
  b<-cmeans(kmean_matrix,k$centers,method="ufcl")# ,m=2.2,rate.par=1)
 }else{b<-cmeans(rbind(kmean_matrix,1),rbind(k$centers,1),method="ufcl")}
  #b$cluster=b$cluster+1
  
  ###START Shift scenario checking
  if(OneNoisyClusterPlusOne==1)
  {
  #  cat("shiftDownOne\n");for (j in 0:20) {b$cluster[b$cluster==(j+1)]<-j}
	cat("shiftUpOne\n");b$cluster=b$cluster+1
  }
  else
  {if(zz==2 | zz==1)
  {cat("shiftDownOneByzz\n");for (j in 0:20) {b$cluster[b$cluster==(j+1)]<-j} }
  if(shiftDownOne==1)  ### MODE
  {cat("shiftDownOne\n");for (j in 0:20) {b$cluster[b$cluster==(j+1)]<-j} }
  }
  
  ux <- unique(b$cluster)
  if(ux[which.max(tabulate(match(b$cluster, ux)))] == 3)  ### MODE
  {
  maxGroupSize=max(tabulate(match(b$cluster, ux)))
  x<-tabulate(match(b$cluster, ux))
  n <- length(x)
  if(n>1){secondMaxGroupSize=sort(x,partial=n-1)[n-1]}else{secondMaxGroupSize=0}
  ratio=secondMaxGroupSize/maxGroupSize
  if(ratio<0.8)   ### Significantly more samples in another category suggesting it should be diploid
  {
  for (j in 0:20) {b$cluster[b$cluster==(j+1)]<-j};#shiftDownOne=1
  }}
  
  ux <- unique(b$cluster)
	
  zzzz=1
  while((is.na(b$centers[zzzz,1]) | mean(b$centers[zzzz,])<0.01) && nrow(b$centers)>zzzz)
  {zzzz=zzzz+1}
  
  if((is.na(b$centers[zzzz,1]) | mean(b$centers[zzzz,])<0.01) && nrow(b$centers)>zzzz)
  {
	if(min(b$centers[zzzz+1,]) < 0.25){if(min(ux)==1){b$cluster[b$cluster==1]<-0;cat("1->0")}}
	else{if(min(ux)==0){b$cluster[b$cluster==0]<-1;cat("0->1")}}
  }else
  {
	if(min(b$centers[zzzz,]) < 0.25){if(min(ux)==1){b$cluster[b$cluster==1]<-0;cat("1->0")}}
	else{if(min(ux)==0){b$cluster[b$cluster==0]<-1;cat("0->1")}}
  }
  
  uxMinC <- min(unique(b$cluster))
  if(uxMinK != uxMinC)
  {cat("DifferentMinBetweenKandC\n")
  if(uxMinK>uxMinC){for(zzzzz in 1:(uxMinK-uxMinC)){b$cluster=b$cluster+1}}
  if(uxMinK<uxMinC){for(zzzzz in 1:(uxMinC-uxMinK)){b$cluster=b$cluster-1}}
  }
  ###END Shift scenario checking
  plot(as.zoo(res1), plot.type="single",col=b$cluster+1,main=paste("cmeans cluster",sep=" "),ylab="Intensity") ######1:ncol(res2)
  plot(as.zoo(res1), plot.type="single",col=apply(b$membership, 1, max)*20,main=paste("cmeans membership",sep=" "),ylab="Intensity") ######1:ncol(res2)
  write.table(apply(b$membership, 1, max),paste(chr,"_",start,"_",end,"_",batch,"_Membership.txt",sep=""),quote=FALSE,col.names=FALSE)
  
  write.table(b$membership,paste(chr,"_",start,"_",end,"_",batch,"_MembershipAll.txt",sep=""),quote=FALSE,col.names=FALSE)
  write.table(b$cluster,paste(chr,"_",start,"_",end,"_",batch,"_cmeansClusterGenotype.txt",sep=""),quote=FALSE,col.names=FALSE)

  ###Mark samples of interest
 if(nrow(kmean_matrix)>1)
 {
  samplesPrior=c("SFARI_p6C08")  ## TO DO: load as a file
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
 
  plot(sort(apply(b$membership, 1, max)),ylab="Probability",main="cmeans membership")

  #for(i in 1:nrow(res1)){ text(i, res1[i,]+((rank(res1[i,])/100)^3)*3, samplenames,cex=0.5) }
  
  #hist(as.zoo(res1),prob=1)
  #lines(density(as.zoo(res1)),col="red")
  
  
  ###convert into standard genotype notation## 
  CNV<-matrix(c(ID,chr,start,end),nrow=1)
  colnames(CNV)<-c("ID","Chr","Start","End")
  genotype<-cbind(CNV,t(k$cluster))
  ##Output##
  ###T-scrore###
  if(tscore=="T") 
  {
    ###Find max geno value###	
    genotable<-table(genotype[,5:ncol(genotype)])
    highestgeno<-max(genotable)
    genodata<-data.frame(genotable)
    if (length(which(genotable==highestgeno))==1)
	{
      genoMode<-names(which(genotable==highestgeno))
    } 
	else 
	{   ###if more than one max###
      genoMode<-names(which(genotable==highestgeno))
      genotable1<-data.frame(genotable)
      if (length(which(c(genotable1[genoMode,1])==2))>0) { genoMode<-2 } 
	  else if (length(which(c(genotable1[genoMode,1])==3))>0) { genoMode<-3 } 
	  else if (length(which(c(genotable1[genoMode,1])==1))>0){ genoMode<-1 }
	  else { genoMode<-4 }
    }
  }
  cat(genoMode,file=paste(chr,"_",start,"_",end,"_",batch,"_genoMode.txt",sep=""))
  
  ##Print common allele if overwritten by param##
  if (!is.null(commonallele)){
    commonA<-read.table(commonallele,row.names=1)
    genoMode<-commonA[ID,1]
    cat(genoMode,file=paste(chr,"_",start,"_",end,"_",batch,"_genoMode.txt",sep=""))
  }
  
 if(nrow(kmean_matrix)>1)
 { 
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
  for (i in 0:20) 
  { 
    temp<-paste("genoids",i,sep="")
    if (length(eval(parse(text = temp)))>0 && genoMode>i ) 
	{
      ztest<-function(dataset,mu) {1-pnorm((mean(dataset)-mu)/sd(dataset))} 
	} 
	else if (length(eval(parse(text = temp)))>0 && genoMode<i) 
	{
      ztest<-function(dataset,mu) {pnorm((mean(dataset)-mu)/sd(dataset))}
    }
	if (length(eval(parse(text = temp)))>0 && genoMode!=i)
	{
	  kmeannogeno<-apply(kmean_matrix[setdiff(eval(parse(text = temp)),genoids),,drop=FALSE],1,mean)
      tscore<-matrix(apply(as.matrix(kmeannogeno),2,function(x){ztest(kmeanwgeno,x)}))
      tscore1<-matrix(paste(ID,tscore))
      row.names(tscore)<-names(kmeannogeno)
      row.names(tscore1)<-names(kmeannogeno)
      write.table(tscore1,paste(chr,"_",start,"_",end,"_",batch,"_",tscoreoutput,sep=""),quote=FALSE,append=TRUE,col.names=FALSE)
      fail<-names(tscore[which(tscore[,1]>z.pvaluecutoff),])
      genotype1<-t(genotype)
      genotype1[fail,]<-genoMode
      genotype<-t(genotype1)
	}
  }
  
  plot(as.zoo(res1), plot.type="single",col=(as.numeric(genotype[1,5:ncol(genotype)])+1),main=paste("chr",chr,":",start,"-",end," (hg19) genotype",sep=""),ylab="Intensity") ######1:ncol(res2)
  
  dbDisconnect(TGDB)
  
  write.table(t(genotype),paste(chr,"_",start,"_",end,"_",batch,"_genotype.txt",sep=""))
}  
  #### Hardy Weinberg Equilibrium p-value testing for further cluster solution validation
  HW.test.ps=rep(-.2,6)
  x <- c(MM = sum(k$cluster==0), MN = sum(k$cluster==1), NN = sum(k$cluster==2))
  if(sum(x)>0){HW.test <- HWChisq(x); HW.test.ps[1]=HW.test$pval; write.table(paste(chr,"_",start,"_",end,"_",batch,"_kclusterDel=",HW.test$pval,sep=""),paste(chr,"_",start,"_",end,"_",batch,"_HWE_P.txt",sep="")) }
  x <- c(MM = sum(k$cluster==4), MN = sum(k$cluster==3), NN = sum(k$cluster==2))
  if(sum(x)>0){HW.test <- HWChisq(x); HW.test.ps[2]=HW.test$pval; write.table(paste(chr,"_",start,"_",end,"_",batch,"_kclusterDup=",HW.test$pval,sep=""),paste(chr,"_",start,"_",end,"_",batch,"_HWE_P.txt",sep=""),append=TRUE) }

  #x <- c(MM = sum(a$class==0), MN = sum(a$class==1), NN = sum(a$class==2))
  #if(sum(x)>0){HW.test <- HWChisq(x); HW.test.ps[2]=HW.test$pval; write.table(paste(chr,"_",start,"_",end,"_",batch,"_MclusterDel=",HW.test$pval,sep=""),paste(chr,"_",start,"_",end,"_",batch,"_HWE_P.txt",sep=""),append=TRUE) }
  #x <- c(MM = sum(a$class==4), MN = sum(a$class==3), NN = sum(a$class==2))
  #if(sum(x)>0){HW.test <- HWChisq(x); HW.test.ps[3]=HW.test$pval; write.table(paste(chr,"_",start,"_",end,"_",batch,"_MclusterDup=",HW.test$pval,sep=""),paste(chr,"_",start,"_",end,"_",batch,"_HWE_P.txt",sep=""),append=TRUE) }
  
  x <- c(MM = sum(b$cluster==0), MN = sum(b$cluster==1), NN = sum(b$cluster==2))
  if(sum(x)>0){HW.test <- HWChisq(x); HW.test.ps[3]=HW.test$pval; write.table(paste(chr,"_",start,"_",end,"_",batch,"_cMeansclusterDel=",HW.test$pval,sep=""),paste(chr,"_",start,"_",end,"_",batch,"_HWE_P.txt",sep=""),append=TRUE) }
  x <- c(MM = sum(b$cluster==4), MN = sum(b$cluster==3), NN = sum(b$cluster==2))
  if(sum(x)>0){HW.test <- HWChisq(x); HW.test.ps[4]=HW.test$pval; write.table(paste(chr,"_",start,"_",end,"_",batch,"_cMeansclusterDup=",HW.test$pval,sep=""),paste(chr,"_",start,"_",end,"_",batch,"_HWE_P.txt",sep=""),append=TRUE) }
    
  x <- c(MM = sum(genotype==0), MN = sum(genotype==1), NN = sum(genotype==2))
  if(sum(x)>0){HW.test <- HWChisq(x); HW.test.ps[5]=HW.test$pval; write.table(paste(chr,"_",start,"_",end,"_",batch,"_genotypeDel=",HW.test$pval,sep=""),paste(chr,"_",start,"_",end,"_",batch,"_HWE_P.txt",sep=""),append=TRUE) }
  x <- c(MM = sum(genotype==4), MN = sum(genotype==3), NN = sum(genotype==2))
  if(sum(x)>0){HW.test <- HWChisq(x); HW.test.ps[6]=HW.test$pval; write.table(paste(chr,"_",start,"_",end,"_",batch,"_genotypeDup=",HW.test$pval,sep=""),paste(chr,"_",start,"_",end,"_",batch,"_HWE_P.txt",sep=""),append=TRUE) }
  x=c("kmeansDel","kmeansDup","cmeansDel","cmeansDup","genoDel","genoDup")
  check=count(is.na(-log(HW.test.ps))|is.infinite(-log(HW.test.ps)))
  if(check$freq == 6 && check$x == TRUE)
  {}else{
  plot(-log(HW.test.ps),xaxt='n',xlab="",main="Hardy Weinberg Equilibrium")
  axis(1, at = seq(length(x)), labels = x, las=2,cex.axis=.7)
  }
  
  dev.off()
  
  #if(output=="genotype"){
  #  return(genotype)}
  #else if(output=="kmeans"){
  #	return(k)}
}

coord<-read.table(args[1])  
samplenames<-colnames(genotypeCov(2,24280142,24286682,normCohorts=args[2],tscore="F"))

geno<-apply(coord,1,function(vals){genotypeCov(chr=vals[1],start=as.numeric(vals[2]),end=as.numeric(vals[3]),ID=vals[4],normCohorts=args[2],
                                               normMid=args[3],output=args[4],type=args[5],bins=as.numeric(args[6]),tscoreoutput=args[8],commonallele=args[9],z.pvaluecutoff=args[10])})  
rownames(geno)<-samplenames
write.table(t(geno),args[7],quote=FALSE,row.names=FALSE)