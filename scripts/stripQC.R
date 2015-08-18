#!/usr/bin/Rscript

#Positional argument (1): full path to Holmes QC.metrics file (input)
#Positional argument (2): full path to Holmes run summary file (output)

#Read command-line arguments
args <- commandArgs(TRUE)

##Ensures Standard Decimal Notation##
options(scipen=1000)

#Load file
x <- read.table(args[1],comment.char="#",header=F)
colnames(x) <- c("ID","Total Read Pairs","Read Alignment Rate","Pair Alignment Rate","Proper Pairs",
                 "Chimeras","Read Duplicates","Pair Duplicates","Median Insert Size",
                 "Insert MAD","Haploid Physical Coverage","Haploid Nucleotide Coverage","Reported Sex",
                 "Observed Sex","Absolute Dosage Z-Score")
x[,2] <- x[,2]/2000000
x[,3:8] <- 100*x[,3:8]

#Write metrics
for(i in c(2:12,15)){
  iqr <- summary(x[,i])[5]-summary(x[,i])[2]
  if(i==2){
    cat(paste(" => ",colnames(x)[i],"\n      -Median: ",
              prettyNum(summary(x[,i])[3],big.mark=","),
              "M\n      -IQR: ",
              prettyNum(summary(x[,i])[2],big.mark=","),
              "-",
              prettyNum(summary(x[,i])[5],big.mark=","),
              "M\n      -Range: ",
              prettyNum(summary(x[,i])[1],big.mark=","),
              "-",
              prettyNum(summary(x[,i])[6],big.mark=","),
              "M\n      -Outliers: ",
              length(x[which(x[,i]>=median(x[,i])+(1.5*iqr) |
                               x[,i]<=median(x[,i]-1.5*iqr)),i]),
              "\n",sep=""),
        file=args[2],append=T)
  }else if(i %in% c(3:8)){
    cat(paste(" => ",colnames(x)[i],"\n      -Median: ",        
              round(summary(x[,i])[3],2),
              "%\n      -IQR: ",
              round(summary(x[,i])[2],2),
              "-",
              round(summary(x[,i])[5],2),
              "%\n      -Range: ",
              round(summary(x[,i])[1],2),
              "-",
              round(summary(x[,i])[6],2),
              "%\n      -Outliers: ",
              length(x[which(x[,i]>=median(x[,i])+(1.5*iqr) |
                               x[,i]<=median(x[,i]-1.5*iqr)),i]),
              "\n",sep=""),
        file=args[2],append=T)
  }else if(i %in% c(9,10)){
    cat(paste(" => ",colnames(x)[i],"\n      -Median: ",    
              prettyNum(summary(x[,i])[3],big.mark=","),
              "\n      -IQR: ",
              prettyNum(summary(x[,i])[2],big.mark=","),
              "-",
              prettyNum(summary(x[,i])[5],big.mark=","),
              "\n      -Range: ",
              prettyNum(summary(x[,i])[1],big.mark=","),
              "-",
              prettyNum(summary(x[,i])[6],big.mark=","),
              "\n      -Outliers: ",
              length(x[which(x[,i]>=median(x[,i])+(1.5*iqr) |
                               x[,i]<=median(x[,i]-1.5*iqr)),i]),
              "\n",sep=""),
        file=args[2],append=T)
  }else{
    cat(paste(" => ",colnames(x)[i],"\n      -Median: ",    
              round(summary(x[,i])[3],2),
              "\n      -IQR: ",
              round(summary(x[,i])[2],2),
              "-",
              round(summary(x[,i])[5],2),
              "\n      -Range: ",
              round(summary(x[,i])[1],2),
              "-",
              round(summary(x[,i])[6],2),
              "\n      -Outliers: ",
              length(x[which(x[,i]>=median(x[,i])+(1.5*iqr) |
                               x[,i]<=median(x[,i]-1.5*iqr)),i]),
              "\n",sep=""),
        file=args[2],append=T)
  }
}
