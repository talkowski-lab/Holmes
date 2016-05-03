#! /usr/bin/env Rscript

#Rscript to run cn.MOPS after genome-wide insert coverage has been created

#####CAN BE PARALLELIZED LIKE THIS######
# cov=[[COV FILE GOES HERE]]
# group=[[groupID GOES HERE]]
# list=[[SAMPLE LIST GOES HERE]]
# mkdir ./cnMOPS_chrsplit
# for contig in $( seq 1 22 ) X Y; do
#  mkdir ./cnMOPS_chrsplit/${contig}
#  cat <( head -n1 ${cov} ) <( awk -v chr=chr${contig} '{ if ($1==chr) print }' ${cov} ) > ./cnMOPS_chrsplit/${contig}/${group}.rawCov.chr${contig}.bed
#  for binsize in 1 3 10 30; do
#   bsub -q big -M 20000 -sla miket_sc -R 'rusage[mem=20000]' -M 20000 -v 26000 -J ${group}_cnMOPS_${binsize}_chr${contig} "Rscript /data/talkowski/rlc47/code/SV/cnMOPS_postcoverage.R -m insert -r ${binsize} -b ${binsize}000 -I ${group} ./cnMOPS_chrsplit/${contig}/${group}.rawCov.chr${contig}.bed ./cnMOPS_chrsplit/${contig}/"
#  done
# done
# GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${group}_cnMOPS_" | wc -l)
#   GATEwait=0
#   until [[ $GATEcount == 0 ]]; do
#    sleep 20s
#    GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${group}_cnMOPS_" | wc -l)
#    GATEwait=$[$GATEwait +1]
#    if [[ $GATEwait == 15 ]]; then
#     echo "$(date): INCOMPLETE"
#     echo "$(date): Waiting on ${GATEcount} jobs to complete"
#     GATEwait=0
#    fi
#   done
# for contig in $( seq 1 22 ) X Y; do
#  for binsize in 1 3 10 30; do
#   fgrep -v "#" ./cnMOPS_chrsplit/${contig}/${group}.${binsize}000kbBins.cnMOPS.gff | awk -v OFS="\t" '{ print $1, $4, $5, $9, $10, $11, $12 }' | sed 's/^chr//g' >> ./cnMOPS_chrsplit/${group}.cnMOPS_master.cnMOPS.gff
#  done
# done
# mkdir ./cnMOPS_calls
# while read bam ID; do
#   ID=$( echo ${ID} | sed 's/\-/\./g' )
#   echo ${ID}
#   mkdir ./cnMOPS_calls/${ID}
#   fgrep -w "${ID}" ./cnMOPS_chrsplit/${group}.cnMOPS_master.cnMOPS.gff | grep 'CN[0-1]' | sed 's/median\=//g' | sed 's/mean\=//g' | sed 's/CN\=//g' | sed 's/\;//g' > ./cnMOPS_chrsplit/${ID}.cnMOPS.preMerge.dels.bed
#   fgrep -w "${ID}" ./cnMOPS_chrsplit/${group}.cnMOPS_master.cnMOPS.gff | grep 'CN[3-9]' | sed 's/median\=//g' | sed 's/mean\=//g' | sed 's/CN\=//g' | sed 's/\;//g' > ./cnMOPS_chrsplit/${ID}.cnMOPS.preMerge.dups.bed
#   bedtools merge -d 1 -c 4,5,6,7 -o distinct,mean,mean,distinct -i <( sed -e 's/^X/23/g' -e 's/^Y/24/g' ./cnMOPS_chrsplit/${ID}.cnMOPS.preMerge.dels.bed | sort -nk1,1 -k2,2 | sed -e 's/^23/X/g' -e 's/^24/Y/g' ) > ./cnMOPS_calls/${ID}/${ID}.cnMOPS.dels.bed
#   bedtools merge -d 1 -c 4,5,6,7 -o distinct,mean,mean,distinct -i <( sed -e 's/^X/23/g' -e 's/^Y/24/g' ./cnMOPS_chrsplit/${ID}.cnMOPS.preMerge.dups.bed | sort -nk1,1 -k2,2 | sed -e 's/^23/X/g' -e 's/^24/Y/g' ) > ./cnMOPS_calls/${ID}/${ID}.cnMOPS.dups.bed
# done < ${list}
# rm -rf ./cnMOPS_chrsplit


#################
#RScript Input
#################

#installs optparse if not already installed
if("optparse" %in% rownames(installed.packages()) == FALSE)
{install.packages("optparse",repos="http://cran.rstudio.com")}
suppressWarnings(suppressPackageStartupMessages(library(optparse)))

#checks for cn.MOPS installation
if("cn.mops" %in% rownames(installed.packages()) == FALSE)
{warning("ERROR: cn.MOPS package not installed.  Attempting to install from Bioconductor...")
source("http://bioconductor.org/biocLite.R")
biocLite("cn.mops")}
suppressPackageStartupMessages(library(cn.mops))

#checks for rtracklayer installation
if("rtracklayer" %in% rownames(installed.packages()) == FALSE)
{warning("ERROR: rtracklayer package not installed.  Attempting to install from Bioconductor...")
source("http://bioconductor.org/biocLite.R")
biocLite("rtracklayer")}
suppressPackageStartupMessages(library(rtracklayer))

#checks for BioRC installation
if("BioRC" %in% rownames(installed.packages()) == FALSE)
{stop("FATAL ERROR: BioRC package not properly installed.  Please install package to current library path manually.")}
suppressPackageStartupMessages(library(BioRC))

#list of command-line options
option_list <- list(
  make_option(c("-m", "--mode"), type="character", default="physical",
              help="specify coverage mode for evaluation from either 'physical' or 'nucleotide' [default %default]", 
              metavar="character"),
  make_option(c("-b","--binsize"), type="integer", default=1000,
              help="bin size for coverage evaluation; this corresponds to one-third of the minimum call resolution [default %default]",
              metavar="integer"),
  make_option(c("-r","--rebin"), type="integer", default=1,
              help="dictates bin compression factor pre-cnMOPS [default %default]",
              metavar="integer"),
  make_option(c("-I", "--ID"), type="character", default="CNMOPS_unknown",
              help="sample group ID [default %default]", 
              metavar="character"),
  make_option(c("-o", "--output_format"), type="character", default="multisample",
              help="output format; choose from either 'list' or 'multisample' [default %default]", 
              metavar="character")
 )

#Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog [options] covMatrix.bed OUTDIR", option_list=option_list),positional_arguments=TRUE)
opts <- args$options

#checks for appropriate positional arguments
if(length(args$args) != 2) 
  {cat("Incorrect number of required positional arguments\n\n")
   stop()}

#writes args & opts to vars
path.to.matrix <- args$args[1]
OUTDIR <- args$args[2]
mode <- opts$mode
binsize <- opts$binsize
ID <- opts$ID
out_format <- opts$output_format
rebinsize <- opts$rebin

#Rebin FX
rebin <- function(df,compression){
  Chr <- df[1,1]
  Start <- df[1,2]
  End <- df[compression,3]
  for(i in 2:(floor(nrow(df)/compression))) {
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

#Loads coverage matrix & coerces to GRanges
if(rebinsize>1){
  cov <- rebin(read.table(path.to.matrix,header=T,sep="\t",stringsAsFactors=F),rebinsize)
} else {
  cov <- read.table(path.to.matrix,header=T,sep="\t",stringsAsFactors=F)
}
cov[2:ncol(cov)] <- apply(cov[2:ncol(cov)],2,function(vals){return(as.numeric(as.character(vals)))})
colnames(cov)[2:3] <- c("Start","Stop")
cov[is.na(cov)] <- 0
cov <- dataframe2GRanges(cov)

#Runs cn.mops
res <- cn.mops(cov)
res <- calcIntegerCopyNumbers(res)

#Exports CNVs
{export(cnvs(res), paste(OUTDIR,"/",ID,".",binsize,"bpBins.cnMOPS.gff",sep=""), "GFF3")}
