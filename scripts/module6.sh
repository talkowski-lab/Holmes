#!/bin/bash

#liWGS-SV Pipeline: Module 6 (Consensus CNV Merging)
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Create master 