#!/bin/bash

##MASTER liWGS-SV PIPELINE SCRIPT
#Contact: rcollins@chgr.mgh.harvard.edu

#Submit this script to long queue on ERISOne to run entire pipeline

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Set up output tree

