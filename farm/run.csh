#!/bin/tcsh 

set ENV_FILE    = "/u/home/dmriser/analysis-projects-clas12/environment/farm.csh"
set GROOVY_FILE = "/u/home/dmriser/analysis-projects-clas12/projects/elastic/monitor-rad.groovy"
#set OUTPUT_DIR  = "/volatile/clas12/dmriser/farm_out/elastic_jan24/data/03/"
#set DATA_DIR    = "/work/clas12/rg-a/trains/v16_v2/skim8_ep"
#set DATA_DIR = "/lustre/expphy/volatile/clas12/rg-a/production/recon/pass0/calib/v35/recon/005038/"

#mkdir -p $OUTPUT_DIR
source $ENV_FILE
#run-groovy $GROOVY_FILE $DATA_DIR/* 
#cp *.hipo $OUTPUT_DIR

# One job per node
run-groovy $GROOVY_FILE input.hipo 

