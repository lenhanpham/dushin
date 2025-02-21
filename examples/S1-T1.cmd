#!/bin/bash 
#PBS -P xxx      
#PBS -q normal  
#PBS -l walltime=1:00:00 
#PBS -l mem=8GB 
#PBS -l ncpus=2     
#PBS -l jobfs=200GB 
#PBS -l wd
#PBS -l storage=scratch/xxx+gdata/xxx   

currdir=$(pwd) 

export PATH=$PATH:/home/xxx/xxx/xxx/dushin/dushin/bin  


INPUT=${PBS_JOBNAME%.*} 

dushin $INPUT.in > $INPUT.log 2> ${INPUT%.*}.err 


