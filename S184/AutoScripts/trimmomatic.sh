# Autogenerated script from write_trimmomatic_script.R
# date Fri Jun 15 14:10:49 2018
# make sure directory paths exist before running script
#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 1
#$ -cwd
#$ -N trim-array
#$ -j y
#$ -l h_vmem=19G
#$ -t 1-2



module load general/JRE/1.8.0_65
module load bioinformatics/trimmomatic/0.36



# Parse parameter file to get variables.
number=$SGE_TASK_ID
paramfile=trimmomatic_param.txt
 
inr1=`sed -n ${number}p $paramfile | awk '{print $1}'`
inr2=`sed -n ${number}p $paramfile | awk '{print $2}'`
outr1=`sed -n ${number}p $paramfile | awk '{print $3}'`
outr1up=`sed -n ${number}p $paramfile | awk '{print $4}'`
outr2=`sed -n ${number}p $paramfile | awk '{print $5}'`
outr2up=`sed -n ${number}p $paramfile | awk '{print $6}'`

# 9. Run the program.
java -jar /opt/apps/bioinformatics/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33 $inr1 $inr2 $outr1 $outr1up $outr2 $outr2up ILLUMINACLIP:/users/k1625253/brc_scratch/Data/MetaData/trimmomatic_adapters.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36



