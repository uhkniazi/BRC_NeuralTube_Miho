# Autogenerated script from 08_methylation_extractor_array_job.R
# date Wed Aug 23 15:15:56 2017
# make sure directory paths exist before running script
#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 8
#$ -cwd
#$ -N meth-ex-array
#$ -j y
#$ -l h_vmem=19G
#$ -t 1-8



module load bioinformatics/samtools/1.3.1



# Parse parameter file to get variables.
number=$SGE_TASK_ID
paramfile=methylation_extractor_param.txt
 
inr1=`sed -n ${number}p $paramfile | awk '{print $1}'`

# 9. Run the program.
/users/k1625253/brc_scratch/Programs/bismark_v0.18.1/bismark_methylation_extractor -p --gzip --parallel 8 --report -o output/methylationExtractor/ $inr1



