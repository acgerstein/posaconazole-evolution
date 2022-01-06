#!/bin/bash
#SBATCH --account=def-acgerste
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=05:00:00
#SBATCH --mem=0
#SBATCH --mail-user=adamubua@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=Trim
#SBATCH --output=%x-%j.out

module load  StdEnv/2020 trimmomatic/0.39

input="/home/abdul/scratch/POS_Analyses/data_in"
output="/home/abdul/scratch/POS_Analyses/data_out/trimmed_reads"

for i in $input/*_R1.fastq.gz;
do
withpath="${i}" filename=${withpath##*/}
base="${filename%*_*R1.fastq.gz}"
sample_name=`echo "${base}" | awk -F "_R1.fastq.gz" '{print $1}'`

time java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $SLURM_CPUS_PER_TASK -trimlog $output/"${base}".log $input/"${base}"_R1.fastq.gz $input/"${base}"_R2.fastq.gz $output/"${base}"_R1.trimmed_PE.fastq.gz $output/"${base}"_R1.trimmed_SE.fastq.gz $output/"${base}"_R2.trimmed_PE.fastq.gz $output/"${base}"_R2.trimmed_SE.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33

done

