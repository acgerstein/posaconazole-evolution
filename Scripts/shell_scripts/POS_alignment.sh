#!/bin/bash
#SBATCH --account=def-acgerste
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=0
#SBATCH --mail-user=adamubua@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=alignment_POS
#SBATCH --output=%x-%j.out


module load nixpkgs/16.09 gcc/7.3.0 intel/2018.3
module load bwa/0.7.17
module load picard/2.23.2 gatk/4.1.2.0

output="/home/abdul/scratch/POS_Analyses/data_out/alignment"
input="/home/abdul/scratch/POS_Analyses/data_out/trimmed_reads"
ref="/home/abdul/scratch/POS_Analyses/reference"

for i in $input/*_R1.trimmed_PE.fastq.gz;
do
withpath="${i}" 
filename=${withpath##*/}
base="${filename%*_*.trimmed_PE.fastq.gz}"
sample_name=`echo "${base}" | awk -F ".trimmed_PE.fastq.gz" '{print $1}'`


mkdir $output/"${base}"

time bwa mem -t $SLURM_CPUS_PER_TASK $ref/C_albicans_A21.fasta $input/"${base}"_R1.trimmed_PE.fastq.gz $input/"${base}"_R2.trimmed_PE.fastq.gz -o $output/"${base}"/"${base}".sam 
echo ""${base}" alignment done!"

time java -jar $EBROOTPICARD/picard.jar SortSam I=$output/"${base}"/"${base}".sam O=$output/"${base}"/"${base}".sorted.sam SORT_ORDER=coordinate 
echo ""${base}"  sort sam done!"

time java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics R=$ref/C_albicans_A21.fasta I=$output/"${base}"/"${base}".sorted.sam O=$output/"${base}"/"${base}".alignment_summary.txt
echo ""${base}"  alignment_summary sort done!"

time java -jar $EBROOTPICARD/picard.jar SamFormatConverter I=$output/"${base}"/"${base}".sorted.sam O=$output/"${base}"/"${base}".sorted.bam R=$ref/C_albicans_A21.fasta
echo ""${base}"  SamFormatConverter done!"

time java -jar $EBROOTPICARD/picard.jar ValidateSamFile I=$output/"${base}"/"${base}".sorted.bam R=$ref/C_albicans_A21.fasta

echo ""${base}"  ValidateSamFile done!"

time java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=$output/"${base}"/"${base}".sorted.bam R=$ref/C_albicans_A21.fasta O=$output/"${base}"/"${base}".sorted.RG.bam RGID=Cell1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM="${base}"

echo ""${base}"   Read Group done!"


#Validate the BAM file.

time java -jar $EBROOTPICARD/picard.jar ValidateSamFile I=$output/"${base}"/"${base}".sorted.RG.bam R=$ref/C_albicans_A21.fasta IGNORE=INVALID_TAG_NM

#Mark potential reads duplicates

time java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$output/"${base}"/"${base}".sorted.RG.bam O=$output/"${base}"/"${base}".sorted_RG_dedup.bam R=$ref/C_albicans_A21.fasta M=$output/"${base}"/"${base}".sorted_RG_dedup.txt CREATE_INDEX=true

echo ""${base}"  Mark duplicates done!"


#Correct possible info differences in the alignmed paired end reads. Readmore on this step

time java -jar $EBROOTPICARD/picard.jar FixMateInformation I=$output/"${base}"/"${base}".sorted_RG_dedup.bam R=$ref/C_albicans_A21.fasta O=$output/"${base}"/"${base}".sorted_RG_dedup_Fixmate.bam ADD_MATE_CIGAR=true CREATE_INDEX=true

echo ""${base}"  Fixmate information done!"



# generate coverage files
bedtools genomecov -d  -ibam $output/"${base}"/"${base}".sorted_RG_dedup_Fixmate.bam > "${base}".txt

#Insert headers
echo -e "Chr\\tlocus\\t"${base}"" | cat - "${base}".txt > "${base}".coverage.txt

done
