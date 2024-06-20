#!/bin/bash

#SBATCH --partition=main             # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=RNA_1_r2          # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                  # Total # of tasks across all nodes
#SBATCH --cpus-per-task=8            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                 # Real memory (RAM) required (MB)
#SBATCH --array=1-2                # Array range

#SBATCH --time=72:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=RNAseq.%A_%a.out     # STDOUT output file
#SBATCH --error=RNAseq.%A_%a.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

# #Trim RNA reads
module load java
java -jar /home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -trimlog ${SLURM_ARRAY_TASK_ID}.r1_trim.log ${SLURM_ARRAY_TASK_ID}.r1_1.fq.gz ${SLURM_ARRAY_TASK_ID}.r1_2.fq.gz ${SLURM_ARRAY_TASK_ID}.r1_1.paired.fastq ${SLURM_ARRAY_TASK_ID}.r1_1.unpaired.fastq ${SLURM_ARRAY_TASK_ID}.r1_2.paired.fastq ${SLURM_ARRAY_TASK_ID}.r1_2.unpaired.fastq  ILLUMINACLIP:/home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:3:10 MINLEN:36

java -jar /home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -trimlog ${SLURM_ARRAY_TASK_ID}.r2_trim.log ${SLURM_ARRAY_TASK_ID}.r2_1.fq.gz ${SLURM_ARRAY_TASK_ID}.r2_2.fq.gz ${SLURM_ARRAY_TASK_ID}.r2_1.paired.fastq ${SLURM_ARRAY_TASK_ID}.r2_1.unpaired.fastq ${SLURM_ARRAY_TASK_ID}.r2_2.paired.fastq ${SLURM_ARRAY_TASK_ID}.r2_2.unpaired.fastq  ILLUMINACLIP:/home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:3:10 MINLEN:36

#Bowtie2 to align reads to TEs
module load bowtie2
bowtie2 --no-mixed --no-discordant --dovetail --very-sensitive --local --phred33 -X 1000 -q -x ../DRIP/combined2.fasta -1 ${SLURM_ARRAY_TASK_ID}.r1_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}.r1_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.r1.combined2.sam

bowtie2 --no-mixed --no-discordant --dovetail --very-sensitive --local --phred33 -X 1000 -q -x ../DRIP/combined2.fasta -1 ${SLURM_ARRAY_TASK_ID}.r2_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}.r2_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.r2.combined2.sam

module load samtools
samtools view -S -b ${SLURM_ARRAY_TASK_ID}.r1.combined2.sam | samtools view -q 20 -b | samtools sort > ${SLURM_ARRAY_TASK_ID}.r1.combined2.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.r1.combined2.sorted.bam

samtools view -S -b ${SLURM_ARRAY_TASK_ID}.r2.combined2.sam | samtools view -q 20 -b | samtools sort > ${SLURM_ARRAY_TASK_ID}.r2.combined2.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.r2.combined2.sorted.bam

#For genic transcripts, run HISAT2...
module load HISAT2/2.1.0

#...first to remove rRNA reads...
hisat2 -p 4 -x /home/ts928/Reference/bdgp6_rrna/genome_rrna -1 ${SLURM_ARRAY_TASK_ID}.r1_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}.r1_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.r1.rrna.sam --un-conc ${SLURM_ARRAY_TASK_ID}.r1-rrna.fastq

hisat2 -p 4 -x /home/ts928/Reference/bdgp6_rrna/genome_rrna -1 ${SLURM_ARRAY_TASK_ID}.r2_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}.r2_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.r2.rrna.sam --un-conc ${SLURM_ARRAY_TASK_ID}.r2-rrna.fastq

#...then to align to transcriptome.
hisat2 -p 4 -x /home/ts928/Reference/bdgp6_tran/genome_tran --novel-splicesite-outfile ${SLURM_ARRAY_TASK_ID}.r1_novel.splice.out.txt --rna-strandness RF -1 ${SLURM_ARRAY_TASK_ID}.r1-rrna.1.fastq -2 ${SLURM_ARRAY_TASK_ID}.r1-rrna.2.fastq -S ${SLURM_ARRAY_TASK_ID}.r1_aligned.sam

hisat2 -p 4 -x /home/ts928/Reference/bdgp6_tran/genome_tran --novel-splicesite-outfile ${SLURM_ARRAY_TASK_ID}.r2_novel.splice.out.txt --rna-strandness RF -1 ${SLURM_ARRAY_TASK_ID}.r2-rrna.1.fastq -2 ${SLURM_ARRAY_TASK_ID}.r2-rrna.2.fastq -S ${SLURM_ARRAY_TASK_ID}.r2_aligned.sam

#PICARD
java -jar /home/ts928/Genomics_Workshop/Programs/picard-tools-1.119/SortSam.jar SO=coordinate INPUT=${SLURM_ARRAY_TASK_ID}.r1_aligned.sam OUTPUT=${SLURM_ARRAY_TASK_ID}.r1_aligned.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

java -jar /home/ts928/Genomics_Workshop/Programs/picard-tools-1.119/SortSam.jar SO=coordinate INPUT=${SLURM_ARRAY_TASK_ID}.r2_aligned.sam OUTPUT=${SLURM_ARRAY_TASK_ID}.r2_aligned.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

# #htseq-ct
module load samtools intel/17.0.2 python/2.7.12

samtools sort -n  ${SLURM_ARRAY_TASK_ID}_r1_aligned.bam | samtools view | htseq-count -m intersect-nonempty -t exon -i gene_id -s reverse --additional-attr=gene_symbol  - /projects/genetics/ellison_lab/genomes/bed/dmel-all-r6.27.gtf > ${SLURM_ARRAY_TASK_ID}_r1.rev.txt

samtools sort -n  ${SLURM_ARRAY_TASK_ID}_r1_aligned.bam | samtools view | htseq-count -m intersect-nonempty -t exon -i gene_id -s reverse --additional-attr=gene_symbol  - /projects/genetics/ellison_lab/genomes/bed/dmel-all-r6.27.gtf > ${SLURM_ARRAY_TASK_ID}_r2.rev.txt

#Pass to R for rlog-normalization


#For visualization in IGV, filter only reads on TEs
samtools view -b -L transposon.sorted.bed -o ${SLURM_ARRAY_TASK_ID}.r1.te.RNA.bam ${SLURM_ARRAY_TASK_ID}.r1_aligned.bam
samtools index ${SLURM_ARRAY_TASK_ID}.r1.te.RNA.bam

samtools view -b -L transposon.sorted.bed -o ${SLURM_ARRAY_TASK_ID}.r2.te.RNA.bam ${SLURM_ARRAY_TASK_ID}.r2_aligned.bam
samtools index ${SLURM_ARRAY_TASK_ID}.r2.te.RNA.bam
