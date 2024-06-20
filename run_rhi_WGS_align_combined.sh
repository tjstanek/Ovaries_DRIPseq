#!/bin/bash

#SBATCH --partition=genetics_1       # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=combined          # Assign a short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=28            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=192G                   # Real memory (RAM) required (MB)
#SBATCH --array=1-2                  # Array range
#SBATCH --time=96:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=wgs.combined.align.%A_%a.out   # STDOUT output file
#SBATCH --error=wgs.combined.align.%A_%a.err    # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env


##Align trimmed reads with bowtie
#echo "Aligning reads..."
module load bowtie2

#Align to combined.fasta to get TE counts immediately
echo "Aligning WGS reads to dm6-combined genome..."
bowtie2 --no-mixed --no-discordant --dovetail --very-sensitive --local --phred33 -X 1000 -q -x combined2.fasta -1 ${SLURM_ARRAY_TASK_ID}_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.combined2.sam
echo "Alignment done!"

##Convert sam to bam
echo "Converting sam to bam..."
module load samtools
samtools view -S -b ${SLURM_ARRAY_TASK_ID}.combined2.sam | samtools view -q 20 -b | samtools sort > ${SLURM_ARRAY_TASK_ID}.combined2.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.combined2.sorted.bam
echo "Conversion to bam done!"
