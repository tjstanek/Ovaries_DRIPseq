#!/bin/bash

#SBATCH --partition=main       # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=drip-trim          # Assign a short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=28            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                   # Real memory (RAM) required (MB)
#SBATCH --array=1-2                  # Array range
#SBATCH --time=72:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=drip.rh.align.%A_%a.out   # STDOUT output file
#SBATCH --error=drip.rh.align.%A_%a.err    # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

##Align trimmed reads with bowtie
#...align to tidal.fasta to get TE counts immediately
echo "Aligning RNaseH reads to dm6-TIDAL, GRCh38 genomes..."
bowtie2 --no-mixed --no-discordant --dovetail --very-sensitive --phred33 -X 1000 -q -x tidal.fasta -1 ${SLURM_ARRAY_TASK_ID}.rh_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}.rh_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.rh.tidal.sam --un-conc ${SLURM_ARRAY_TASK_ID}.rh.dm.unaligned.fastq --al-conc ${SLURM_ARRAY_TASK_ID}.rh.dm.aligned.fastq

echo "Aligning DRIP-seq reads to dm6-TIDAL genome..."
bowtie2 --no-mixed --no-discordant --dovetail --very-sensitive --phred33 -X 1000 -q -x tidal.fasta -1 ${SLURM_ARRAY_TASK_ID}.r1_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}.r1_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.r1.tidal.sam
bowtie2 --no-mixed --no-discordant --dovetail --very-sensitive --phred33 -X 1000 -q -x tidal.fasta -1 ${SLURM_ARRAY_TASK_ID}.r2_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}.r2_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.r2.tidal.sam
bowtie2 --no-mixed --no-discordant --dovetail --very-sensitive --phred33 -X 1000 -q -x tidal.fasta -1 ${SLURM_ARRAY_TASK_ID}.Input_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}.Input_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.Input.tidal.sam

##Convert sam to bam
echo "Converting sam to bam..."
module load samtools
samtools view -S -b ${SLURM_ARRAY_TASK_ID}.rh.tidal.sam | samtools view -q 20 -b | samtools sort > ${SLURM_ARRAY_TASK_ID}.rh.tidal.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.rh.tidal.sorted.bam

samtools view -S -b ${SLURM_ARRAY_TASK_ID}.r1.tidal.sam | samtools view -q 20 -b | samtools sort > ${SLURM_ARRAY_TASK_ID}.r1.tidal.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.r1.tidal.sorted.bam

samtools view -S -b ${SLURM_ARRAY_TASK_ID}.r2.tidal.sam | samtools view -q 20 -b | samtools sort > ${SLURM_ARRAY_TASK_ID}.r2.tidal.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.r2.tidal.sorted.bam

samtools view -S -b ${SLURM_ARRAY_TASK_ID}.Input.tidal.sam | samtools view -q 20 -b | samtools sort > ${SLURM_ARRAY_TASK_ID}.Input.tidal.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.Input.tidal.sorted.bam

echo "Done!"
