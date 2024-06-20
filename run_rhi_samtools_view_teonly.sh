#!/bin/bash

#SBATCH --partition=main# Partition (job queue)
#SBATCH --requeue                  # Return job to queue if preempted
#SBATCH --job-name=samtools        # Assign a short name to your job
#SBATCH --nodes=1                  # Numbers of nodes you require
#SBATCH --ntasks=1                 # Total # of tasks across all nodes
#SBATCH --cpus-per-task=12          # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                 # Real memory (RAM) required (MB)
#SBATCH --array=2               # Array range

#SBATCH --time=48:00:00            # Total run time limit (HH:MM:SS)
#SBATCH --output=rhi_samtools_tidalonly.%A_%a.out       #STDOUT output file
#SBATCH --error=rhi_samtools_tidalonly.%A_%a.err        #STDERR output file (optional)
#SBATCH --export=ALL               # Export you current env to the job env
#SBATCH --mail-type=END,ARRAY_TASKS   # Notify user by mail when certain event types occur
#SBATCH --mail-user=ts928@hginj.rutgers.edu  #User email

module load samtools

#filter only reads on TEs
echo "Filtering primary chroms..."
samtools view -b -L transposon.sorted.bed -o ${SLURM_ARRAY_TASK_ID}.r1.te.only.bam ${SLURM_ARRAY_TASK_ID}.r1.combined2.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.r1.te.only.bam

samtools view -b -L transposon.sorted.bed -o ${SLURM_ARRAY_TASK_ID}.r2.te.only.bam ${SLURM_ARRAY_TASK_ID}.r2.combined2.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.r2.te.only.bam

samtools view -b -L transposon.sorted.bed -o ${SLURM_ARRAY_TASK_ID}.rh.te.only.bam ${SLURM_ARRAY_TASK_ID}.rh.combined2.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.rh.te.only.bam

samtools view -b -L transposon.sorted.bed -o ${SLURM_ARRAY_TASK_ID}.Input.te.only.bam ${SLURM_ARRAY_TASK_ID}.Input.combined2.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.Input.te.only.bam


echo "Done!"
