#!/bin/bash

#SBATCH --partition=main     # Partition (job queue)
#SBATCH --requeue                  # Return job to queue if preempted
#SBATCH --job-name=kseek           # Assign a short name to your job
#SBATCH --nodes=1                  # Numbers of nodes you require
#SBATCH --ntasks=1                 # Total # of tasks across all nodes
#SBATCH --cpus-per-task=4          # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                 # Real memory (RAM) required (MB)
#SBATCH --array=1-3                  # Array range
#SBATCH --time=72:00:00            # Total run time limit (HH:MM:SS)
#SBATCH --output=kseek.%A_%a.out   #STDOUT output file
#SBATCH --error=kseek.%A_%a.err    #STDERR output file (optional)
#SBATCH --export=ALL               # Export you current env to the job env

##WGS
perl /projects/genetics/ellison_lab/rohan/k-seek/k_seek.pl ../WGS/${SLURM_ARRAY_TASK_ID}_1.paired.fastq WGS.${SLURM_ARRAY_TASK_ID}_1

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ../WGS/${SLURM_ARRAY_TASK_ID}_2.paired.fastq WGS.${SLURM_ARRAY_TASK_ID}_2

##DRIP
perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ${SLURM_ARRAY_TASK_ID}.r1_1.paired.fastq DRIP.${SLURM_ARRAY_TASK_ID}.r1_1

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ${SLURM_ARRAY_TASK_ID}.r1_2.paired.fastq DRIP.${SLURM_ARRAY_TASK_ID}.r1_2

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ${SLURM_ARRAY_TASK_ID}.r2_1.paired.fastq DRIP.${SLURM_ARRAY_TASK_ID}.r2_1

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ${SLURM_ARRAY_TASK_ID}.r2_2.paired.fastq DRIP.${SLURM_ARRAY_TASK_ID}.r2_2

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ${SLURM_ARRAY_TASK_ID}.Input_1.paired.fastq DRIP.${SLURM_ARRAY_TASK_ID}.Input_1

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ${SLURM_ARRAY_TASK_ID}.Input_2.paired.fastq DRIP.${SLURM_ARRAY_TASK_ID}.Input_2

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ${SLURM_ARRAY_TASK_ID}.rh.dm.aligned.1.fastq DRIP.${SLURM_ARRAY_TASK_ID}.rh_1

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ${SLURM_ARRAY_TASK_ID}.rh.dm.aligned.2.fastq DRIP.${SLURM_ARRAY_TASK_ID}.rh_2

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ${SLURM_ARRAY_TASK_ID}.rh.dm.hs.aligned.1.fastq DRIP.${SLURM_ARRAY_TASK_ID}.rh_1.dmhs

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ${SLURM_ARRAY_TASK_ID}.rh.dm.hs.aligned.2.fastq DRIP.${SLURM_ARRAY_TASK_ID}.rh_2.dmhs

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ${SLURM_ARRAY_TASK_ID}.rh.noaligned.1.fastq DRIP.${SLURM_ARRAY_TASK_ID}.rh_1.none

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ${SLURM_ARRAY_TASK_ID}.rh.noaligned.2.fastq DRIP.${SLURM_ARRAY_TASK_ID}.rh_2.none

##RNA
perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ../RNA/${SLURM_ARRAY_TASK_ID}.r1_1.paired.fastq RNA.${SLURM_ARRAY_TASK_ID}.r1_1

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ../RNA/${SLURM_ARRAY_TASK_ID}.r1_2.paired.fastq RNA.${SLURM_ARRAY_TASK_ID}.r1_2

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ../RNA/${SLURM_ARRAY_TASK_ID}.r2_1.paired.fastq RNA.${SLURM_ARRAY_TASK_ID}.r2_1

perl /projects/genetics/ellison_lab/tim/k-seek/k_seek.pl ../RNA/${SLURM_ARRAY_TASK_ID}.r2_2.paired.fastq RNA.${SLURM_ARRAY_TASK_ID}.r2_2
