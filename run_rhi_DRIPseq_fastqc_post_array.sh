#!/bin/bash

#SBATCH --partition=genetics_1             # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=fastqc_post_array          # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                  # Total # of tasks across all nodes
#SBATCH --cpus-per-task=28            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=80G                 # Real memory (RAM) required (MB)
#SBATCH --array=1-2                # Array range

#SBATCH --time=16:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=fastqc.%A_%a.out     # STDOUT output file
#SBATCH --error=fastqc.%A_%a.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

module load java
fastqc -o DRIP_fastqc_posttrim ${SLURM_ARRAY_TASK_ID}.r1_1.paired.fastq ${SLURM_ARRAY_TASK_ID}.r1_2.paired.fastq
fastqc -o DRIP_fastqc_posttrim ${SLURM_ARRAY_TASK_ID}.r2_1.paired.fastq ${SLURM_ARRAY_TASK_ID}.r2_2.paired.fastq
fastqc -o DRIP_fastqc_posttrim ${SLURM_ARRAY_TASK_ID}.Input_1.paired.fastq ${SLURM_ARRAY_TASK_ID}.Input_2.paired.fastq
fastqc -o DRIP_fastqc_posttrim ${SLURM_ARRAY_TASK_ID}.rh_1.paired.fastq ${SLURM_ARRAY_TASK_ID}.rh_2.paired.fastq
