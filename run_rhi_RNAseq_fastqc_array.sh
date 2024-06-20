#!/bin/bash
#SBATCH --partition=genetics_1             # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=fastqc          # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                  # Total # of tasks across all nodes
#SBATCH --cpus-per-task=8            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=40G                 # Real memory (RAM) required (MB)
#SBATCH --array=1-2                # Array range

#SBATCH --time=2:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=fastqc.%A_%a.out     # STDOUT output file
#SBATCH --error=fastqc.%A_%a.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

module load java
fastqc -o RNA_fastqc_pretrim ${SLURM_ARRAY_TASK_ID}.r1_1.fq.gz ${SLURM_ARRAY_TASK_ID}.r1_2.fq.gz
fastqc -o RNA_fastqc_pretrim ${SLURM_ARRAY_TASK_ID}.r2_1.fq.gz ${SLURM_ARRAY_TASK_ID}.r2_2.fq.gz
