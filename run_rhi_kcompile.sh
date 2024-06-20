#!/bin/bash
#SBATCH --partition=genetics_1     # Partition (job queue)
#SBATCH --requeue                  # Return job to queue if preempted
#SBATCH --job-name=kcompile        # Assign a short name to your job
#SBATCH --nodes=1                  # Numbers of nodes you require
#SBATCH --ntasks=1                 # Total # of tasks across all nodes
#SBATCH --cpus-per-task=8          # Cores per task (>1 if multithread tasks)
#SBATCH --mem=80G                 # Real memory (RAM) required (MB)

#SBATCH --time=1:00:00            # Total run time limit (HH:MM:SS)
#SBATCH --output=kcompile.%N.%j.out       #STDOUT output file
#SBATCH --error=kcompile.%N.%j.err        #STDERR output file (optional)
#SBATCH --export=ALL               # Export you current env to the job env

#Compile k-seek analyses for WGS, DRIPseq, & RNAseq
perl /projects/genetics/ellison_lab/rohan/k-seek/k_compiler.pl kseek/ rhi.wgs.dripseq.rnaseq.total.tab

perl transpose.pl rhi.wgs.dripseq.rnaseq.total.tab.rep.compiled rhi.wgs.dripseq.rnaseq.transposed
