#!/bin/bash

#SBATCH --partition=main    # Partition (job queue)
#SBATCH --job-name=deseq2        # Assign an 8-character name to your job, no spaces
#SBATCH --nodes=1                # Number of compute nodes
#SBATCH --ntasks=1               # Processes (usually = cores) on each node
#SBATCH --cpus-per-task=28       # Threads per process (or per core)
#SBATCH --export=ALL             # Export you current environment settings to the job environment
#SBATCH --time=2:00:00
#SBATCH --mem=100G
#SBATCH --output=rlognorm.%N.%j.out

#rlog-normalize htseq-cts
R --no-save < rhi.RNAseq.rlog.norm.R >& rhi.DRIPseq.RNAseq.kseek.rlog.norm.R.log
