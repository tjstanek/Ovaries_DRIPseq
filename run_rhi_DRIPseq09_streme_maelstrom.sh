#!/bin/bash

#SBATCH --partition=genetics_1        # Partition (job queue)
#SBATCH --requeue                     # Return job to the queue if preempted
#SBATCH --job-name=streme              # Assign an short name to your job
#SBATCH --nodes=1                     # Number of nodes you require
#SBATCH --ntasks=1                    # Total # of tasks across all nodes
#SBATCH --cpus-per-task=4            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                    # Real memory (RAM) required (MB)

#SBATCH --time=1:00:00               # Total run time limit (HH:MM:SS)
#SBATCH --output=Dcon.mloop.motif.%N.%j.out  # STDOUT output file
#SBATCH --error=Dcon.mloop.motif.%N.%j.err   # STDERR output file (optional)
#SBATCH --export=ALL                  # Export you current env to the job env


#Search for low-complexity motifs in ovarian TE R-loop peakset using MEME (anr mode)
meme Dcon.mloop.fasta -objfun de -neg Ctrl.mloop.split.fasta -hsfrac 0.67 -dna -mod anr -nmotifs 50 -oc Dcon.de_anr

#Test for enrichment of adult female-enriched genic R-loop motifs in ovarian genic R-loop motifs
~/meme/bin/sea sea --verbosity 1 --thresh 10.0 --p Diffbind.note.fasta --m ../FE_1e-2/streme.txt --oc rhinonDE.AFE_sea

#Test for enrichment of ovarian R-loop TE motifs in adult genic R-loop peakset
~/meme/bin/sea --verbosity 1 --thresh 10.0  --p Adults_streme/nonDE.all.fasta --m Dcon.de_anr/meme.txt --oc AnonDE.Dcon.mloop_de_anr

#Test for enrichment of adult genic common R-loop motifs in ovarian TE R-loop peakset
~/meme/bin/sea --verbosity 1 --thresh 10.0  --p Dcon.mloop.fasta --m ../AnonDE_1e-10/streme.txt --oc Dcon.mloop.AnonDE_sea


