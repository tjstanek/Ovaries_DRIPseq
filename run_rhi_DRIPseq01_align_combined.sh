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
#SBATCH --output=drip.combined2.align.%A_%a.out   # STDOUT output file
#SBATCH --error=drip.combined2.align.%A_%a.err    # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

#Trim reads with trimmomatic
echo "Trimming reads..."
module load java
java -jar /home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -trimlog ${SLURM_ARRAY_TASK_ID}.r1_trim.log ${SLURM_ARRAY_TASK_ID}.r1_1.fq.gz ${SLURM_ARRAY_TASK_ID}.r1_2.fq.gz ${SLURM_ARRAY_TASK_ID}.r1_1.paired.fastq ${SLURM_ARRAY_TASK_ID}.r1_1.unpaired.fastq ${SLURM_ARRAY_TASK_ID}.r1_2.paired.fastq ${SLURM_ARRAY_TASK_ID}.r1_2.unpaired.fastq  ILLUMINACLIP:/home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:15 TRAILING:15 SLIDINGWINDOW:3:10 MINLEN:36

java -jar /home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -trimlog ${SLURM_ARRAY_TASK_ID}.r2_trim.log ${SLURM_ARRAY_TASK_ID}.r2_1.fq.gz ${SLURM_ARRAY_TASK_ID}.r2_2.fq.gz ${SLURM_ARRAY_TASK_ID}.r2_1.paired.fastq ${SLURM_ARRAY_TASK_ID}.r2_1.unpaired.fastq ${SLURM_ARRAY_TASK_ID}.r2_2.paired.fastq ${SLURM_ARRAY_TASK_ID}.r2_2.unpaired.fastq  ILLUMINACLIP:/home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:15 TRAILING:15 SLIDINGWINDOW:3:10 MINLEN:36

java -jar /home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -trimlog ${SLURM_ARRAY_TASK_ID}.Input_trim.log ${SLURM_ARRAY_TASK_ID}.Input_1.fq.gz ${SLURM_ARRAY_TASK_ID}.Input_2.fq.gz ${SLURM_ARRAY_TASK_ID}.Input_1.paired.fastq ${SLURM_ARRAY_TASK_ID}.Input_1.unpaired.fastq ${SLURM_ARRAY_TASK_ID}.Input_2.paired.fastq ${SLURM_ARRAY_TASK_ID}.Input_2.unpaired.fastq  ILLUMINACLIP:/home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:15 TRAILING:15 SLIDINGWINDOW:3:10 MINLEN:36

java -jar /home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -trimlog ${SLURM_ARRAY_TASK_ID}.rh_trim.log ${SLURM_ARRAY_TASK_ID}.rh_1.fq.gz ${SLURM_ARRAY_TASK_ID}.rh_2.fq.gz ${SLURM_ARRAY_TASK_ID}.rh_1.paired.fastq ${SLURM_ARRAY_TASK_ID}.rh_1.unpaired.fastq ${SLURM_ARRAY_TASK_ID}.rh_2.paired.fastq ${SLURM_ARRAY_TASK_ID}.rh_2.unpaired.fastq  ILLUMINACLIP:/home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:15 TRAILING:15 SLIDINGWINDOW:3:10 MINLEN:36

#Align trimmed reads with bowtie to combined.fasta to get TE counts immediately
echo "Aligning RNase-H-treated reads to dm6-combined genome..."
bowtie2 --no-mixed --no-discordant --dovetail --very-sensitive --local --phred33 -X 1000 -q -x combined2.fasta -1 ${SLURM_ARRAY_TASK_ID}.rh.dm.aligned.1.fastq -2 ${SLURM_ARRAY_TASK_ID}.rh.dm.aligned.2.fastq -S ${SLURM_ARRAY_TASK_ID}.rh.combined2.sam
bowtie2 --no-mixed --no-discordant --dovetail --phred33 -X 1000 -q -x GRCh38_noalt_as/GRCh38_noalt_as -1 ${SLURM_ARRAY_TASK_ID}.rh.dm.unaligned.1.fastq -2 ${SLURM_ARRAY_TASK_ID}.rh.dm.unaligned.2.fastq -S ${SLURM_ARRAY_TASK_ID}.rh.hs.aligned.sam --un-conc ${SLURM_ARRAY_TASK_ID}.rh.noaligned.fastq --al-conc ${SLURM_ARRAY_TASK_ID}.rh.hs.aligned.fastq

cat ${SLURM_ARRAY_TASK_ID}.rh.dm.aligned.1.fastq ${SLURM_ARRAY_TASK_ID}.rh.hs.aligned.1.fastq > ${SLURM_ARRAY_TASK_ID}.rh.dm.hs.aligned.1.fastq
cat ${SLURM_ARRAY_TASK_ID}.rh.dm.aligned.2.fastq ${SLURM_ARRAY_TASK_ID}.rh.hs.aligned.2.fastq > ${SLURM_ARRAY_TASK_ID}.rh.dm.hs.aligned.2.fastq

echo "Aligning DRIP-seq reads to dm6-combined genome..."
bowtie2 --no-mixed --no-discordant --dovetail --very-sensitive --local --phred33 -X 1000 -q -x combined2.fasta -1 ${SLURM_ARRAY_TASK_ID}.r1_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}.r1_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.r1.combined2.sam
bowtie2 --no-mixed --no-discordant --dovetail --very-sensitive --local --phred33 -X 1000 -q -x combined2.fasta -1 ${SLURM_ARRAY_TASK_ID}.r2_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}.r2_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.r2.combined2.sam
bowtie2 --no-mixed --no-discordant --dovetail --very-sensitive --local --phred33 -X 1000 -q -x combined2.fasta -1 ${SLURM_ARRAY_TASK_ID}.Input_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}.Input_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.Input.combined2.sam

##Convert sam to bam
echo "Converting sam to bam..."
module load samtools
samtools view -S -b ${SLURM_ARRAY_TASK_ID}.rh.combined2.sam | samtools view -q 20 -b | samtools sort > ${SLURM_ARRAY_TASK_ID}.rh.combined2.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.rh.combined2.sorted.bam

samtools view -S -b ${SLURM_ARRAY_TASK_ID}.r1.combined2.sam | samtools view -q 20 -b | samtools sort > ${SLURM_ARRAY_TASK_ID}.r1.combined2.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.r1.combined2.sorted.bam

samtools view -S -b ${SLURM_ARRAY_TASK_ID}.r2.combined2.sam | samtools view -q 20 -b | samtools sort > ${SLURM_ARRAY_TASK_ID}.r2.combined2.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.r2.combined2.sorted.bam

samtools view -S -b ${SLURM_ARRAY_TASK_ID}.Input.combined2.sam | samtools view -q 20 -b | samtools sort > ${SLURM_ARRAY_TASK_ID}.Input.combined2.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.Input.combined2.sorted.bam


##Peakcalling
echo "Calling peaks with MACS2..."
macs2 callpeak --outdir idr -t ${SLURM_ARRAY_TASK_ID}.r1.mq20.sorted.bam -c ${SLURM_ARRAY_TASK_ID}.Input.mq20.sorted.bam -f BAMPE -g dm -n ${SLURM_ARRAY_TASK_ID}.r1 -B -p 1e-3 2> idr/macs/${SLURM_ARRAY_TASK_ID}.r1_macs2.log
sort -k8,8nr idr/${SLURM_ARRAY_TASK_ID}.r1_peaks.narrowPeak > idr/macs/${SLURM_ARRAY_TASK_ID}.r1_peaks.narrowPeak

macs2 callpeak --outdir idr -t ${SLURM_ARRAY_TASK_ID}.r2.mq20.sorted.bam -c ${SLURM_ARRAY_TASK_ID}.Input.mq20.sorted.bam -f BAMPE -g dm -n ${SLURM_ARRAY_TASK_ID}.r2 -B -p 1e-3 2> idr/macs/${SLURM_ARRAY_TASK_ID}.r2_macs2.log
sort -k8,8nr idr/${SLURM_ARRAY_TASK_ID}.r2_peaks.narrowPeak > idr/macs/${SLURM_ARRAY_TASK_ID}.r2_peaks.narrowPeak

#...call peaks on te-containing bams
macs2 callpeak --outdir idr -t ${SLURM_ARRAY_TASK_ID}.rh.te.mq20.sorted.bam -c ${SLURM_ARRAY_TASK_ID}.Input.te.mq20.sorted.bam -f BAMPE -g dm -n ${SLURM_ARRAY_TASK_ID}.rh -B -p 1e-3 2> idr/macs/${SLURM_ARRAY_TASK_ID}.rh_macs2.log
sort -k8,8nr idr/${SLURM_ARRAY_TASK_ID}.rh_peaks.narrowPeak > idr/macs/${SLURM_ARRAY_TASK_ID}.rh_peaks.narrowPeak

#...call peaks using RNaseH+ bams as -c
macs2 callpeak --outdir idr -t ${SLURM_ARRAY_TASK_ID}.r1.te.mq20.sorted.bam -c ${SLURM_ARRAY_TASK_ID}.rh.te.mq20.sorted.bam -f BAMPE -g dm -n ${SLURM_ARRAY_TASK_ID}.r1c -B -p 1e-3 2> idr/macs/${SLURM_ARRAY_TASK_ID}.r1c_macs2.log
sort -k8,8nr idr/${SLURM_ARRAY_TASK_ID}.r1c_peaks.narrowPeak > idr/macs/${SLURM_ARRAY_TASK_ID}.r1c_peaks.narrowPeak

macs2 callpeak --outdir idr -t ${SLURM_ARRAY_TASK_ID}.r2.te.mq20.sorted.bam -c ${SLURM_ARRAY_TASK_ID}.rh.te.mq20.sorted.bam -f BAMPE -g dm -n ${SLURM_ARRAY_TASK_ID}.r2c -B -p 1e-3 2> idr/macs/${SLURM_ARRAY_TASK_ID}.r2c_macs2.log
sort -k8,8nr idr/${SLURM_ARRAY_TASK_ID}.r2c_peaks.narrowPeak > idr/macs/${SLURM_ARRAY_TASK_ID}.r2c_peaks.narrowPeak

macs2 callpeak --outdir idr -t ${SLURM_ARRAY_TASK_ID}.rh2.dm.aligned.sorted.bam -c ${SLURM_ARRAY_TASK_ID}.Input.te.mq20.sorted.bam -f BAMPE -g dm -n ${SLURM_ARRAY_TASK_ID}.rh2.dm.aligned -B -p 1e-3 2> idr/macs/${SLURM_ARRAY_TASK_ID}.rh2.dm.aligned_macs2.log

macs2 callpeak --outdir idr -t ${SLURM_ARRAY_TASK_ID}.rh2.hs.unaligned.sorted.bam -c ${SLURM_ARRAY_TASK_ID}.Input.te.mq20.sorted.bam -f BAMPE -g dm -n ${SLURM_ARRAY_TASK_ID}.rh2.hs.unaligned -B -p 1e-3 2> idr/macs/${SLURM_ARRAY_TASK_ID}.rh2.hs.unaligned_macs2.log


##IDR analysis
idr --samples idr/${SLURM_ARRAY_TASK_ID}.r1_peaks.narrowPeak idr/${SLURM_ARRAY_TASK_ID}.r2_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${SLURM_ARRAY_TASK_ID}.idr -idr --plot --log-output-file ${SLURM_ARRAY_TASK_ID}.idr.log

cat ${SLURM_ARRAY_TASK_ID}.idr | awk '{if($5 >= 540) print $0}' > ${SLURM_ARRAY_TASK_ID}.idr.peaks.sig
wc -l ${SLURM_ARRAY_TASK_ID}.idr >> idr.commonpeaks.all
cat ${SLURM_ARRAY_TASK_ID}.idr | awk '{if($5 >= 540) print $0}' | wc -l >> idr.commonpeaks.sig

#Convert narrowPeak files into bed for diffbind
echo "Creating BED files..."
cut -f 1-6 ${SLURM_ARRAY_TASK_ID}.idr > macs2/${SLURM_ARRAY_TASK_ID}-idr.bed
echo "Done creating BED files!"

date

##IDR pseudorep on downsampled r1 and r2
#Set paths
macsDir=/scratch/ts928/novogene_DRIPseq/rhi/DRIP/chipseq-project/macs
outputDir=/scratch/ts928/novogene_DRIPseq/rhi/DRIP/chipseq-project/pooled_pseudoreps
tmpDir=/scratch/ts928/novogene_DRIPseq/rhi/DRIP/chipseq-project/tmp

echo "Merging BAM files for pseudoreplicates..."
samtools merge -u ${tmpDir}/${SLURM_ARRAY_TASK_ID}.merged.bam ${SLURM_ARRAY_TASK_ID}.r1.mq20.sorted.bam ${SLURM_ARRAY_TASK_ID}.r2.mq20.sorted.bam
samtools view -H ${tmpDir}/${SLURM_ARRAY_TASK_ID}.merged.bam > ${tmpDir}/${SLURM_ARRAY_TASK_ID}_header.sam

#Split merged treatments
nlines=$(samtools view ${tmpDir}/${SLURM_ARRAY_TASK_ID}.merged.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
samtools view ${tmpDir}/${SLURM_ARRAY_TASK_ID}.merged.bam | shuf - | split -d -l ${nlines} - "${tmpDir}/${SLURM_ARRAY_TASK_ID}" # This will shuffle the lines in the file and split it into two SAM files
cat ${tmpDir}/${SLURM_ARRAY_TASK_ID}_header.sam ${tmpDir}/${SLURM_ARRAY_TASK_ID}00 | samtools view -bS - > ${outputDir}/${SLURM_ARRAY_TASK_ID}00.bam
cat ${tmpDir}/${SLURM_ARRAY_TASK_ID}_header.sam ${tmpDir}/${SLURM_ARRAY_TASK_ID}01 | samtools view -bS - > ${outputDir}/${SLURM_ARRAY_TASK_ID}01.bam

#Peak calling on pseudoreplicates
echo "Calling peaks for pseudoreplicate1"
macs2 callpeak -t ${outputDir}/${SLURM_ARRAY_TASK_ID}00.bam -c ${SLURM_ARRAY_TASK_ID}.Input.mq20.sorted.bam -f BAMPE -g dm -n $macsDir/${SLURM_ARRAY_TASK_ID}.r1_pr -B -p 1e-3 2> $macsDir/${SLURM_ARRAY_TASK_ID}.r1_pr_macs2.log

echo "Calling peaks for pseudoreplicate2"
macs2 callpeak -t ${outputDir}/${SLURM_ARRAY_TASK_ID}01.bam -c ${SLURM_ARRAY_TASK_ID}.Input.mq20.sorted.bam -f BAMPE -g dm -n $macsDir/.${SLURM_ARRAY_TASK_ID}.r2_pr -B -p 1e-3 2> $macsDir/${SLURM_ARRAY_TASK_ID}.r2_pr_macs2.log

#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr $macsDir/${SLURM_ARRAY_TASK_ID}.r1_pr_peaks.narrowPeak | head -n 100000 > $macsDir/${SLURM_ARRAY_TASK_ID}.r1_pr_sorted.narrowPeak
sort -k8,8nr $macsDir/${SLURM_ARRAY_TASK_ID}.r2_pr_peaks.narrowPeak | head -n 100000 > $macsDir/${SLURM_ARRAY_TASK_ID}.r2_pr_sorted.narrowPeak

#Independent replicate IDR
echo "Running IDR on pseudoreplicates..."
idr --samples $macsDir/${SLURM_ARRAY_TASK_ID}.r1_pr_sorted.narrowPeak $macsDir/${SLURM_ARRAY_TASK_ID}.r2_pr_sorted.narrowPeak --input-file-type narrowPeak --output-file ${SLURM_ARRAY_TASK_ID}.pseudorep-idr_p1e-3 --rank p.value --plot --log-output-file ${SLURM_ARRAY_TASK_ID}.pr.idr.p1e-3.log

##IDR self-pseudorep on r1 and r2
#Same paths as IDR pseudorep above
#Merge treatment BAMS
echo "Creating header files for self-pseudoreplicates..."
samtools view -H ${SLURM_ARRAY_TASK_ID}.r1.mq20.sorted.bam > ${tmpDir}/${SLURM_ARRAY_TASK_ID}_r1_header.sam
samtools view -H ${SLURM_ARRAY_TASK_ID}.r2.mq20.sorted.bam > ${tmpDir}/${SLURM_ARRAY_TASK_ID}_r2_header.sam

#Split individual replicated treatments r1
nlines=$(samtools view ${SLURM_ARRAY_TASK_ID}.r1.mq20.sorted.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
samtools view ${SLURM_ARRAY_TASK_ID}.r1.mq20.sorted.bam | shuf - | split -d -l ${nlines} - "${tmpDir}/${SLURM_ARRAY_TASK_ID}_r1" # This will shuffle the lines in the file and split it into two SAM files
cat ${tmpDir}/${SLURM_ARRAY_TASK_ID}_r1_header.sam ${tmpDir}/${SLURM_ARRAY_TASK_ID}_r100 | samtools view -bS - > ${outputDir}/${SLURM_ARRAY_TASK_ID}_r100.bam
cat ${tmpDir}/${SLURM_ARRAY_TASK_ID}_r1_header.sam ${tmpDir}/${SLURM_ARRAY_TASK_ID}_r101 | samtools view -bS - > ${outputDir}/${SLURM_ARRAY_TASK_ID}_r101.bam

#Split individual replicated treatments new r2
nlines=$(samtools view ${SLURM_ARRAY_TASK_ID}.r2.mq20.sorted.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
samtools view ${SLURM_ARRAY_TASK_ID}.r2.mq20.sorted.bam | shuf - | split -d -l ${nlines} - "${tmpDir}/${SLURM_ARRAY_TASK_ID}_r2" # This will shuffle the lines in the file and split it into two SAM files
cat ${tmpDir}/${SLURM_ARRAY_TASK_ID}_r2_header.sam ${tmpDir}/${SLURM_ARRAY_TASK_ID}_r200 | samtools view -bS - > ${outputDir}/${SLURM_ARRAY_TASK_ID}_r200.bam
cat ${tmpDir}/${SLURM_ARRAY_TASK_ID}_r2_header.sam ${tmpDir}/${SLURM_ARRAY_TASK_ID}_r201 | samtools view -bS - > ${outputDir}/${SLURM_ARRAY_TASK_ID}_r201.bam

#Peak calling on self-pseudoreplicates
echo "Calling peaks for pseudoreplicate1"
macs2 callpeak -t ${outputDir}/${SLURM_ARRAY_TASK_ID}_r100.bam -c ${SLURM_ARRAY_TASK_ID}.Input.mq20.sorted.bam -f BAMPE -g dm -n $macsDir/${SLURM_ARRAY_TASK_ID}.r100_pr -B -p 1e-3 2> $macsDir/${SLURM_ARRAY_TASK_ID}.r100_pr_macs2.log

echo "Calling peaks for pseudoreplicate2"
macs2 callpeak -t ${outputDir}/${SLURM_ARRAY_TASK_ID}_r101.bam -c ${SLURM_ARRAY_TASK_ID}.Input.mq20.sorted.bam -f BAMPE -g dm -n $macsDir/${SLURM_ARRAY_TASK_ID}.r101_pr -B -p 1e-3 2> $macsDir/${SLURM_ARRAY_TASK_ID}.r101_pr_macs2.log

echo "Calling peaks for pseudoreplicate1"
macs2 callpeak -t ${outputDir}/${SLURM_ARRAY_TASK_ID}_r200.bam -c ${SLURM_ARRAY_TASK_ID}.Input.mq20.sorted.bam -f BAMPE -g dm -n $macsDir/${SLURM_ARRAY_TASK_ID}.r200_pr -B -p 1e-3 2> $macsDir/${SLURM_ARRAY_TASK_ID}.r200_pr_macs2.log

echo "Calling peaks for pseudoreplicate2"
macs2 callpeak -t ${outputDir}/${SLURM_ARRAY_TASK_ID}_r201.bam -c ${SLURM_ARRAY_TASK_ID}.Input.mq20.sorted.bam -f BAMPE -g dm -n $macsDir/${SLURM_ARRAY_TASK_ID}.r201_pr -B -p 1e-3 2> $macsDir/${SLURM_ARRAY_TASK_ID}.r201_pr_macs2.log

#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr $macsDir/${SLURM_ARRAY_TASK_ID}.r100_pr_peaks.narrowPeak | head -n 100000 > $macsDir/${SLURM_ARRAY_TASK_ID}.r100_pr_sorted.narrowPeak
sort -k8,8nr $macsDir/${SLURM_ARRAY_TASK_ID}.r101_pr_peaks.narrowPeak | head -n 100000 > $macsDir/${SLURM_ARRAY_TASK_ID}.r101_pr_sorted.narrowPeak
sort -k8,8nr $macsDir/${SLURM_ARRAY_TASK_ID}.r200_pr_peaks.narrowPeak | head -n 100000 > $macsDir/${SLURM_ARRAY_TASK_ID}.r200_pr_sorted.narrowPeak
sort -k8,8nr $macsDir/${SLURM_ARRAY_TASK_ID}.r201_pr_peaks.narrowPeak | head -n 100000 > $macsDir/${SLURM_ARRAY_TASK_ID}.r201_pr_sorted.narrowPeak

#Independent replicate IDR
echo "Running IDR on pseudoreplicates..."
idr --samples $macsDir/${SLURM_ARRAY_TASK_ID}.r100_pr_sorted.narrowPeak $macsDir/${SLURM_ARRAY_TASK_ID}.r101_pr_sorted.narrowPeak --input-file-type narrowPeak --output-file ${SLURM_ARRAY_TASK_ID}.r1_self_pseudorep-idr_p1e-3 --rank p.value --plot --log-output-file ${SLURM_ARRAY_TASK_ID}.r1.self.pr.idr.p1e-3.log

idr --samples $macsDir/${SLURM_ARRAY_TASK_ID}.r200_pr_sorted.narrowPeak $macsDir/${SLURM_ARRAY_TASK_ID}.r201_pr_sorted.narrowPeak --input-file-type narrowPeak --output-file ${SLURM_ARRAY_TASK_ID}.r2_self_pseudorep-idr_p1e-3 --rank p.value --plot --log-output-file ${SLURM_ARRAY_TASK_ID}.r2.self.pr.idr.p1e-3.log

echo "Done!"
