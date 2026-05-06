#!/bin/sh
# *******************************************
# Script to perform ClinSV to call CNV from WGS
# *******************************************

set -e
ulimit -s 10240

module load clinsv/1.1

workdir_remote="/lscratch/$SLURM_JOB_ID/${sample}"
workdir=/data/buirabornlt/Melanoma_family/WGS/Calling/ClinSV
bamdir=/data/buirabornlt/WGS_MelaNostrum/WGS/Bam_file/
refdata_path=/usr/local/apps/clinsv/ref_hg38/clinsv/refdata-b38

mkdir -p $workdir/$sample
mkdir -p $workdir_remote

cd $workdir/$sample

# Copy cram files and convert to bam
#rsync $cramfile ./
#samtools view -T /data/buirabornlt/Reference/Homo_sapiens_assembly38.fasta -b -o $sample.bam $cramfile
#samtools index $sample.bam

# Run ClinSV
# clinsv -p $workdir/$sample -i $sample.bam
singularity run --bind $refdata_path:/app/ref-data/refdata-b38 \
	    --bind $project_folder:$workdir/$sample \
	    --bind $input_path:$bamdir \
	    /usr/local/apps/clinsv/1.1/libexec/1.1.sif \
        /app/clinsv/bin/clinsv \
	    -r all \
	    -p $workdir/$sample  \
	    -i $bamdir/$sample.bam \
	    -ref /app/ref-data/refdata-b38 \
	    -eval yes
