#!/bin/sh
# *******************************************
# Script to perform ClinSV to call CNV from WGS
# *******************************************

set -e
ulimit -s 10240

module load clinsv/1.1

workdir_remote="/lscratch/$SLURM_JOB_ID/"
workdir='/path/to/working/dir'
bamdir='/path/to/bam/file/dir'
bamfile='name/bamfile'
refdata_path=/usr/local/apps/clinsv/ref_hg38/clinsv/refdata-b38

mkdir -p $workdir_remote

cd $workdir

# Run ClinSV
singularity run --bind $refdata_path:/app/ref-data/refdata-b38 \
	    --bind $project_folder:$workdir \
	    --bind $input_path:$bamdir \
	    /usr/local/apps/clinsv/1.1/libexec/1.1.sif \
        /app/clinsv/bin/clinsv \
	    -r all \
	    -p $workdir  \
	    -i $bamdir/$bamfile \
	    -ref /app/ref-data/refdata-b38 \
	    -eval yes
