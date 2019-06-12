#!/bin/bash
# This script is intended to streamline generation of blobplots
# for a given reference sequence and set of paired reads,
# it checks it generates plots for the whole dataset,
# the reads mapping to the reference, and the reads not
# mapping to the reference
#
# 2017-11-24
# Nick Waters
usage() { echo "Usage: $0   [-a assembly.fasta] [-t num_threads] [-d blast_database_name ] [-o ./output/dir/ ] " 1>&2; exit 1; }

set -euo pipefail
IFS=$'\n\t'
while getopts ":c:t:d:o:h" opt; do
    case $opt in
	h) usage
	   ;;
	c) ASSEMBLY="$OPTARG"
	 ;;
	t) THREADS="$OPTARG"
	   ;;
	d) BDB="$OPTARG"
	   ;;
	o) OUTDIR="$OPTARG"
	   ;;
	\?) echo "Invalid option -$OPTARG" >&2
	    ;;
    esac
done
#  lets ensure that the environment variable BLASTDB has been set, as per https://www.ncbi.nlm.nih.gov/books/NBK52640/
echo "BLAST database search path: $BLASTDB"
python --version
echo $SHELL

blobtools create -h
# ensure that metaspades. samtools, and blastn are present
# we have set -e., so if any of these calls fail, the script fails

# samtools --version | grep samtools 
# blastn -version | grep "Package: blast" 

LOG=${OUTDIR}blob.log
  
mkdir $OUTDIR
PATH="/home/nw42839/GitHub/blobtools/ncbi-blast-2.8.1+/bin:$PATH"
HITS_ALL=${OUTDIR}all_blastn_nt.out
echo "Running blast for all"
blastn -db $BDB -query ${ASSEMBLY} -out $HITS_ALL -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -num_threads $THREADS 

echo "running blobtools cmds"
BLOB_ALL=${OUTDIR}blob_all/
mkdir $BLOB_ALL
blobtools create -i ${ASSEMBLY} -y spades -o $BLOB_ALL -t $HITS_ALL  

echo "Generating plots at the Species level" 
blobtools blobplot -i ${BLOB_ALL}blobDB.json -o ${BLOB_ALL}all --rank species
