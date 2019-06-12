#!/bin/bash
# This script is intended to streamline generation of blobplots
# for a given reference sequence and set of paired reads,
# it checks it generates plots for the whole dataset,
# the reads mapping to the reference, and the reads not
# mapping to the reference
#
# 2017-11-24
# Nick Waters
usage() { echo "Usage: $0 [-r reference.fasta] [-F forward.fq ] [-R reverse.fastq] [-d blast_database_name ] [-o ./output/dir/ ] " 1>&2; exit 1; }
source  activate blob2
set -euo pipefail
IFS=$'\n\t'
while getopts ":r:F:R:o:d:t:h" opt; do
    case $opt in
	h) usage
	   ;;	
	r) REF="$OPTARG"
	   ;;
	F) READS1="$OPTARG"
	 ;;
	R) READS2="$OPTARG"
	   ;;
	d) BDB="$OPTARG"
	   ;;
	o) OUTDIR="$OPTARG"
	   ;;
	t) THREADS="$OPTARG"
	   ;;
	\?) echo "Invalid option -$OPTARG" >&2
	    ;;
    esac
done
#  lets ensure that the environment variable BLASTDB has been set, as per https://www.ncbi.nlm.nih.gov/books/NBK52640/
echo "BLAST database search path: $BLASTDB"

# ensure that metaspades. samtools, and blastn are present
# we have set -e., so if any of these calls fail, the script fails
samtools --version | grep samtools 
blastn -version | grep "Package: blast" 
echo "checking bwa"
#bwa |& grep "Version:"
echo "setting log"
LOG=${OUTDIR}blob.log

echo "making outdir"
if [ -d "$OUTDIR" ]
then
   echo "using existing assembly"
else
    echo "creating output directory"
    mkdir $OUTDIR
fi

# echo "index and map reads to reference" 
# bwa index $REF >> $LOG 2>&1
# BWASAM=${OUTDIR}mapped_reads.sam
# bwa mem -t $THREADS $REF $READS1 $READS2 > $BWASAM

echo " extract reads mapping" 
READS1_MAPPED=${OUTDIR}mapped_1.fastq
READS2_MAPPED=${OUTDIR}mapped_2.fastq
READS1_UNMAPPED=${OUTDIR}unmapped_1.fastq
READS2_UNMAPPED=${OUTDIR}unmapped_2.fastq
# samtools flagstat $BWASAM >> $LOG 2>&1
# samtools fastq -F 12 $BWASAM -1 $READS1_MAPPED -2 $READS2_MAPPED >> $LOG 2>&1
# samtools fastq -f 12 $BWASAM -1 $READS1_UNMAPPED -2 $READS2_UNMAPPED >> $LOG 2>&1

echo " assembling raw reads and subset reads with metaspades"
SPADES_ALL=${OUTDIR}metaspades_all/
# SPADES_MAPPED=${OUTDIR}metaspades_mapped/
# SPADES_UNMAPPED=${OUTDIR}metaspades_unmapped/

if [ -d "$SPADES_ALL" ]
then
    echo "using prior assembly"
else
    mkdir $SPADES_ALL
    skesa --cores $THREADS --fastq $READS1 $READS2 --contigs_out $SPADES_ALL/contigs.fasta
fi

# mkdir $SPADES_MAPPED
# mkdir $SPADES_UNMAPPED
#skesa --cores $THREADS --fastq $READS1_MAPPED $READS2_MAPPED --contigs_out $SPADES_MAPPED/contigs.fasta
#skesa --cores $THREADS --fastq $READS1_UNMAPPED $READS2_UNMAPPED --contigs_out $SPADES_UNMAPPED/contigs.fasta

echo "Running blast for mapped"
HITS_MAPPED=${OUTDIR}mapped_blastn_nt.out
HITS_UNMAPPED=${OUTDIR}unmapped_blastn_nt.out
HITS_ALL=${OUTDIR}all_blastn_nt.out
# blastn -db $BDB -query ${SPADES_MAPPED}contigs.fasta -out $HITS_MAPPED -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -num_threads 3
# echo "Running blast for unmapped"
# blastn -db $BDB -query ${SPADES_UNMAPPED}contigs.fasta -out $HITS_UNMAPPED -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -num_threads $THREADS
echo "Running blast for all"
blastn -db $BDB -query ${SPADES_ALL}contigs.fasta -out $HITS_ALL -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -num_threads $THREADS

echo "running blobtools cmds"
BLOB_ALL=${OUTDIR}blob_all/
BLOB_MAPPED=${OUTDIR}blob_mapped/
BLOB_UNMAPPED=${OUTDIR}blob_unmapped/
mkdir $BLOB_ALL
mkdir $BLOB_MAPPED
mkdir $BLOB_UNMAPPED

# remap to assembly
bwa index ${SPADES_ALL}contigs.fasta
# bwa index ${SPADES_MAPPED}contigs.fasta
# bwa index ${SPADES_UNMAPPED}contigs.fasta

# 
BWASAM_ALL=${BLOB_ALL}mapping.sam
# BWASAM_ALL=${BLOB_MAPPED}mapping.sam
# BWASAM_ALL=${BLOB_UNMAPPED}mapping.sam
bwa mem -t $THREADS ${SPADES_ALL}contigs.fasta $READS1 $READS2 > $BWASAM_ALL
# bwa mem ${SPADES_MAPPED}contigs.fasta $READS1 $READS2 > $BWASAM_MAPPED
# bwa mem ${SPADES_UNMAPPED}contigs.fasta $READS1 $READS2 > $BWASAM_UNMAPPED

# blobtools create -i ${SPADES_MAPPED}contigs.fasta  -o $BLOB_MAPPED -t $HITS_MAPPED  -s $BWASAM_MAPPED 
# blobtools create -i ${SPADES_UNMAPPED}contigs.fasta  -o $BLOB_UNMAPPED -t $HITS_UNMAPPED  -s $BWASAM_MAPPED 
singularity exec ./blobtools.simg /blobtools/blobtools create -i ${SPADES_ALL}contigs.fasta  -o $BLOB_ALL -t $HITS_ALL  -s $BWASAM_ALL 

echo "Generating plots at the Species level" 
# blobtools blobplot -i ${BLOB_MAPPED}blobDB.json -o ${BLOB_MAPPED}mapped --rank species
# blobtools blobplot -i ${BLOB_UNMAPPED}blobDB.json -o ${BLOB_UNMAPPED}unmapped --rank species
singularity exec ./blobtools.simg /blobtools/blobtools blobplot -i ${BLOB_ALL}blobDB.json -o ${BLOB_ALL}all --rank species
