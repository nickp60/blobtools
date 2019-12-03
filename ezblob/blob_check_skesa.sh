#!/bin/bash
# This script is intended to streamline generation of blobplots
# for a given reference sequence and set of paired reads,
# it checks it generates plots for the whole dataset,
# the reads mapping to the reference, and the reads not
# mapping to the reference
#
# 2017-11-24
# Nick Waters
usage() { echo -e "Usage: $0 [-r reference.fasta] [-F forward.fq ] [-R reverse.fastq] [-A assembly.fasta] [-d blast_database_name ] [-o ./output/dir/ ] [-t #threads] [-m #memory] \n you can submit either single-end reads with -F, paired-end reads with -F and -R, or a premade assembly with -A " 1>&2; exit 1; }
set -euox pipefail
IFS=$'\n\t'
READS2=""
ASSEMBLY=""
while getopts ":r:F:R:o:d:t:m:A:h" opt; do
    case $opt in
	h) usage
	   ;;
	r) REF="$OPTARG"
	   echo "reference: $REF"
	   ;;
	F) READS1="$OPTARG"
	   echo "Forward reads $READS1"
	   ;;
	R) READS2="$OPTARG"
	   echo "Reverse reads $READS2"
	   ;;
	A) ASSEMBLY="$OPTARG"
	   echo "Assembly $ASSEMBLY"
	   ;;
	d) BDB="$OPTARG"
	   echo "BLAST database $BDB"
	   ;;
	o) OUTDIR="$OPTARG"
	   echo "Output directory: $OUTDIR"
	   ;;
	t) THREADS="$OPTARG"
	   echo "Threads: $THREADS"
	   ;;
	m) MEMORY="$OPTARG"
	   echo "MEMORY: $MEMORY"
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

# # echo "index and map reads to reference"
# bwa index $REF >> $LOG 2>&1
# BWASAM=${OUTDIR}mapped_reads.sam
# bwa mem -t $THREADS $REF $READS1 $READS2 > $BWASAM

# echo " extract reads mapping"
# READS1_UNMAPPED=${OUTDIR}unmapped_1.fastq
# READS2_UNMAPPED=${OUTDIR}unmapped_2.fastq
# samtools flagstat $BWASAM >> $LOG 2>&1
# samtools fastq -f 12 $BWASAM -1 $READS1_UNMAPPED -2 $READS2_UNMAPPED >> $LOG 2>&1
SPADES_ALL=${OUTDIR}assembly_all/
mkdir $SPADES_ALL
if [ "$ASSEMBLY" != "" ]
then
    echo "Assembly provided; skipping SKESA"
    cp $ASSEMBLY $SPADES_ALL/contigs.fasta
else
    echo " assembling raw reads and subset reads with SKESA"
    #SPADES_UNMAPPED=${OUTDIR}assembly_unmapped/

    skesa --cores $THREADS --memory $MEMORY --fastq $READS1 $READS2 --contigs_out $SPADES_ALL/contigs.fasta
    #mkdir $SPADES_UNMAPPED
    #skesa --cores $THREADS --fastq $READS1_UNMAPPED $READS2_UNMAPPED --contigs_out $SPADES_UNMAPPED/contigs.fasta
fi

echo "Running blast for mapped"
HITS_UNMAPPED=${OUTDIR}unmapped_blastn_nt.out
HITS_ALL=${OUTDIR}all_blastn_nt.out

blastn -db /db/$BDB -query ${SPADES_ALL}contigs.fasta -out $HITS_ALL -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -num_threads $THREADS
#echo "Running blast for unmapped"
#blastn -db /db/$BDB -query ${SPADES_UNMAPPED}contigs.fasta -out $HITS_UNMAPPED -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -num_threads $THREADS

echo "running blobtools cmds"
BLOB_ALL=${OUTDIR}blob_all/
#BLOB_UNMAPPED=${OUTDIR}blob_unmapped/
mkdir $BLOB_ALL
#mkdir $BLOB_UNMAPPED

# remap to assembly
bwa index ${SPADES_ALL}contigs.fasta
# bwa index ${SPADES_MAPPED}contigs.fasta
# bwa index ${SPADES_UNMAPPED}contigs.fasta

#
if [ "$ASSEMBLY" == '' ]
then

    BWASAM_ALL=${BLOB_ALL}mapping.bam
    BWASAM_ALL_SORT=${BLOB_ALL}mapping.sorted.bam

    # BWASAM_UNMAPPED=${BLOB_UNMAPPED}mapping.bam
    # BWASAM_UNMAPPED_SORT=${BLOB_UNMAPPED}mapping.sorted.bam

    bwa mem -t $THREADS ${SPADES_ALL}contigs.fasta $READS1 $READS2 | samtools view -hb - > $BWASAM_ALL
    samtools sort -T /tmp/aln.sorted -o $BWASAM_ALL_SORT $BWASAM_ALL
    samtools index $BWASAM_ALL_SORT

    # bwa mem -t $THREADS ${SPADES_UNMAPPED}contigs.fasta $READS1_UNMAPPED $READS2_UNMAPPED | samtools view -hb - > $BWASAM_UNMAPPED
    # samtools sort -T /tmp/aln.sorted -o $BWASAM_UNMAPPED_SORT $BWASAM_UNMAPPED
    # samtools index $BWASAM_UNMAPPED_SORT



    #/blobtools/blobtools create -i ${SPADES_UNMAPPED}contigs.fasta  -o $BLOB_UNMAPPED -t $HITS_UNMAPPED  -b $BWASAM_UNMAPPED_SORT
    /blobtools/blobtools create -i ${SPADES_ALL}contigs.fasta  -o $BLOB_ALL -t $HITS_ALL  -b $BWASAM_ALL_SORT

else
    # the fasta type from SKESA is similar enough to velvet: https://github.com/nickp60/blobtools/blob/runner_script/lib/BtIO.py
    /blobtools/blobtools create -i ${SPADES_ALL}contigs.fasta -y velvet  -o $BLOB_ALL -t $HITS_ALL
fi

echo "Generating plots at the Species level"
#/blobtools/blobtools plot -i ${BLOB_UNMAPPED}blobDB.json -o ${BLOB_UNMAPPED}unmapped --rank species
/blobtools/blobtools plot -i ${BLOB_ALL}blobDB.json -o ${BLOB_ALL}all --rank species
