FROM continuumio/miniconda3
MAINTAINER Nick Waters <nickp60@gmail.com>
# set up dependencies
RUN conda install -c bioconda skesa blast samtools bwa
RUN conda install -c anaconda matplotlib docopt tqdm wget pyyaml git
RUN conda install -c bioconda pysam --update-deps

# Download blobtools, set up the nodes db
RUN git clone https://github.com/DRL/blobtools.git
WORKDIR blobtools
RUN ./blobtools -h
RUN wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
RUN tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
RUN ./blobtools nodesdb --nodes data/nodes.dmp --names data/names.dmp

# download taxonomy for blastdb
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz -P /db/
RUN tar zxf /db/taxdb.tar.gz -C /db/
RUN ls /db

# either download or add in the DB
#RUN update_blastdb.pl --passive --decompress ref_prok_rep_genomes
# For easier debugging build, download first.
# I had issues with a mysterious bad exit code after the DB was unpacked when
# running update_blast.pl from the container, and I dont have time to figure out why
ADD ref_prok_rep_* /db/
ENV BLASTDB=/db/

# add in the script for running assembly, mapping, etc
ADD blob_check_skesa.sh /bin/blob_check_skesa.sh
RUN chmod +x /bin/blob_check_skesa.sh
ENTRYPOINT [ "/bin/blob_check_skesa.sh" ]
