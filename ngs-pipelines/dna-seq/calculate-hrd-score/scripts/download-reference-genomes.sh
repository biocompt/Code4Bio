#!/bin/sh

# GENCODE GRCh38.p14 primary assembly (FASTA)
if [ ! -f "utilities/reference-genomes/GRCh38.primary_assembly.genome.fa" ]; then
    echo "Downloading GRCh38 FASTA..."
    curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz \
         --create-dirs -o utilities/reference-genomes/GRCh38.primary_assembly.genome.fa.gz
    bgzip -d utilities/reference-genomes/GRCh38.primary_assembly.genome.fa.gz
else
    echo "GRCh38 FASTA already exists. Skipping..."
fi


# GATK Resources Bundle for alignment
if [ ! -f "utilities/reference-genomes/GRCh38.known_indels.vcf.gz" ]; then
    echo "Downloading known indels VCF..."
    curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz \
         --create-dirs -o utilities/reference-genomes/GRCh38.known_indels.vcf.gz
else
    echo "Known indels VCF already exists. Skipping..."
fi

if [ ! -f "utilities/reference-genomes/GRCh38.known_indels.vcf.gz.tbi" ]; then
    echo "Downloading known indels index..."
    curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi \
         --create-dirs -o utilities/reference-genomes/GRCh38.known_indels.vcf.gz.tbi
else
    echo "Known indels index already exists. Skipping..."
fi

if [ ! -f "utilities/reference-genomes/GRCh38.Mills_1000G_standard_indels.vcf.gz" ]; then
    echo "Downloading Mills & 1000G indels..."
    curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
         --create-dirs -o utilities/reference-genomes/GRCh38.Mills_1000G_standard_indels.vcf.gz
else
    echo "Mills & 1000G indels VCF already exists. Skipping..."
fi

if [ ! -f "utilities/reference-genomes/GRCh38.Mills_1000G_standard_indels.vcf.gz.tbi" ]; then
    echo "Downloading Mills & 1000G index..."
    curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
         --create-dirs -o utilities/reference-genomes/GRCh38.Mills_1000G_standard_indels.vcf.gz.tbi
else
    echo "Mills & 1000G index already exists. Skipping..."
fi