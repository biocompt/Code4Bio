# Solving the Initial Doubts
## 1. Which Genome Should I Use? Where Can I Find It?
When starting an RNA-Seq analysis, one of the first challenges is selecting the reference genome and its gene annotations. These files (typically in FASTA format for the genome sequence and GTF/GFF3 format for annotations) are foundational for aligning sequencing reads and quantifying gene expression. Below are the most widely used databases and their key differences:

- **GENCODE** is often the top choice for human and mouse studies. Unlike other databases, GENCODE combines manual curation by expert biologists with computational pipelines to ensure high accuracy. It includes not only protein-coding genes but also non-coding RNAs, pseudogenes, and regulatory elements. This makes it ideal for RNA-Seq analyses where detecting alternative splicing or lncRNAs is critical. Another advantage is its compatibility with popular alignment tools like STAR or Salmon, as it standardizes chromosome naming (e.g., "chr1" instead of "1").
- **Ensembl** shares many annotations with GENCODE due to their collaborative relationship but covers a broader range of species. While Ensembl is comprehensive, its automated pipelines may include more speculative gene models compared to GENCODE’s curated approach.
- **NCBI RefSeq** focuses on high-confidence, experimentally validated transcripts. It’s smaller in scope than GENCODE and better suited for projects prioritizing well-characterized genes. However, it may miss novel transcripts or poorly expressed genes.
- **UCSC Genome Browser** provides genome assemblies and annotation tracks, but its files often lack the detail required for modern RNA-Seq workflows. Historically, UCSC’s "hg19" (GRCh37) was widely used, but its legacy status makes it less relevant for new projects.

### **My Choice: GENCODE**
After evaluating multiple databases, I consistently use GENCODE for RNA-Seq analyses. The decision hinges on its unparalleled balance of accuracy, biological relevance, and practical usability. GENCODE stands out for its manually curated annotations, a feature lacking in purely computational pipelines like Ensembl’s automated tracks. Another advantage is its transcript diversity. Unlike RefSeq—which prioritizes "representative" isoforms—GENCODE catalogs all credible splice variants. This is critical for studies in tissues like the liver or immune system, where alternative splicing is rampant.

GENCODE also simplifies compatibility with bioinformatics tools. Its GTF files include standardized chromosome names (e.g., chr1 instead of 1), which align seamlessly with tools like STAR or Salmon. In contrast, UCSC’s mixed naming conventions, wasting hours of troubleshooting.

Finally, GENCODE’s regular updates integrate the latest research, such as newly discovered miRNAs or corrected pseudogene annotations. While legacy databases like UCSC’s freeze outdated data, GENCODE evolves with the field, ensuring reproducibility in peer-reviewed studies.

## 2. GRCh37 vs. GRCh38: Why Version Matters
The human genome has two primary assembly versions: **GRCh37 (hg19)** and **GRCh38 (hg38)**. While GRCh37 dominated for over a decade due to its stability and widespread adoption, it has critical limitations. For example, it contains assembly errors in regions like centromeres and telomeres, and it omits many genes discovered after 2009.

In contrast, GRCh38 (specifically the latest patches like GRCh38.p14) addresses these gaps and errors. It includes over 400 corrected regions, alternative haplotypes for complex loci (e.g., the major histocompatibility complex), and updated gene annotations reflecting recent research. For RNA-Seq, this means fewer alignment artifacts and better detection of immune-related genes, cancer biomarkers, and non-coding RNAs. **Unless you’re working with legacy data tied to older tools, GRCh38 is the clear choice for accuracy and future-proofing.**

## 3. Understanding GENCODE File Types: What to Choose?
GENCODE offers multiple file types, which can be overwhelming for new users. For bulk RNA-Seq, the choice simplifies significantly by prioritizing essential data and avoiding redundancy.

For **annotation files (GTF/GFF3), the most practical option is the Basic gene annotation (CHR)**. This set includes only experimentally validated "basic" annotations on reference chromosomes (e.g., chr1 to chr22, X, Y, MT), excluding ambiguous regions like scaffolds or alternate haplotypes. For example, if you’re analyzing a gene like BRCA1, this file retains only the best-supported isoforms, reducing false positives during quantification. In contrast, the Comprehensive annotation includes all predictions (even unvalidated ones), which can introduce noise in standard projects.

For the **reference genome (FASTA), select the Genome sequence, primary assembly (PRI)**. This file contains the sequence of primary chromosomes and unplaced scaffolds, covering 99% of alignment needs. Options like Genome sequence (ALL) add patches and alternate haplotypes (e.g., HLA system variants), which are useful for population or immunogenomics studies but unnecessary for most gene expression analyses.

```bash
# Genome (FASTA)  
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz  

# Annotations (GTF)  
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz
```

## 4. Quantification Tools: RSEM vs. Salmon
Once you have your reference genome and annotations, the next step is selecting a quantification tool to estimate gene/transcript expression. While htseq-count was historically popular, modern tools like Salmon and RSEM offer advantages such as bias correction and faster processing, while directly generating gene-level counts. Below is a practical comparison:

#### RSEM (RNA-Seq by Expectation-Maximization)
RSEM is an alignment-based tool that works with aligners like Bowtie2 or STAR. It quantifies expression by resolving ambiguously mapped reads across isoforms using an expectation-maximization (EM) algorithm. While accurate, RSEM requires a two-step workflow:

1. Align reads to the genome (e.g., generating a BAM file).
2. Run RSEM on the BAM file to estimate counts.

This makes it computationally heavy, especially for large datasets. Additionally, RSEM does not natively correct for GC bias or sequence-specific biases, which can affect accuracy in low-expression genes.

#### Salmon
Salmon is an alignment-free (or quasi-mapping) tool that skips the alignment step entirely. It uses k-mer indexing to map reads directly to the transcriptome, making it 10–20x faster than RSEM. Key advantages include:

- **Bias correction**: Salmon adjusts for GC content, sequence bias, and positional fragment bias, improving accuracy in quantification.
- **Dual outputs**: It provides both transcript-level counts (ideal for isoform analysis) and gene-level counts (via tximport post-processing).
- **Resource efficiency**: Requires less RAM and disk space compared to RSEM.

### My Choice: Salmon
For bulk RNA-Seq, Salmon is the preferred choice for most workflows:
1. **Speed and scalability**: Critical for classrooms or labs with limited computational resources.
2. **Bias-aware models**: Reduces false expression estimates in genes with extreme GC content.
3. **Simpler workflow**: No intermediate BAM files to manage.
