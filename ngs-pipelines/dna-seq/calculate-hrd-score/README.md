# HRD Score Computational Pipeline
A standardized pipeline for calculating Homologous Recombination Deficiency (HRD) scores from paired tumor/normal sequencing data.

## 📊 HRD Score Components
The HRD score is calculated through integration of three genomic instability metrics:
- **HRD-LOH (Loss of Heterozygosity in HR-deficient tumors)**: Count of LOH regions >15 Mb that do not span an entire chromosome.  
- **LST (Large-Scale Transitions)**: Count of chromosomal breaks >10 Mb, separated by ≤3 Mb.  
- **Telomeric AI (Telomeric Allelic Imbalance)**: Count of allelic imbalances (AI) that extend to the telomeric end of a chromosome.  

**Final HRD Score** = HRD-LOH + LST + Telomeric AI

## 🧬 Input Requirements  
### 🔹 Required Data  
- **Paired-end sequencing data** (tumor/normal pairs).  
- **Supported data types**:  
  - Whole Exome Sequencing (**WES**)  
  - Whole Genome Sequencing (**WGS**)
  - RNA-Seq

### 📂 File Naming Convention  
- **Tumor samples**: `{ID}_C_1.fastq.gz`, `{ID}_C_2.fastq.gz`  
- **Normal samples**: `{ID}_N_1.fastq.gz`, `{ID}_N_2.fastq.gz`  

### 📝 Example Input Format  
The input must be provided as a structured table with the following columns:  

| ID  | Type | FASTQ Filename            |
|-----|------|---------------------------|
| LC_1 | C    | LC_1_C_1.fastq.gz         |
| LC_1 | C    | LC_1_C_2.fastq.gz         |
| LC_1 | N    | LC_1_N_1.fastq.gz         |
| LC_1 | N    | LC_1_N_2.fastq.gz         |
| LC_2 | C    | LC_2_C_1.fastq.gz         |
| LC_2 | C    | LC_2_C_2.fastq.gz         |
| LC_2 | N    | LC_2_N_1.fastq.gz         |
| LC_2 | N    | LC_2_N_2.fastq.gz         |

**Legend:**  
- **ID**: Sample identifier.  
- **Type**: `"C"` for tumor (cancer) samples, `"N"` for normal samples.  
- **FASTQ Filename**: Corresponding FASTQ file name for sequencing data.

## 🚀 Pipeline Execution
**1. Build Docker Image**
```
docker build -t hrd-calculator .
```
**2. Run Analysis**
```
docker run -it --rm \
  -v /path/to/data:/data \
  -v /path/to/results:/results \
  hrd-calculator \
  --cores 8 \
  --config samples="/data/samples.csv" \
  --use-conda
```

## 📄 Output File Format  
### 📁 Output Filename: `hrd_score.txt`  
The output file contains HRD-related scores for each sample in a tab-separated format with the following columns:  

| Sample_ID | HRD_LOH | LST  | Telomeric_AI | HRD_sum |
|-----------|----------|------|--------------|---------|
| LC_1      | 58       | 19   | 22           | 17      |
| LC_2      | 47       | 15   | 18           | 14      |

### 🔍 Column Descriptions  
- **Sample_ID**: Unique identifier for each sample.  
- **HRD_LOH**: Count of LOH regions >15 Mb.  
- **LST (Large-Scale Transitions)**: Number of chromosomal breaks >10 Mb, separated by ≤3 Mb.  
- **Telomeric_AI (Telomeric Allelic Imbalance)**: Number of allelic imbalances extending to the telomeric end.  
- **HRD_sum**: Overall homologous recombination deficiency (HRD) score.  
