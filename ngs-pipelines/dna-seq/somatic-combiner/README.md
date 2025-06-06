# Somatic-Combiner
This Python script merges somatic VCF files from different callers (e.g., Mutect2, MuSE, Strelka, VarScan) into a single VCF file containing variants that are present in at least two of the provided callers. It also filters the input VCFs, processes variant information, and outputs a merged VCF with the relevant data for each variant.

### Requirements
- `Python 3.x`
- Required libraries: `os`, `argparse`, `collections`, `datetime`, `tqdm`, `numpy`

### Arguments
- `--mutect2`. Path to Mutect2 VCF file (optional)
- `--muse`. Path to MuSE VCF file (optional)
- `--strelka-snvs`. Path to Strelka SNVs VCF file (optional)
- `--strelka-indels`. Path to Strelka Indels VCF file (optional)
- `--varscan-snvs`. Path to VarScan SNVs VCF file (optional)
- `--varscan-indels`. Path to VarScan Indels VCF file (optional)
- `--output_dir`. Directory where the filtered and merged VCF will be saved
- `--id`. Name for the merged VCF output (without path)

### Features
- **VCF Filtering**. Filters the input VCF files to exclude irrelevant or incomplete variants based on the caller.
- **Merging Variants**. Merges variants that are present in at least two different callers.
- **Variant Consensus**. Calculates consensus genotypes (GT) and other fields like allele frequency (AF), depth (DP), genotype quality (GQ), etc.
- **Custom Headers**. Adds custom headers for the merged VCF, including information on callers and variant characteristics.

### Output
The script generates a merged VCF file containing only the variants that are present in ≥2 callers, with combined information from the different callers (GT, AF, DP, etc.).
The output file will be saved in the specified output_dir with the name provided in id.

### Example
```bash
python somatic-combiner/somatic-combiner.py --mutect2 mutect2.vcf.gz --muse muse.vcf.gz --strelka-indels strelka2.indels.vcf.gz --strelka-snvs strelka2.snvs.vcf.gz --varscan-indels varscan.indels.vcf.gz --varscan-snvs varscan.snvs.vcf.gz --output_dir output/ --id somatic-combiner
```

### License
This script is released under the MIT License. See the LICENSE file for more details.
