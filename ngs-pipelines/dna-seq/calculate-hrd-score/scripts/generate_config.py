import yaml
import os
import pandas as pd
import argparse

# 🎯 Configurar argumentos
parser = argparse.ArgumentParser(description="Genera snakeconfig.yaml para Snakemake")
parser.add_argument("--csv_file", type=str, required=True, help="Archivo CSV con muestras")
parser.add_argument("--data_dir", type=str, required=True, help="Directorio donde están los FASTQs")
parser.add_argument("--user_slurm", type=str, required=True, help="SLURM user account (ejemplo: bihcxcpm)")
parser.add_argument("--config_file", type=str, default="utilities/snakeconfig.yaml", help="Archivo de configuración de Snakemake")
args = parser.parse_args()

# 📌 Función para leer el CSV y estructurar las muestras
def get_samples_from_csv(csv_file, data_dir):
    df = pd.read_csv(csv_file)
    samples = {}
    for _, row in df.iterrows():
        sample_id = row["ID"]
        category = row["Type"]
        fastq_file = row["FASTQ"]

        # Verifica si el archivo existe en el directorio data
        fastq_path = os.path.join(data_dir, fastq_file)
        if not os.path.exists(fastq_path):
            print(f"⚠️ Advertencia: No se encontró {fastq_path}")

        # Organizar en la estructura esperada
        if sample_id not in samples:
            samples[sample_id] = {"C": [], "N": []}
        if category in ["C", "N"]:
            samples[sample_id][category].append(fastq_path)
    return samples

# 📦 Cargar muestras desde CSV
samples = get_samples_from_csv(args.csv_file, args.data_dir)

# 📌 Estructura base para snakeconfig.yaml
config = {
    "samples": samples,
    "ref": {
        "ref_fa": "utilities/reference-genomes/GRCh38.primary_assembly.genome.fa",
        "known_indels": "utilities/reference-genomes/GRCh38.known_indels.vcf.gz",
        "standard_indels": "utilities/reference-genomes/GRCh38.Mills_1000G_standard_indels.vcf.gz"
    },
    "resources": {
        # Aquí establecemos el usuario SLURM pasado por línea de comandos
        "user_cesga": args.user_slurm
    },
    "scripts": {
        "hrd-score": "utilities/scripts/run_sequenza.R"
    },
    "software": {
        "fastp": "utilities/software/fastp.yaml",
        "bwa": "utilities/software/bwa.yaml",
        "sequenza": "utilities/software/sequenza.yaml",
        "hrd-score": "utilities/software/hrd-score.yaml"
    }
}

# 📝 Guardar el YAML resultante
with open(args.config_file, "w") as f:
    yaml.dump(config, f, default_flow_style=False)

print(f"✅ Archivo {args.config_file} actualizado con {len(samples)} muestras. SLURM user: {args.user_slurm}")