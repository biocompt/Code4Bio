#!/usr/bin/env python3
"""
filtrar_strelka_vcf.py

Este script carga un VCF de Strelka y, sin cambiar el esquema de columnas
de VCF, reescribe únicamente los subcampos que nos interesan:

  • En INFO: conservar solo “SOMATIC” (flag)
  • En FORMAT: conservar “DP”, “MQ” y “MQ0”:
      - DP se extrae de FORMAT original (por muestra)
      - MQ y MQ0 se extraen de INFO (valores por variante) y se incluyen
        en cada muestra bajo FORMAT (misma cifra para todas)

Imprime por pantalla todas las variantes en un VCF “reescrito”,
con las mismas columnas VCF (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, muestras…),
pero en INFO únicamente “SOMATIC” si existe (o “.” si no), 
y en FORMAT “DP:MQ:MQ0” con los valores correspondientes, evitando imprimir “[]” en ningún caso.

Uso:
    python filtrar_strelka_vcf.py -v ruta/al/Strelka.vcf.gz
"""

import argparse
import numpy as np
from cyvcf2 import VCF

def format_scalar(val):
    """
    Convierte valores escalares (incluyendo numpy scalars) a string.
    Si val es None o NaN, devuelve ".".
    """
    if val is None:
        return "."
    if isinstance(val, (np.integer, np.floating)):
        if np.isnan(val):
            return "."
        # np.integer → str(int(val)); np.floating → str(val) sin corchetes
        return str(int(val)) if isinstance(val, np.integer) else str(val)
    if isinstance(val, (int, float)):
        return str(val)
    # Si es lista/tupla/ndarray, unimos elementos con coma (sin corchetes)
    if isinstance(val, (list, tuple, np.ndarray)):
        elems = []
        for x in val:
            if x is None:
                elems.append(".")
            elif isinstance(x, (np.integer, np.floating)):
                if np.isnan(x):
                    elems.append(".")
                else:
                    elems.append(str(int(x)) if isinstance(x, np.integer) else str(x))
            else:
                elems.append(str(x))
        return ",".join(elems)
    # Cualquier otro tipo, convertir a str
    return str(val)

def main():
    parser = argparse.ArgumentParser(
        description="Reescribir un VCF de Strelka: INFO→solo SOMATIC, FORMAT→DP:MQ:MQ0"
    )
    parser.add_argument(
        "-v", "--vcf",
        required=True,
        help="Ruta al VCF de Strelka que se va a procesar"
    )
    args = parser.parse_args()
    vcf_path = args.vcf

    # Intentar abrir el VCF
    try:
        vcf_reader = VCF(vcf_path)
    except Exception as e:
        print(f"Error al abrir '{vcf_path}': {e}")
        return

    # Lista de muestras en orden
    sample_names = vcf_reader.samples

    # Imprimir línea de columnas VCF (sin metadatos ##)
    header_cols = [
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"
    ] + sample_names
    print("\t".join(header_cols))

    for record in vcf_reader:
        # 1) Campos VCF básicos
        chrom = record.CHROM
        pos   = record.POS
        id_field = record.ID if record.ID else "."
        ref  = record.REF
        alt  = ",".join(record.ALT) if record.ALT else "."
        qual = str(record.QUAL) if record.QUAL is not None else "."

        # 2) FILTER puede ser lista o cadena
        filt_field = record.FILTER
        if isinstance(filt_field, list):
            filt_str = ";".join(filt_field) if filt_field else "."
        elif isinstance(filt_field, str):
            filt_str = filt_field
        else:
            filt_str = "."

        # 3) INFO: conservar solo SOMATIC (flag)
        info_dict = dict(record.INFO)
        info_str = "SOMATIC" if "SOMATIC" in info_dict else "."

        # 4) Recuperar MQ y MQ0 desde INFO (valores por variante)
        mq_val = info_dict.get("MQ")
        mq0_val = info_dict.get("MQ0")
        mq_str = format_scalar(mq_val)
        mq0_str = format_scalar(mq0_val)

        # 5) FORMAT fijo: “DP:MQ:MQ0”
        fmt_str = "DP:MQ:MQ0"

        # 6) Para cada muestra: extraer DP de FORMAT original
        dp_array = record.format("DP")  # numpy array (n_muestras,) o None

        # Construir fila de salida
        row = [
            chrom, str(pos), id_field, ref, alt,
            qual, filt_str, info_str, fmt_str
        ]

        for i in range(len(sample_names)):
            # 6.a) DP de muestra
            dp_str = "."
            if dp_array is not None:
                try:
                    dp_str = format_scalar(dp_array[i])
                except Exception:
                    dp_str = "."

            # 6.b) MQ y MQ0 vienen de la variante (same for all samples)
            #     Ya formateados en mq_str y mq0_str

            sample_field = f"{dp_str}:{mq_str}:{mq0_str}"
            row.append(sample_field)

        # Imprimir línea completa de la variante
        print("\t".join(row))

if __name__ == "__main__":
    main()
