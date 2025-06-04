#!/usr/bin/env python3
"""
filtrar_muse_vcf.py

Este script carga un VCF de MuSE y, sin cambiar el esquema de columnas
de VCF, reescribe únicamente los subcampos que nos interesan:

  • En INFO: conservar solo “SOMATIC” (flag)
  • En FORMAT: conservar “GT”, “DP” y “BQ” (por muestra)

En el caso de BQ, si el valor es un array (o lista), se unen los elementos
con comas y sin corchetes (por ejemplo, “30,25”).

Imprime por pantalla todas las variantes en un VCF “reescrito”,
con las mismas columnas VCF (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, muestras…),
pero en INFO únicamente “SOMATIC” si existe (o “.” si no), y en FORMAT “GT:DP:BQ”
con el genotipo, la profundidad de lectura y el valor de BQ por muestra.

Uso:
    python filtrar_muse_vcf.py -v ruta/al/MuSE.vcf.gz
"""

import argparse
import numpy as np
from cyvcf2 import VCF

def format_scalar_or_array(val):
    """
    Dado un valor que puede ser escalar o array/tupla/np.ndarray, 
    devuelve una cadena sin corchetes y con comas si hay varios elementos.
    """
    if val is None:
        return "."
    # Si es numpy scalar
    if isinstance(val, (np.integer, np.floating)):
        return str(int(val)) if isinstance(val, np.integer) else str(val)
    # Si es array numpy o lista/tupla
    if isinstance(val, (list, tuple, np.ndarray)):
        elementos = []
        for x in val:
            if x is None or (isinstance(x, (int, float, np.integer, np.floating)) and np.isnan(x)):
                elementos.append(".")
            else:
                # Convertir a entero si es entero numpy, si no, a str
                if isinstance(x, np.integer):
                    elementos.append(str(int(x)))
                else:
                    elementos.append(str(x))
        return ",".join(elementos)
    # Caso escalar estándar (int, float, str, etc.)
    return str(val)

def main():
    parser = argparse.ArgumentParser(
        description="Reescribir un VCF de MuSE: INFO→solo SOMATIC, FORMAT→GT,DP,BQ"
    )
    parser.add_argument(
        "-v", "--vcf",
        required=True,
        help="Ruta al VCF de MuSE que se va a procesar"
    )
    args = parser.parse_args()
    vcf_path = args.vcf

    # Intentamos abrir el VCF
    try:
        vcf_reader = VCF(vcf_path)
    except Exception as e:
        print(f"Error al abrir '{vcf_path}': {e}")
        return

    # Obtenemos la lista de nombres de muestras (en orden)
    sample_names = vcf_reader.samples

    # Construimos e imprimimos la cabecera VCF (columnas fijas + muestras)
    # NOTA: No imprimimos líneas de metadatos (##); solo la línea de columnas
    header_cols = [
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"
    ] + sample_names
    print("\t".join(header_cols))

    # Iterar sobre todas las variantes
    for record in vcf_reader:
        # 1) Campos VCF básicos
        chrom = record.CHROM
        pos   = record.POS
        id_field = record.ID if record.ID else "."
        ref  = record.REF
        alt  = ",".join(record.ALT) if record.ALT else "."
        qual = str(record.QUAL) if record.QUAL is not None else "."

        # FILTER puede venir como lista o cadena
        filt_field = record.FILTER
        if isinstance(filt_field, list):
            filt_str = ";".join(filt_field) if filt_field else "."
        elif isinstance(filt_field, str):
            filt_str = filt_field
        else:
            filt_str = "."

        # 2) INFO: conservar solo SOMATIC (flag). Si no existe, "."
        info_dict = dict(record.INFO)
        info_str = "SOMATIC" if "SOMATIC" in info_dict else "."

        # 3) FORMAT fijo: “GT:DP:BQ”
        fmt_str = "GT:DP:BQ"

        # 4) Por cada muestra: extraer GT, DP y BQ
        genotypes = record.genotypes      # lista de [a1, a2, phased_flag]
        dp_array   = record.format("DP")  # numpy array con DP por muestra (o None)
        bq_array   = record.format("BQ")  # numpy array con BQ por muestra (o None)

        # Construimos la fila de salida
        row = [
            chrom, str(pos), id_field, ref, alt,
            qual, filt_str, info_str, fmt_str
        ]

        for i in range(len(sample_names)):
            # 4.a) Genotipo (GT)
            g = genotypes[i]
            a1, a2 = g[0], g[1]
            if a1 is None or (isinstance(a1, int) and a1 < 0):
                gt_str = "."
            else:
                gt_str = f"{int(a1)}/{int(a2)}"

            # 4.b) DP de muestra
            dp_val = None
            if dp_array is not None:
                try:
                    dp_val = dp_array[i]
                except Exception:
                    dp_val = None
            dp_str = format_scalar_or_array(dp_val)

            # 4.c) BQ de muestra
            bq_val = None
            if bq_array is not None:
                try:
                    bq_val = bq_array[i]
                except Exception:
                    bq_val = None
            bq_str = format_scalar_or_array(bq_val)

            # 4.d) Concatenamos GT, DP y BQ con ":"
            sample_field = f"{gt_str}:{dp_str}:{bq_str}"
            row.append(sample_field)

        # Imprimimos la línea completa de la variante
        print("\t".join(row))

if __name__ == "__main__":
    main()
