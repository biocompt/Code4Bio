#!/usr/bin/env python3
"""
filtrar_mutect2_vcf.py

Este script carga un VCF de Mutect2 y, sin cambiar el esquema de columnas
de VCF, reescribe únicamente los subcampos que nos interesan:

  • En INFO:  conservar solo “DP”
  • En FORMAT: conservar “GT”, “DP” (tomado siempre de INFO:DP), “AF” y “SB”

Imprime por pantalla todas las variantes en un VCF “reescrito”,
con las mismas columnas VCF (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, muestras…),
pero en INFO únicamente “DP=<valor>” y en FORMAT “GT:DP:AF:SB” para cada muestra,
donde DP viene de INFO, y AF y SB se toman del FORMAT original.

Uso:
    python filtrar_mutect2_vcf.py -v ruta/al/Mutect2.vcf.gz
"""

import argparse
import numpy as np
from cyvcf2 import VCF

def format_scalar_or_array(val):
    """
    Dado un valor que puede ser escalar, array o numpy scalar,
    devuelve una cadena sin corchetes y con comas si hay varios elementos.
    Si es None o contiene NaN, devuelve “.”.
    """
    if val is None:
        return "."
    # numpy scalar (int64 o float64)
    if isinstance(val, (np.integer, np.floating)):
        # Para float, lo convertimos a str tal cual; para entero, a int->str
        if isinstance(val, np.integer):
            return str(int(val))
        else:
            return str(val)
    # numpy.ndarray o lista/tupla
    if isinstance(val, (list, tuple, np.ndarray)):
        elems = []
        for x in val:
            if x is None:
                elems.append(".")
            elif isinstance(x, (np.integer, np.floating)):
                if isinstance(x, np.integer):
                    elems.append(str(int(x)))
                else:
                    elems.append(str(x))
            else:
                elems.append(str(x))
        return ",".join(elems)
    # Caso escalar puro (int, float, str)
    return str(val)

def main():
    parser = argparse.ArgumentParser(
        description="Reescribir un VCF de Mutect2: INFO→solo DP, FORMAT→GT,DP,AF,SB"
    )
    parser.add_argument(
        "-v", "--vcf",
        required=True,
        help="Ruta al VCF de Mutect2 que se va a procesar"
    )
    args = parser.parse_args()
    vcf_path = args.vcf

    # Intentamos abrir el VCF
    try:
        vcf_reader = VCF(vcf_path)
    except Exception as e:
        print(f"Error al abrir '{vcf_path}': {e}")
        return

    # Lista de muestras (en orden)
    sample_names = vcf_reader.samples

    # Construir e imprimir la línea de columnas VCF sin los metacampos (##)
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

        # FILTER puede venir como lista o cadena
        filt_field = record.FILTER
        if isinstance(filt_field, list):
            filt_str = ";".join(filt_field) if filt_field else "."
        elif isinstance(filt_field, str):
            filt_str = filt_field
        else:
            filt_str = "."

        # 2) INFO: conservar solo DP
        info_dict = dict(record.INFO)
        dp_info = info_dict.get("DP")
        if dp_info is not None:
            info_str = f"DP={dp_info}"
        else:
            info_str = "DP=."

        # 3) FORMAT fijo: “GT:DP:AF:SB”
        fmt_str = "GT:DP:AF:SB"

        # 4) Recuperar datos por muestra
        genotypes = record.genotypes      # lista de [a1, a2, phased_flag]
        # DP por muestra no lo usamos; tomamos dp_info siempre
        af_array = record.format("AF")    # numpy array (n_muestras, ...) o None
        sb_array = record.format("SB")    # numpy array (n_muestras, ...) o None

        # Construir la fila de salida
        row = [
            chrom, str(pos), id_field, ref, alt,
            qual, filt_str, info_str, fmt_str
        ]

        for i in range(len(sample_names)):
            # 4.a) GT
            g = genotypes[i]
            a1, a2 = g[0], g[1]
            if a1 is None or (isinstance(a1, int) and a1 < 0):
                gt_str = "."
            else:
                gt_str = f"{int(a1)}/{int(a2)}"

            # 4.b) DP: siempre dp_info
            dp_str = str(dp_info) if dp_info is not None else "."

            # 4.c) AF si existe en FORMAT
            af_val = None
            if af_array is not None:
                try:
                    # Puede ser array 1-D o 2-D, pero cyvcf2 suele devolver (n_muestras,) para AF
                    af_val = af_array[i]
                except Exception:
                    af_val = None
            af_str = format_scalar_or_array(af_val)

            # 4.d) SB si existe en FORMAT
            sb_val = None
            if sb_array is not None:
                try:
                    sb_val = sb_array[i]
                except Exception:
                    sb_val = None
            sb_str = format_scalar_or_array(sb_val)

            # 4.e) Concatenar en “GT:DP:AF:SB”
            sample_field = f"{gt_str}:{dp_str}:{af_str}:{sb_str}"
            row.append(sample_field)

        # Imprimir la línea de la variante
        print("\t".join(row))

if __name__ == "__main__":
    main()
