#!/usr/bin/env python3
"""
filtrar_varscan_vcf.py

Este script carga un VCF de VarScan (SNVs o indels) y
reescribe el VCF manteniendo únicamente:

  • En INFO: conservar “DP” y “SOMATIC” (flag), si existen.
  • En FORMAT: conservar “GT”, “GQ”, “DP” y “FREQ” (renombrado como AF),
    extrayendo estos valores directamente de cada muestra, sin corchetes,
    tratando correctamente “.” y redondeando porcentajes a 3 decimales
    como valores decimales (por ejemplo, “35.71%” → “0.357”).

Reglas de extracción:
  - GT: a partir de record.genotypes.
  - GQ: a partir de record.format("GQ"); si es un entero negativo (sentinela), se considera “.”.
  - DP: a partir de record.format("DP"); si es entero negativo, se considera “.”.
  - FREQ (AF): a partir de record.format("FREQ”), que devuelve bytes o numpy.bytes_,
    se decodifica la cadena; si acaba en ‘%’, se convierte a decimal y se redondea a 3 decimales.
    Si es “.” se deja “.”.

Imprime todas las variantes en un VCF reescrito por pantalla,
con las mismas columnas originales, pero:
  - INFO solo “DP=<valor>” y “SOMATIC” (si existe).
  - FORMAT “GT:GQ:DP:AF” por muestra.

Uso:
    python filtrar_varscan_vcf.py -v ruta/al/VarScan.vcf.gz
"""

import argparse
import numpy as np
from cyvcf2 import VCF

def get_format_value(array, index):
    """
    Dado un array devuelto por record.format("FIELD") y un índice de muestra,
    devuelve un valor procesado:
      - Si array es None, retorna None.
      - raw = array[index]
      - Si raw es np.ndarray, lista o tupla, extraer primer elemento.
      - Si raw es np.integer y raw < 0, retorna None (sentinela).
      - Si raw es np.integer >= 0, retorna int(raw).
      - Si raw es np.floating y np.isnan(raw), retorna None.
      - Si raw es np.floating, retorna float(raw).
      - Si raw es int nativo:
           • si raw < 0, retorna None; sino, retorna int.
      - Si raw es float nativo y np.isnan(raw), retorna None; sino, retorna float.
      - Si raw es bytes o np.bytes_, decodificar a str:
           • si str == ".", retorna None
           • si str acaba en "%", convertir porcentaje a decimal, redondear a 3 decimales y retornar float
           • sino, retornar str.
      - Si raw es str:
           • si str == ".", retorna None
           • si str acaba en "%", convertir porcentaje a decimal, redondear a 3 decimales y retornar float
           • sino, retornar str.
    """
    if array is None:
        return None
    try:
        raw = array[index]
    except Exception:
        return None

    # Si es array/ndarray/lista/tupla, extraer primer elemento
    if isinstance(raw, (list, tuple, np.ndarray)):
        if len(raw) == 0:
            return None
        raw = raw[0]

    # NumPy integer
    if isinstance(raw, np.integer):
        val = int(raw)
        if val < 0:
            return None
        return val

    # NumPy floating
    if isinstance(raw, np.floating):
        if np.isnan(raw):
            return None
        return float(raw)

    # Int nativo
    if isinstance(raw, int):
        if raw < 0:
            return None
        return raw

    # Float nativo
    if isinstance(raw, float):
        if np.isnan(raw):
            return None
        return raw

    # Bytes o np.bytes_
    if isinstance(raw, (bytes, np.bytes_)):
        try:
            s = raw.decode("utf-8")
        except Exception:
            s = str(raw)
        s = s.strip()
        if s == ".":
            return None
        if s.endswith("%"):
            try:
                num = float(s.rstrip("%")) / 100.0
                return round(num, 3)
            except Exception:
                return None
        return s

    # Str
    if isinstance(raw, str):
        s = raw.strip()
        if s == ".":
            return None
        if s.endswith("%"):
            try:
                num = float(s.rstrip("%")) / 100.0
                return round(num, 3)
            except Exception:
                return None
        return s

    # Cualquier otro caso
    return str(raw)

def format_field(val):
    """
    Toma el valor devuelto por get_format_value y lo convierte a string:
      - Si val es None, retorna "."
      - Si es int, retorna str(val)
      - Si es float, retorna str(val) (mostrando decimales)
      - Si es str, retorna val
      - Cualquier otro caso, convertir a str
    """
    if val is None:
        return "."
    if isinstance(val, int):
        return str(val)
    if isinstance(val, float):
        return str(val)
    if isinstance(val, str):
        return val
    return str(val)

def main():
    parser = argparse.ArgumentParser(
        description="Reescribir un VCF de VarScan: INFO→solo DP,SOMATIC; FORMAT→GT:GQ:DP:AF"
    )
    parser.add_argument(
        "-v", "--vcf",
        required=True,
        help="Ruta al VCF de VarScan (SNVs o indels) que se va a procesar"
    )
    args = parser.parse_args()
    vcf_path = args.vcf

    try:
        vcf_reader = VCF(vcf_path)
    except Exception as e:
        print(f"Error al abrir '{vcf_path}': {e}")
        return

    sample_names = vcf_reader.samples  # lista de muestras

    # Imprimir línea de columnas VCF (sin metacampos ##)
    header_cols = [
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"
    ] + sample_names
    print("\t".join(header_cols))

    for record in vcf_reader:
        # 1) Campos básicos
        chrom = record.CHROM
        pos   = record.POS
        id_field = record.ID if record.ID else "."
        ref  = record.REF
        alt  = ",".join(record.ALT) if record.ALT else "."
        qual = str(record.QUAL) if record.QUAL is not None else "."

        # 2) FILTER: lista o cadena
        filt = record.FILTER
        if isinstance(filt, list):
            filt_str = ";".join(filt) if filt else "."
        elif isinstance(filt, str):
            filt_str = filt
        else:
            filt_str = "."

        # 3) INFO: conservar DP y SOMATIC
        info_dict = dict(record.INFO)
        info_parts = []
        if "DP" in info_dict:
            dp_info = info_dict.get("DP")
            info_parts.append(f"DP={dp_info}")
        if "SOMATIC" in info_dict:
            info_parts.append("SOMATIC")
        info_str = ";".join(info_parts) if info_parts else "."

        # 4) FORMAT fijo
        fmt_str = "GT:GQ:DP:AF"

        # 5) Extraer datos por muestra
        genotypes  = record.genotypes       # lista de [a1, a2, phased]
        gq_array    = record.format("GQ")   # array o None
        dp_array    = record.format("DP")   # array o None
        freq_array  = record.format("FREQ") # array o None

        # 6) Construir fila base
        row = [
            chrom, str(pos), id_field, ref, alt,
            qual, filt_str, info_str, fmt_str
        ]

        # 7) Por cada muestra
        if sample_names:
            for i in range(len(sample_names)):
                # 7.a) GT
                g = genotypes[i]
                a1, a2 = g[0], g[1]
                if a1 is None or (isinstance(a1, int) and a1 < 0):
                    gt_str = "."
                else:
                    gt_str = f"{int(a1)}/{int(a2)}"

                # 7.b) GQ
                raw_gq = get_format_value(gq_array, i)
                gq_str = format_field(raw_gq)

                # 7.c) DP
                raw_dp = get_format_value(dp_array, i)
                dp_str = format_field(raw_dp)

                # 7.d) FREQ → AF (con porcentaje a decimal y redondeado)
                raw_af = get_format_value(freq_array, i)
                af_str = format_field(raw_af)

                row.append(f"{gt_str}:{gq_str}:{dp_str}:{af_str}")
        else:
            # Si no hay muestras, agregamos un solo "."
            row.append(".")

        # 8) Imprimir línea completa
        print("\t".join(row))

if __name__ == "__main__":
    main()
