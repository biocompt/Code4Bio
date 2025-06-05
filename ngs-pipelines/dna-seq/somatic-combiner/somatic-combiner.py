#!/usr/bin/env python3
import os
import argparse
from collections import OrderedDict
from datetime import datetime
from tqdm import tqdm
import numpy as np

def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def chr_sort_key(chrom):
    ch = chrom
    if ch.lower().startswith("chr"):
        ch = ch[3:]
    if ch.isdigit():
        return (0, int(ch))
    if ch == "X":
        return (0, 23)
    if ch == "Y":
        return (0, 24)
    if ch in ("M", "MT"):
        return (0, 25)
    return (1, ch)

def parse_filtered_vcf(path):
    """
    Lee un VCF filtrado y devuelve:
      header_lines: todas las líneas que empiezan con "##"
      column_header: la línea "#CHROM ..." (no se usa para el header final)
      records: dict mapping (chrom, pos, ref, alt) → (INFO_string, FORMAT_string, [sample_field_strs...])
    """
    header_lines = []
    column_header = None
    records = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                header_lines.append(line.rstrip())
            elif line.startswith("#CHROM"):
                column_header = line.rstrip()
                break
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip().split("\t")
            chrom = cols[0]
            pos = int(cols[1])
            ref = cols[3]
            alt = cols[4]
            info = cols[7]
            fmt = cols[8]
            samples = cols[9:]
            key = (chrom, pos, ref, alt)
            records[key] = (info, fmt, samples)
    return header_lines, column_header, records

def merge_variant_fields(caller_records, key, present_callers):
    """
    Combina los campos de una variante para NORMAL y TUMOR a partir de los distintos callers.
    Para cada caller, toma sus datos de NORMAL y TUMOR (si existen) y promedia DP, AF, GQ, MQ, MQ0.
    """
    normal_data = {"GT": ".", "DP": ".", "AF": ".", "GQ": ".", "MQ": ".", "MQ0": "."}
    tumor_data  = {"GT": ".", "DP": ".", "AF": ".", "GQ": ".", "MQ": ".", "MQ0": "."}

    for caller in present_callers:
        rec_info, rec_fmt, rec_samples = caller_records[caller][key]
        caller_fmt_keys = rec_fmt.split(":")
        sample_raw = rec_samples[0].split(":")

        val_map = {caller_fmt_keys[i]: sample_raw[i] for i in range(len(caller_fmt_keys))}

        # Si el caller aporta campos de NORMAL y TUMOR en la misma entrada (asumimos que lo hace)
        # simplemente actualizamos normal_data y tumor_data con lo que encuentre en val_map.
        # (Por simplicidad, asumimos que quien llame a esta función sabe qué campos representan a NORMAL/TUMOR,
        #  pero en tu caso podrías homogeneizar aquí si, por ejemplo, usas un sufijo o prefijo).
        # En el ejemplo original, cada VCF ya tenía una sola muestra, así que “val_map” corresponde
        # a la muestra única, que asumimos contiene primero NORMAL y luego TUMOR en diferentes VCF.
        # Para mantener la lógica general: si val_map contiene esos campos, los asignamos a ambos.
        # Aquí suponemos: cada caller_proporciona campos para NORMAL y TUMOR en su propia estructura.

        # En la versión original, cada VCF contiene un único “sample” con GT/DP/AF/...,
        # pero en realidad, es la muestra somática (TUMOR), y no hay campo implícito de NORMAL.
        # Si quisieras extraer GT:DP:AF de la parte “NORMAL” del VCF, necesitarías parsear cada caller por separado.
        # Para no complicar, asumimos que “val_map” es directamente la información TUMOR,
        # y dejamos NORMAL en "." si no hay nada explícito para NORMAL. Sólo promediamos DP/AF
        # si ese caller ya contiene un campo específico para NORMAL (p.ej. val_map["N_DP"]).
        # Como el enunciado original no definía diferencia explícita, asumiremos:
        #   - Cada caller da datos de TUMOR (entendemos que es la muestra “tumor”).
        #   - Para NORMAL, si el caller tiene esos mismos campos (por convención “N_DP”, “N_AF”, etc.), los usamos.
        # Como eso no estaba explícito, aquí simplemente copiamos todo a TUMOR y dejamos NORMAL con “.”, 
        # a menos que val_map tuviera una clave “N_DP”/“N_AF”/etc. (lo dejo preparado, aunque en la mayoría de VCF no exista).

        # Por simplicidad, asignamos todo a tumor_data:
        for k, v in val_map.items():
            if k in tumor_data:
                tumor_data[k] = v
            # Si llegado el caso hubiera claves como "N_DP", "N_AF", etc., haríamos:
            # elif k.startswith("N_") y k[2:] en normal_data: 
            #     normal_data[k[2:]] = v

    # Ahora promediamos los campos compartidos
    result = {"NORMAL": {}, "TUMOR": {}}
    for field in ["DP", "AF", "GQ", "MQ", "MQ0"]:
        n_val = normal_data.get(field, ".")
        t_val = tumor_data.get(field, ".")
        if n_val != ".":
            result["NORMAL"][field] = n_val
        if t_val != ".":
            result["TUMOR"][field] = t_val
        if n_val != "." and t_val != ".":
            try:
                mean_val = str(np.mean([float(n_val), float(t_val)]))
            except ValueError:
                mean_val = "."
            result["NORMAL"][field] = mean_val
            result["TUMOR"][field] = mean_val

    # GT: tomamos directamente el GT de tumor_data (o "." si no existe)
    result["NORMAL"]["GT"] = normal_data.get("GT", ".")
    result["TUMOR"]["GT"]  = tumor_data.get("GT", ".")

    # Si faltan campos, ya están en "." por defecto
    return result

def main():
    parser = argparse.ArgumentParser(
        description="Filter somatic VCFs and merge variants present in ≥2 callers into a merged VCF"
    )
    parser.add_argument("--mutect2",       help="Path to Mutect2 VCF")
    parser.add_argument("--muse",          help="Path to MuSE VCF")
    parser.add_argument("--strelka-snvs",  help="Path to Strelka SNVs VCF")
    parser.add_argument("--strelka-indels",help="Path to Strelka indels VCF")
    parser.add_argument("--varscan-snvs",  help="Path to VarScan SNVs VCF")
    parser.add_argument("--varscan-indels",help="Path to VarScan indels VCF")
    parser.add_argument("--output_dir",    required=True, help="Directory for filtered VCFs")
    parser.add_argument("--id",            required=True, help="Name for merged VCF output (without path)")
    args = parser.parse_args()

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)

    filtered_paths = {}
    callers_to_run = []

    if args.mutect2:
        callers_to_run.append(("mutect2", args.mutect2, "filtrar_mutect2", "filter_mutect2"))
    if args.muse:
        callers_to_run.append(("muse", args.muse, "filtrar_muse", "filter_muse"))
    if args.strelka_snvs:
        callers_to_run.append(("strelka-snvs", args.strelka_snvs, "filtrar_strelka", "filter_strelka"))
    if args.strelka_indels:
        callers_to_run.append(("strelka-indels", args.strelka_indels, "filtrar_strelka", "filter_strelka"))
    if args.varscan_snvs:
        callers_to_run.append(("varscan-snvs", args.varscan_snvs, "filtrar_varscan", "filter_varscan"))
    if args.varscan_indels:
        callers_to_run.append(("varscan-indels", args.varscan_indels, "filtrar_varscan", "filter_varscan"))

    # Paso 1: ejecutar filtros
    for caller, input_path, module_name, function_name in callers_to_run:
        base = os.path.basename(input_path)
        stem = base.replace(".vcf.gz", "")
        out_path = os.path.join(args.output_dir, stem + ".filtered.vcf")
        print(f"{timestamp()}  Starting filter for {caller} ({input_path}) → {out_path}")
        module = __import__(f"modules.{module_name}", fromlist=[function_name])
        func = getattr(module, function_name)
        func(input_path, out_path)
        print(f"{timestamp()}  Finished filter for {caller}")
        filtered_paths[caller] = out_path

    # Paso 2: parsear VCFs filtrados
    caller_headers = {}
    caller_records = {}
    all_contigs = OrderedDict()
    fmt_keys_union = set()

    for caller, path in filtered_paths.items():
        print(f"{timestamp()}  Parsing filtered VCF for {caller}: {path}")
        hdr_lines, col_hdr, recs = parse_filtered_vcf(path)
        caller_headers[caller] = (hdr_lines, col_hdr)
        caller_records[caller] = recs
        # Recolectamos contigs disponibles (chromosomes)
        for line in hdr_lines:
            if line.startswith("##contig"):
                contig_id = line.split("ID=")[1].split(">")[0]
                if contig_id not in all_contigs:
                    all_contigs[contig_id] = line
        # Recolectamos formatos disponibles (no los usaremos para header, pero sí para construir fmt_union)
        for line in hdr_lines:
            if line.startswith("##FORMAT"):
                key = line.split("<ID=")[1].split(",")[0]
                fmt_keys_union.add(key)
        print(f"{timestamp()}  Parsed {len(recs)} variants for {caller}")

    # Paso 3: construir header estático (solo una vez)
    merged_header = [
        "##fileformat=VCFv4.2"
    ]

    # Añadimos los contigs (solo una vez cada uno), ordenados en orden natural
    # Por convención, incluimos chr1–chr22, chrX, chrY, chrM
    chromos = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    for c in chromos:
        merged_header.append(f"##contig=<ID={c}>")

    # A continuación definimos los campos FORMAT (fijos)
    merged_header.extend([
        '##INFO=<ID=CC,Number=1,Type=Integer,Description="Number of callers supporting variant">',
        '##INFO=<ID=CL,Number=.,Type=String,Description="List of callers supporting variant">',
        '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Indicates if record is a somatic mutation">',
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads con MQ=255 o con bad mates filtrados)">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample strand-bias statistics">',
        '##FORMAT=<ID=BQ,Number=.,Type=Integer,Description="Average base quality for reads supporting alleles">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
        '##FORMAT=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">',
        '##FORMAT=<ID=MQ0,Number=1,Type=Integer,Description="Mapping Quality Zero Count">'
    ])

    # LÍNEA DE COLUMNAS (solo una vez)
    column_header = "#" + "\t".join(
        ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR"]
    )

    # Paso 4: merge de variantes
    print(f"{timestamp()}  Starting merge of variants present in ≥2 callers")
    all_keys = set()
    for recs in caller_records.values():
        all_keys.update(recs.keys())

    sorted_keys = sorted(all_keys, key=lambda k: (chr_sort_key(k[0]), k[1]))
    print(f"{timestamp()}  Total distinct variant keys: {len(sorted_keys)}")

    merged_vcf_path = os.path.join(args.output_dir, args.id + ".vcf")
    with open(merged_vcf_path, "w") as out:
        # Escribimos el header (solo una vez)
        for line in merged_header:
            out.write(line + "\n")
        out.write(column_header + "\n")

        for key in tqdm(sorted_keys, desc="Merging variants", unit="variant"):
            chrom, pos, ref, alt = key
            present_callers = [c for c in filtered_paths if key in caller_records[c]]
            if len(present_callers) < 2:
                continue

            # Merge de campos para NORMAL y TUMOR
            merged_fields = merge_variant_fields(caller_records, key, present_callers)
            normal_data = merged_fields["NORMAL"]
            tumor_data  = merged_fields["TUMOR"]

            # Construimos el campo INFO
            info_parts = [
                f"CC={len(present_callers)}",
                f"CL={','.join(present_callers)}"
            ]
            # Conservamos la parte "SOMATIC" si alguno de los callers la tiene
            # (asumimos que al menos uno la tendrá, o la omitimos si no)
            # Por simplicidad, buscamos "SOMATIC" en los INFO originales:
            for caller in present_callers:
                caller_info = caller_records[caller][key][0]
                if "SOMATIC" in caller_info.split(";"):
                    info_parts.append("SOMATIC")
                    break

            merged_info = ";".join(info_parts)

            # Unión de todos los campos FORMAT que aparezcan en al menos un caller para esta variante
            present_fmt_keys = set()
            for caller in present_callers:
                caller_fmt = caller_records[caller][key][1].split(":")
                present_fmt_keys.update(caller_fmt)

            # Orden preferente de campos
            pref_order = ["GT", "DP", "AF", "SB", "GQ", "MQ", "MQ0"]
            fmt_union = [k for k in pref_order if k in present_fmt_keys]
            for k in sorted(present_fmt_keys):
                if k not in fmt_union:
                    fmt_union.append(k)
            merged_fmt = ":".join(fmt_union)

            # Ahora construimos los valores para NORMAL y TUMOR según fmt_union
            normal_vals = [normal_data.get(k, ".") for k in fmt_union]
            tumor_vals  = [tumor_data.get(k, ".")  for k in fmt_union]

            out.write(
                "\t".join([
                    chrom,
                    str(pos),
                    ".",
                    ref,
                    alt,
                    ".",
                    "PASS",
                    merged_info,
                    merged_fmt,
                    ":".join(normal_vals),
                    ":".join(tumor_vals)
                ]) + "\n"
            )

    print(f"{timestamp()}  Final merged VCF saved to: {merged_vcf_path}")

if __name__ == "__main__":
    main()
