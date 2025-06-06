#!/usr/bin/env python3
import os
import argparse
from collections import OrderedDict, defaultdict
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

def simplify_caller_name(caller):
    """Simplifica nombres de callers quitando todo después del guion"""
    return caller.split("-")[0]

def parse_filtered_vcf(path):
    header_lines = []
    column_header = None
    records = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                header_lines.append(line.rstrip())
            elif line.startswith("#CHROM"):
                column_header = line.rstrip()
                samples = line.rstrip().split("\t")[9:]
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
            sample_data = cols[9:]
            
            fmt_keys = fmt.split(":")
            normal_dict = dict(zip(fmt_keys, sample_data[0].split(":")))
            tumor_dict = dict(zip(fmt_keys, sample_data[1].split(":"))) if len(sample_data) > 1 else {}
            
            key = (chrom, pos, ref, alt)
            records[key] = (info, fmt, normal_dict, tumor_dict)
    return header_lines, column_header, records

def get_consensus_gt(gt_list):
    if not gt_list:
        return "."
    valid_gts = [gt for gt in gt_list if gt not in (".", "./.", ".|.")]
    if not valid_gts:
        return "."
    if all(gt == valid_gts[0] for gt in valid_gts):
        return valid_gts[0]
    return "."

def merge_variant_fields(caller_records, key, present_callers):
    normal_data = defaultdict(list)
    tumor_data = defaultdict(list)
    
    numeric_fields = ["DP", "AF", "GQ", "MQ", "MQ0"]
    all_fields = ["GT"] + numeric_fields
    
    for caller in present_callers:
        info, fmt, normal_dict, tumor_dict = caller_records[caller][key]
        
        for field in all_fields:
            if field in normal_dict and normal_dict[field] != ".":
                if field == "GT":
                    normal_data["GT"].append(normal_dict[field])
                else:
                    try:
                        val = float(normal_dict[field])
                        normal_data[field].append(val)
                    except ValueError:
                        pass
        
        for field in all_fields:
            if field in tumor_dict and tumor_dict[field] != ".":
                if field == "GT":
                    tumor_data["GT"].append(tumor_dict[field])
                else:
                    try:
                        val = float(tumor_dict[field])
                        tumor_data[field].append(val)
                    except ValueError:
                        pass
    
    def process_sample_data(data_dict):
        result = {}
        result["GT"] = get_consensus_gt(data_dict.get("GT", []))
        
        for field in numeric_fields:
            values = data_dict.get(field, [])
            if not values:
                result[field] = "."
            else:
                if field == "AF":
                    avg = round(np.mean(values), 2)
                    result[field] = f"{avg:.2f}"
                else:
                    avg = int(round(np.mean(values)))
                    result[field] = str(avg)
        return result
    
    return {
        "NORMAL": process_sample_data(normal_data),
        "TUMOR": process_sample_data(tumor_data)
    }

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
        callers_to_run.append(("mutect2", args.mutect2, "filter_mutect2", "filter_mutect2"))
    if args.muse:
        callers_to_run.append(("muse", args.muse, "filter_muse", "filter_muse"))
    if args.strelka_snvs:
        callers_to_run.append(("strelka-snvs", args.strelka_snvs, "filter_strelka", "filter_strelka"))
    if args.strelka_indels:
        callers_to_run.append(("strelka-indels", args.strelka_indels, "filter_strelka", "filter_strelka"))
    if args.varscan_snvs:
        callers_to_run.append(("varscan-snvs", args.varscan_snvs, "filter_varscan", "filter_varscan"))
    if args.varscan_indels:
        callers_to_run.append(("varscan-indels", args.varscan_indels, "filter_varscan", "filter_varscan"))

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

    caller_headers = {}
    caller_records = {}
    all_contigs = OrderedDict()
    fmt_keys_union = set()

    for caller, path in filtered_paths.items():
        print(f"{timestamp()}  Parsing filtered VCF for {caller}: {path}")
        hdr_lines, col_hdr, recs = parse_filtered_vcf(path)
        caller_headers[caller] = (hdr_lines, col_hdr)
        caller_records[caller] = recs
        for line in hdr_lines:
            if line.startswith("##contig"):
                contig_id = line.split("ID=")[1].split(">")[0]
                if contig_id not in all_contigs:
                    all_contigs[contig_id] = line
        for line in hdr_lines:
            if line.startswith("##FORMAT"):
                key = line.split("<ID=")[1].split(",")[0]
                fmt_keys_union.add(key)
        print(f"{timestamp()}  Parsed {len(recs)} variants for {caller}")

    merged_header = [
        "##fileformat=VCFv4.2"
    ]
    chromos = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    for c in chromos:
        merged_header.append(f"##contig=<ID={c}>")

    merged_header.extend([
        '##INFO=<ID=CC,Number=1,Type=Integer,Description="Number of callers supporting variant">',
        '##INFO=<ID=CL,Number=.,Type=String,Description="List of callers supporting variant">',
        '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Indicates if record is a somatic mutation">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FORMAT=<ID=MQ,Number=1,Type=Integer,Description="Mapping quality">',
        '##FORMAT=<ID=MQ0,Number=1,Type=Integer,Description="Count of MAPQ=0 reads">'
    ])

    column_header = "#" + "\t".join(
        ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR"]
    )

    print(f"{timestamp()}  Starting merge of variants present in ≥2 callers")
    all_keys = set()
    for recs in caller_records.values():
        all_keys.update(recs.keys())

    sorted_keys = sorted(all_keys, key=lambda k: (chr_sort_key(k[0]), k[1]))
    print(f"{timestamp()}  Total distinct variant keys: {len(sorted_keys)}")

    merged_vcf_path = os.path.join(args.output_dir, args.id + ".vcf")
    with open(merged_vcf_path, "w") as out:
        for line in merged_header:
            out.write(line + "\n")
        out.write(column_header + "\n")

        for key in tqdm(sorted_keys, desc="Merging variants", unit="variant"):
            chrom, pos, ref, alt = key
            present_callers = [c for c in filtered_paths if key in caller_records[c]]
            if len(present_callers) < 2:
                continue

            merged_fields = merge_variant_fields(caller_records, key, present_callers)
            normal_data = merged_fields["NORMAL"]
            tumor_data = merged_fields["TUMOR"]

            simplified_callers = [simplify_caller_name(c) for c in present_callers]
            unique_callers = sorted(list(set(simplified_callers)))  # Eliminar duplicados
            
            info_parts = [
                f"CC={len(unique_callers)}",
                f"CL={','.join(unique_callers)}"
            ]
            for caller in present_callers:
                caller_info = caller_records[caller][key][0]
                if "SOMATIC" in caller_info.split(";"):
                    info_parts.append("SOMATIC")
                    break

            merged_info = ";".join(info_parts)

            fmt_union = ["GT", "DP", "AF", "GQ", "MQ", "MQ0"]
            
            normal_vals = [normal_data.get(f, ".") for f in fmt_union]
            tumor_vals = [tumor_data.get(f, ".") for f in fmt_union]

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
                    ":".join(fmt_union),
                    ":".join(normal_vals),
                    ":".join(tumor_vals)
                ]) + "\n"
            )

    print(f"{timestamp()}  Final merged VCF saved to: {merged_vcf_path}")

if __name__ == "__main__":
    main()
