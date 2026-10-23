import pysam
import numpy as np
import matplotlib.pyplot as plt

# Abrir BAM
bam_file = "modulo6.bam"
bam = pysam.AlignmentFile(bam_file, "rb")

prof_dict= {}

for lectura in bam.fetch():
    if lectura.is_unmapped:
        continue
    for pos in range (lectura.reference_start, lectura.reference_end):
        prof_dict[pos]= prof_dict.get(pos, 0) + 1



Num_pos_cubiertas= len(prof_dict)
primer_pos_leida = min(prof_dict.keys())
ultima_pos_leida = max(prof_dict.keys())

print(f"Cantidad de Posiciones Cubiertas: {Num_pos_cubiertas}")
print(f"Desde: {primer_pos_leida} a {ultima_pos_leida}")

positions = sorted(prof_dict.keys())
depths = np.array([prof_dict[p] for p in positions])
plt.figure(figsize=(12, 4))
plt.plot(positions, depths)
plt.xlabel("Posicion")
plt.ylabel("Profundidad")
plt.title(f"Profundidad para c/ Posicion en {bam_file}")
plt.tight_layout()
plt.show()

mean_depth = depths.mean()
sd_depth = depths.std()
print("Profundidad Promedio:", mean_depth)
print("DS (Profundidad):", sd_depth)

read_lengths = []

for read in bam.fetch():
    if read.is_unmapped:
        continue
    read_lengths.append(read.query_length)

mean_len = np.mean(read_lengths)
sd_len   = np.std(read_lengths)
print(f"Largo Promedio por lectura: {mean_len}")
print(f"DS(largo de lectura): {sd_len}")

thresholds = [1,5,10,20,30,50,100]

for value in thresholds:
    covered = np.sum(depths >= value)
    total   = len(depths)
    coverage_percent = covered / total * 100
    print(f"Cobertura para umbral {value}× : {coverage_percent:.2f}%")


import pysam
from collections import defaultdict

bam_file = "modulo6.bam"
ref_fasta = "Homo_sapiens_assembly38.chr1.fasta"

bam = pysam.AlignmentFile(bam_file, "rb")
fasta = pysam.FastaFile(ref_fasta)

# Dictionary:
# variants[pos] = { "REF": "A", "TOTAL_DEPTH": X, "ALTS": { "G": 5, "T": 1, "INS": 2, "DEL": 3 } }
variants = {}

for pileupcolumn in bam.pileup(stepper='samtools', truncate=True):

    pos = pileupcolumn.reference_pos   # 0-based reference coordinate
    ref_base = fasta.fetch(bam.get_reference_name(0), pos, pos+1).upper()

    total_depth = pileupcolumn.nsegments
    alt_counts = defaultdict(int)

    # Inspect each read covering this position
    for pileupread in pileupcolumn.pileups:

        if pileupread.is_del:
            # deletion relative to reference
            alt_counts["DEL"] += 1
            continue

        if pileupread.is_refskip:
            continue

        base = pileupread.alignment.query_sequence[pileupread.query_position].upper()

        if base != ref_base:
            alt_counts[base] += 1

    # Record only positions where any difference exists
    if len(alt_counts) > 0:
        variants[pos] = {
            "REF": ref_base,
            "TOTAL_DEPTH": total_depth,
            "ALTS": dict(alt_counts)
        }



# primeras posiciones
for p in list(variants.keys())[:10]:
    print(p, variants[p])

# ------------------------
# Thresholds (fractions + absolute minimum)
# ------------------------
DP_MIN = 10           # minimal total depth
AD_ABS_MIN = 3        # alt must be supported by at least this many reads
VAF_HET_MIN = 0.20    # lower bound for heterozygous VAF
VAF_HET_MAX = 0.80    # upper bound for heterozygous VAF
VAF_HOM_MIN = 0.80    # lower bound for homozygous alt VAF

output_vcf = "modulo6_mycalls_adaptive.vcf"

with open(output_vcf, "w") as vcf:
    vcf.write("##fileformat=VCFv4.2\n")
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    chrom = bam.get_reference_name(0)

    for pos, data in variants.items():
        ref = data["REF"]
        dp  = data["TOTAL_DEPTH"]

        if dp < DP_MIN:
            continue

        for alt, ad in data["ALTS"].items():
            # absolute alt count guard
            if ad < AD_ABS_MIN:
                continue

            vaf = ad / dp

            # classify by VAF
            if vaf >= VAF_HOM_MIN:
                gt = "1/1"
            elif VAF_HET_MIN <= vaf < VAF_HET_MAX:
                gt = "0/1"
            else:
                # VAF too low or ambiguous (e.g., 0.85 < vaf < 0.8 won't happen with these numbers,
                # but catches low-mid values)
                continue

            info = f"DP={dp};AD={ad};AF={vaf:.3f};GT={gt}"
            vcf.write(f"{chrom}\t{pos+1}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n")

print("VCF generated:", output_vcf)