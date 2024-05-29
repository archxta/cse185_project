import pysam

def index_bam_file(bam_file):
    os.system(f'samtools index {bam_file}')
    print(f"Indexed {bam_file}")

def analyze_reads(bam, reference_name):
    variants = []

    for pileupcolumn in bam.pileup(reference_name):
        ref_pos = pileupcolumn.reference_pos
        base_counts = {}
        total_reads = 0

        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                if base in base_counts:
                    base_counts[base] += 1
                else:
                    base_counts[base] = 1
                total_reads += 1

        for base, count in base_counts.items():
            if count / total_reads > 0.2: 
                variants.append((reference_name, ref_pos, base, count / total_reads))

    return variants

