import pysam
from pysam import FastaFile

def analyze_reads(bam_file, reference_fasta):
    bam = pysam.AlignmentFile(bam_file, "rb")
    reference = FastaFile(reference_fasta)
    variants = []

    for pileupcolumn in bam.pileup():
        ref_pos = pileupcolumn.reference_pos
        ref_name = bam.get_reference_name(pileupcolumn.tid)
        ref_base = reference.fetch(ref_name, ref_pos, ref_pos+1)
        base_counts = {}
        total_reads = 0

        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                base_counts[base] = base_counts.get(base, 0) + 1
                total_reads += 1

        for base, count in base_counts.items():
            if base != ref_base and count / total_reads > 0.2:
                variants.append((ref_name, ref_pos, ref_base, base, count / total_reads))

    return variants
