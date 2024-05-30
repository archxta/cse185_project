import pysam
import os
import sys
import pybwa
from pysam import FastaFile

def index_reference_genome(reference_fasta):
    bwa_index = pybwa.Index(reference_fasta)
    bwa_index.build_index()

def align_reads_to_reference(bwa_index, reads_fastq):
    alignment = bwa_index.align(reads_fastq)
    return alignment

def index_bam_file(bam_file):
    os.system(f'samtools index {bam_file}')
    print(f"Indexed {bam_file}")

def analyze_reads(bam_file, reference_fasta):
    bam = pysam.AlignmentFile(bam_file, "rb")
    reference = FastaFile(reference_fasta)
    variants = []

    for pileupcolumn in bam.pileup():
        ref_pos = pileupcolumn.reference_pos
        ref_base = reference.fetch(reference.get_reference_name(pileupcolumn.tid), ref_pos, ref_pos+1)
        base_counts = {}
        total_reads = 0

        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                base_counts[base] = base_counts.get(base, 0) + 1
                total_reads += 1

        for base, count in base_counts.items():
            if base != ref_base and count / total_reads > 0.2:
                variants.append((reference.get_reference_name(pileupcolumn.tid), ref_pos, ref_base, base, count / total_reads))

    return variants

def write_vcf(variants, reference_fasta, output_vcf):
    vcf_header = f"""##fileformat=VCFv4.2
##reference={reference_fasta}
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""
    with open(output_vcf, 'w') as vcf:
        vcf.write(vcf_header + '\n')
        for variant in variants:
            chrom, pos, ref, alt, qual = variant
            vcf.write(f"{chrom}\t{pos+1}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")

def main():
    if len(sys.argv) != 5:
        print("Usage: python scripts/code.py <reference_fasta> <reads_fastq> <output_bam> <output_vcf>")
        sys.exit(1)

    reference_fasta, reads_fastq, output_bam, output_vcf = sys.argv[1:]

    index_reference_genome(reference_fasta)
    bwa_index = pybwa.Index(reference_fasta)
    alignment = align_reads_to_reference(bwa_index, reads_fastq)

    with open(output_bam, 'wb') as bam_out:
        bam_out.write(alignment)

    index_bam_file(output_bam)

    variants = analyze_reads(output_bam, reference_fasta)
    write_vcf(variants, reference_fasta, output_vcf)

if __name__ == "__main__":
    main()
