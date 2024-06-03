import os
import sys

def analyze_mpileup(mpileup_file):
    variants = []
    with open(mpileup_file) as mpileup:
        for line in mpileup:
            fields = line.split()
            chrom, pos, ref_base, count = fields[0], int(fields[1]) - 1, fields[2], int(fields[3])
            alt_base = fields[4] if len(fields) > 4 else '.'
            variants.append((chrom, pos, ref_base, alt_base, count))
    return variants

def write_vcf(variants, output_vcf):
    vcf_header = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""
    with open(output_vcf, 'w') as vcf:
        vcf.write(vcf_header + '\n')
        for variant in variants:
            chrom, pos, ref, alt, qual = variant
            vcf.write(f"{chrom}\t{pos+1}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_mpileup> <output_vcf>")
        sys.exit(1)

    input_mpileup, output_vcf = sys.argv[1:]

    # Analyze mpileup file and call variants
    variants = analyze_mpileup(input_mpileup)

    # Write VCF file
    write_vcf(variants, output_vcf)
    print(f"Variants written to {output_vcf}")

if __name__ == "__main__":
    main()
