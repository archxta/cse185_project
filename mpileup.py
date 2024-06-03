import os
import sys

def analyze_mpileup(mpileup_file, min_var_freq=0.2, min_hom_freq=0.8, p_value=0.01):
    variants = []
    snp_count = 0  # Initialize SNP count
    with open(mpileup_file) as mpileup:
        for line in mpileup:
            fields = line.split()
            chrom, pos, ref_base, count = fields[0], int(fields[1]) - 1, fields[2], int(fields[3])
            alt_base = fields[4] if len(fields) > 4 else '.'
            if alt_base != '.' and alt_base != ref_base:  # Check for valid alternative base (SNP)
                var_freq = count / (count + int(fields[6]))  # Calculate variant frequency
                if var_freq >= min_var_freq:  # Check minimum variant frequency
                    if var_freq >= min_hom_freq or (1 - var_freq) >= min_hom_freq:  # Check minimum frequency for homozygotes
                        if float(fields[5]) <= p_value:  # Check p-value
                            variants.append((chrom, pos, ref_base, alt_base, count))
                            snp_count += 1
    return variants, snp_count

def write_vcf(variants, output_vcf):
    vcf_header = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""
    with open(output_vcf, 'w') as vcf:
        vcf.write(vcf_header + '\n')
        for variant in variants:
            chrom, pos, ref, alt, qual = variant
            vcf.write(f"{chrom}\t{pos+1}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_mpileup> <output_vcf> <min_var_freq>")
        sys.exit(1)

    input_mpileup, output_vcf, min_var_freq = sys.argv[1:]

    # Analyze mpileup file and call variants
    variants, snp_count = analyze_mpileup(input_mpileup, float(min_var_freq))

    # Write VCF file
    write_vcf(variants, output_vcf)
    print(f"Variants written to {output_vcf}")
    print(f"Number of SNPs detected: {snp_count}")

if __name__ == "__main__":
    main()
