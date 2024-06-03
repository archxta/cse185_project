import os
import sys

def analyze_mpileup(mpileup_file, min_var_freq):
    variants = []
    snp_count = 0

    with open(mpileup_file) as mpileup:
        for line in mpileup:
            fields = line.split()
            
            if len(fields) < 7:
                print(f"Warning: Line ignored: Invalid format for pileup at line: {line}")
                continue
            
            chrom, pos, ref_base, count = fields[0], int(fields[1]), fields[2], int(fields[3])
            alt_base = fields[4] if len(fields) > 4 else '.'

            base_counts = {ref_base: count}
            total_reads = count
            
            if len(fields) > 5:
                alt_bases = fields[5:]
                for alt in alt_bases:
                    base, base_count = alt.split(':')
                    base_count = int(base_count)
                    if base != ref_base:
                        base_counts[base] = base_count
                        total_reads += base_count

            for base, base_count in base_counts.items():
                if base != ref_base and total_reads > 0:
                    var_freq = base_count / total_reads
                    if var_freq >= min_var_freq:
                        variants.append((chrom, pos, ref_base, base, var_freq))
                        snp_count += 1

    return variants, snp_count

def write_vcf(variants, output_vcf):
    vcf_header = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""
    with open(output_vcf, 'w') as vcf:
        vcf.write(vcf_header + '\n')
        for variant in variants:
            chrom, pos, ref, alt, var_freq = variant
            vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 mpileup.py <input_mpileup> <min_var_freq> <output_vcf>")
        sys.exit(1)

    input_mpileup, min_var_freq, output_vcf = sys.argv[1:]

    # Analyze mpileup file and call variants
    variants, snp_count = analyze_mpileup(input_mpileup, float(min_var_freq))

    # Write VCF file
    write_vcf(variants, output_vcf)
    print(f"Variants written to {output_vcf}")
    print(f"Total number of SNPs detected: {snp_count}")

if __name__ == "__main__":
    main()
