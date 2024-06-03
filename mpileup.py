import sys

def analyze_mpileup(mpileup_file, min_var_freq=0.2, min_coverage=8, p_value=0.01):
    variants = []
    with open(mpileup_file) as mpileup:
        for line in mpileup:
            fields = line.strip().split('\t')
            if len(fields) < 6:
                print(f"Warning: Line ignored: Invalid format for pileup at line {line.strip()}")
                continue
            
            chrom, pos, ref_base, read_count, read_bases, base_qualities = fields[:6]
            pos = int(pos)
            read_count = int(read_count)

            if read_count < min_coverage:
                continue
            
            base_counts = { 'A': 0, 'C': 0, 'G': 0, 'T': 0 }
            for base in read_bases:
                if base in base_counts:
                    base_counts[base] += 1
            
            total_reads = sum(base_counts.values())
            for base, count in base_counts.items():
                if base != ref_base:
                    var_freq = count / total_reads
                    if var_freq >= min_var_freq:
                        # Assuming p-value calculation and other filters are applied here
                        variants.append((chrom, pos, ref_base, base, count, var_freq))

    return variants

def write_vcf(variants, output_vcf):
    vcf_header = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""
    with open(output_vcf, 'w') as vcf:
        vcf.write(vcf_header + '\n')
        for variant in variants:
            chrom, pos, ref, alt, count, var_freq = variant
            vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_mpileup> <output_vcf> <min_var_freq>")
        sys.exit(1)

    input_mpileup, output_vcf, min_var_freq = sys.argv[1:]

    # Analyze mpileup file and call variants
    variants = analyze_mpileup(input_mpileup, float(min_var_freq))

    # Write VCF file
    write_vcf(variants, output_vcf)
    print(f"Variants written to {output_vcf}")
    print(f"Number of SNPs detected: {len(variants)}")

if __name__ == "__main__":
    main()
