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
            ref_base = ref_base.upper()  # Convert reference base to uppercase

            if read_count < min_coverage:
                continue
            
            base_counts = { 'A': 0, 'C': 0, 'G': 0, 'T': 0 }
            i = 0
            while i < len(read_bases):
                base = read_bases[i]
                if base in '.,':  # Reference base match
                    base_counts[ref_base] += 1
                elif base in 'ACGTNacgtn':
                    base_counts[base.upper()] += 1
                elif base == '^':  # Start of a read segment, skip next base (mapping quality)
                    i += 1
                elif base == '$':  # End of a read segment
                    pass
                elif base in '+-':  # Indels, skip the indel length and the indel bases
                    indel_len = ''
                    i += 1
                    while i < len(read_bases) and read_bases[i].isdigit():
                        indel_len += read_bases[i]
                        i += 1
                    i += int(indel_len) - 1
                i += 1
            
            total_reads = sum(base_counts.values())
            for base, count in base_counts.items():
                if base != ref_base:
                    var_freq = count / total_reads
                    if var_freq >= min_var_freq:
                        # Assuming p-value calculation and other filters are applied here
                        # Dummy QUAL score as 100 for simplicity
                        qual_score = 100
                        info = f"DP={total_reads};AF={var_freq:.2f}"
                        variants.append((chrom, pos, ref_base, base, qual_score, info))

    return variants

def write_vcf(variants, output_vcf):
    vcf_header = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""
    with open(output_vcf, 'w') as vcf:
        vcf.write(vcf_header + '\n')
        for variant in variants:
            chrom, pos, ref, alt, qual, info = variant
            vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\tPASS\t{info}\n")

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
