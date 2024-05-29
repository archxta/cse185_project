def write_vcf(variants, output_file):
    with open(output_file, 'w') as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for variant in variants:
            chrom, pos, base, confidence = variant
            vcf.write(f"{chrom}\t{pos + 1}\t.\t.\t{base}\t{confidence}\tPASS\t.\n")
