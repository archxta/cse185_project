import os

def write_vcf(variants, reference_fasta, output_vcf):
    vcf_header = f"""##fileformat=VCFv4.2
##reference={reference_fasta}
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""
    os.makedirs(os.path.dirname(output_vcf), exist_ok=True)
    with open(output_vcf, 'w') as vcf:
        vcf.write(vcf_header + '\n')
        for variant in variants:
            chrom, pos, ref, alt, qual = variant
            vcf.write(f"{chrom}\t{pos+1}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")
