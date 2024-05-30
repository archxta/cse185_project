import os

def align_reads_to_reference(reference_fasta, reads_fastq, output_sam):
    os.system(f'bwa mem {reference_fasta} {reads_fastq} > {output_sam}')
    print(f"Aligned {reads_fastq} to {reference_fasta}")
