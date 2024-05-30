import os

def index_reference_genome(reference_fasta):
    os.system(f'bwa index {reference_fasta}')
    os.system(f'samtools faidx {reference_fasta}')
    print(f"Indexed {reference_fasta}")
