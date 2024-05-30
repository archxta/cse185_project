import os

def index_bam_file(bam_file):
    os.system(f'samtools index {bam_file}')
    print(f"Indexed {bam_file}")
