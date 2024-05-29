import pybwa

def align_reads_to_reference(bwa_index, reads_fastq):
    alignment = bwa_index.align(reads_fastq)
    return alignment

