import os
import pybwa

def index_reference_genome(reference_fasta):
    bwa_index = pybwa.Index(reference_fasta)
    bwa_index.build_index()
    return bwa_index

