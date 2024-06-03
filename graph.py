from Bio import SeqIO
import matplotlib.pyplot as plt

# Define the path to the directory containing the variant FASTA files
fasta_directory = 'data/'

variant_files = ['variant1.fasta', 'variant2.fasta', 'variant3.fasta']
variant_sequences = {}

# Read variant FASTA files
for variant_file in variant_files:
    file_path = fasta_directory + variant_file
    with open(file_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            variant_sequences[variant_file] = str(record.seq)

# Identify variants and extract variant information
variant_positions = {}
for variant_file, sequence in variant_sequences.items():
    variant_positions[variant_file] = []
    for i, base in enumerate(sequence):
        # Check if base differs from reference (assumed reference is first sequence)
        if base != variant_sequences[variant_files[0]][i]:
            variant_positions[variant_file].append(i)

# Plot data
plt.figure(figsize=(10, 6))
for variant_file, positions in variant_positions.items():
    allele_frequencies = [len(positions) / len(variant_sequences[variant_file])] * len(positions)
    plt.scatter(positions, allele_frequencies, label=variant_file)

plt.xlabel('Variant Position')
plt.ylabel('Allele Frequency')
plt.title('Allele Frequencies of Variants')
plt.legend()
plt.grid(True)
plt.savefig('allele_frequency_plot.png')
plt.show()

