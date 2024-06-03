import matplotlib.pyplot as plt

def parse_vcf(vcf_file):
    variant_positions = []
    allele_frequencies = []

    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):  # Skip header lines
                continue
            fields = line.strip().split('\t')
            position = int(fields[1])
            info = fields[7].split(';')
            allele_frequency = float([x.split('=')[1] for x in info if x.startswith('AF=')][0])
            variant_positions.append(position)
            allele_frequencies.append(allele_frequency)

    return variant_positions, allele_frequencies

# Parse and plot data from each VCF file
vcf_files = ['output1.vcf', 'output2.vcf', 'output3.vcf']
colors = ['blue', 'green', 'red']

plt.figure(figsize=(10, 6))
for i, vcf_file in enumerate(vcf_files):
    positions, frequencies = parse_vcf(vcf_file)
    plt.plot(positions, frequencies, label=f'File {i+1}', color=colors[i], marker='o')

plt.xlabel('Variant Position')
plt.ylabel('Allele Frequency')
plt.title('Allele Frequencies of Variants')
plt.legend()
plt.grid(True)
plt.savefig('allele_frequency_plot.png')
plt.show()
