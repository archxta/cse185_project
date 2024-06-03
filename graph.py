import vcf
import matplotlib.pyplot as plt

def extract_variant_positions(vcf_file):
    positions = []
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        positions.append(record.POS)
    return positions

# Path to the VCF files
vcf_file1 = 'output.vcf'
vcf_file2 = 'output2.vcf'

# Extract variant positions from VCF files
positions1 = extract_variant_positions(vcf_file1)
positions2 = extract_variant_positions(vcf_file2)

# Ensure both lists have the same length
min_length = min(len(positions1), len(positions2))
positions1 = positions1[:min_length]
positions2 = positions2[:min_length]

# Create scatter plot
plt.figure(figsize=(8, 6))
plt.scatter(positions1, positions2, color='blue', alpha=0.5)
plt.xlabel('Variant positions in output.vcf')
plt.ylabel('Variant positions in output2.vcf')
plt.title('Scatter Plot of Variant Positions')
plt.grid(True)

# Save scatter plot to a file
plt.savefig('scatter_plot.png')

# Show the plot
plt.show()
