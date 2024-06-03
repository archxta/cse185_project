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
positions1 = set(extract_variant_positions(vcf_file1))
positions2 = set(extract_variant_positions(vcf_file2))

# Find variants common to both files
common_positions = positions1.intersection(positions2)

# Create scatter plot
plt.figure(figsize=(8, 6))

# Plot variants found in both files (purple)
plt.scatter(common_positions, common_positions, color='purple', alpha=0.5, label='Variants in both files')

# Plot variants found only in output.vcf (blue)
plt.scatter(positions1.difference(positions2), positions1.difference(positions2), color='blue', alpha=0.5, label='Variants only in output.vcf')

# Plot variants found only in output2.vcf (red)
plt.scatter(positions2.difference(positions1), positions2.difference(positions1), color='red', alpha=0.5, label='Variants only in output2.vcf')

plt.xlabel('Variant positions in output.vcf')
plt.ylabel('Variant positions in output2.vcf')
plt.title('Scatter Plot of Variant Positions')

# Add diagonal line with slope 1
max_position = max(max(positions1), max(positions2))
plt.plot([0, max_position], [0, max_position], color='black', linestyle='--')

# Add legend
plt.legend()

# Save scatter plot to a file
plt.savefig('scatter_plot.png')

# Show the plot
plt.show()
