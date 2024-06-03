import matplotlib.pyplot as plt

# Example data from your VCF file
variant_positions = [53]
allele_frequencies = [0.5]  # Assuming a heterozygous variant with 50% allele frequency

# Plotting
plt.figure(figsize=(8, 6))
plt.bar(variant_positions, allele_frequencies, color='blue', alpha=0.5)
plt.xlabel('Variant Position')
plt.ylabel('Allele Frequency')
plt.title('Allele Frequencies of Variants')
plt.xticks(variant_positions)

# Save the plot as an image file
plt.savefig('allele_frequency_plot.png')
