import matplotlib.pyplot as plt

# Define the sequences
sequence1 = "AACTTCATCCACGTTCACCTTGCCCCACAGGGCAGTAACGGCAGACTTCTCCACAGGAGTC"
sequence2 = "AACTTCATCCACGTTCACCTTGCCCCACAGGGCAGTAACGGCAGACTTCTCCTCAGGAGTC"

# Define x coordinates for each nucleotide
x1 = range(len(sequence1))
x2 = range(len(sequence2))

# Plotting
plt.figure(figsize=(15, 2))
plt.scatter(x1, [1] * len(x1), c='blue', label='Sequence 1', s=100)
plt.scatter(x2, [2] * len(x2), c='red', label='Sequence 2', s=100)
plt.yticks([])
plt.xlabel('Position')
plt.title('Sequence Comparison')
plt.legend()
plt.show()
