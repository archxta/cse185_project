import subprocess
import time
import matplotlib.pyplot as plt

def benchmark_tool(tool_cmd):
    """Benchmark a given command line tool for runtime."""
    start_time = time.time()
    subprocess.run(tool_cmd, shell=True)
    runtime = time.time() - start_time
    return runtime

def run_varscan(input_mpileup, output_vcf):
    """Run VarScan tool."""
    varscan_cmd = (
        f"java -jar VarScan.jar mpileup2snp {input_mpileup} "
        f"--min-var-frequency 0.2 --min-freq-for-hom 0.8 --p-value 0.01 "
        f"--output-vcf 1 --variants-only > {output_vcf}"
    )
    subprocess.run(varscan_cmd, shell=True)

def run_mpileup_py(input_mpileup, output_vcf, min_var_freq=0.2, min_avg_qual=15, p_value_thresh=0.01):
    """Run mpileup.py tool."""
    mpileup_py_cmd = (
        f"python3 mpileup.py {input_mpileup} {output_vcf} "
        f"{min_var_freq} {min_avg_qual} {p_value_thresh}"
    )
    subprocess.run(mpileup_py_cmd, shell=True)

def count_variants(vcf_file):
    """Count the number of variants in a VCF file."""
    try:
        with open(vcf_file, 'r') as f:
            variants = sum(1 for line in f if not line.startswith('#'))
    except FileNotFoundError:
        variants = 0
    return variants

def compare_tools(input_mpileup, output_vcf_varscan, output_vcf_mpileup_py, min_var_freq=0.2, min_avg_qual=15, p_value_thresh=0.01):
    """Compare VarScan and mpileup.py in terms of number of variants."""
    # Run VarScan
    run_varscan(input_mpileup, output_vcf_varscan)

    # Run mpileup.py
    run_mpileup_py(input_mpileup, output_vcf_mpileup_py, min_var_freq, min_avg_qual, p_value_thresh)

    # Count variants detected by each tool
    varscan_variants = count_variants(output_vcf_varscan)
    mpileup_py_variants = count_variants(output_vcf_mpileup_py)

    return {
        'VarScan': {'variants': varscan_variants},
        'mpileup.py': {'variants': mpileup_py_variants}
    }

def save_results_to_file(results, filename):
    """Save results to a file."""
    with open(filename, 'w') as f:
        for tool, data in results.items():
            f.write(f"{tool}:\n")
            f.write(f"  Variants: {data['variants']}\n")
            f.write("\n")

def visualize_results(results):
    """Visualize benchmarking results."""
    tools = list(results.keys())
    variants = [results[tool]['variants'] for tool in tools]

    plt.bar(tools, variants, color=['blue', 'green'])
    plt.title('Number of Variants Detected')
    plt.xlabel('Tool')
    plt.ylabel('Number of Variants')
    plt.show()

if __name__ == "__main__":
    input_mpileup = "NA12878_child.mpileup"
    output_vcf_varscan = "varscan_output.vcf"
    output_vcf_mpileup_py = "mpileup_py_output.vcf"
    results_file = "comparison_results.txt"

    # Compare tools
    results = compare_tools(input_mpileup, output_vcf_varscan, output_vcf_mpileup_py)

    # Save results to a file
    save_results_to_file(results, results_file)

    # Visualize results
    visualize_results(results)
