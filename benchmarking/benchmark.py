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
    varscan_cmd = f"java -jar VarScan.jar mpileup2snp {input_mpileup} --min-var-frequency 0.2 --min-freq-for-hom 0.8 --p-value 0.01 --output-vcf 1 --variants-only > {output_vcf}"
    subprocess.run(varscan_cmd, shell=True)

def run_snv_caller(input_mpileup, output_vcf):
    """Run SNV Caller tool."""
    snv_caller_cmd = f"python3 mpileup.py {input_mpileup} {output_vcf}"
    subprocess.run(snv_caller_cmd, shell=True)

def count_variants(vcf_file):
    """Count the number of variants in a VCF file."""
    try:
        with open(vcf_file, 'r') as f:
            variants = sum(1 for line in f if not line.startswith('#'))
    except FileNotFoundError:
        variants = 0
    return variants

def compare_tools(input_mpileup, output_vcf_varscan, output_vcf_snvcaller):
    """Compare VarScan and SNV Caller in terms of number of variants."""
    # Run VarScan
    run_varscan(input_mpileup, output_vcf_varscan)

    # Run SNV Caller
    run_snv_caller(input_mpileup, output_vcf_snvcaller)

    # Count variants detected by each tool
    varscan_variants = count_variants(output_vcf_varscan)
    snvcaller_variants = count_variants(output_vcf_snvcaller)

    return {
        'VarScan': {'variants': varscan_variants},
        'SNVCaller': {'variants': snvcaller_variants}
    }

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
    output_vcf_snvcaller = "snvcaller_output.vcf"

    # Compare tools
    results = compare_tools(input_mpileup, output_vcf_varscan, output_vcf_snvcaller)

    # Visualize results
    visualize_results(results)
