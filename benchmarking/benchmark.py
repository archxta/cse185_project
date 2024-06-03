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
    varscan_cmd = f"varscan mpileup2snp {input_mpileup} --output-vcf {output_vcf}"
    subprocess.run(varscan_cmd, shell=True)

def run_snv_caller(input_mpileup, output_vcf):
    """Run SNV Caller tool."""
    snv_caller_cmd = f"python snv_caller.py {input_mpileup} {output_vcf}"
    subprocess.run(snv_caller_cmd, shell=True)

def compare_tools(input_mpileup, output_vcf_varscan, output_vcf_snvcaller):
    """Compare VarScan and SNV Caller in terms of runtime."""
    # Benchmark VarScan
    varscan_runtime = benchmark_tool(f"varscan mpileup2snp {input_mpileup} --output-vcf {output_vcf_varscan}")

    # Benchmark SNV Caller
    snvcaller_runtime = benchmark_tool(f"python snv_caller.py {input_mpileup} {output_vcf_snvcaller}")

    return {
        'VarScan': {'runtime': varscan_runtime},
        'SNVCaller': {'runtime': snvcaller_runtime}
    }

def visualize_results(results):
    """Visualize benchmarking results."""
    tools = list(results.keys())
    runtimes = [results[tool]['runtime'] for tool in tools]

    plt.bar(tools, runtimes, color=['blue', 'green'])
    plt.title('Runtime Comparison (seconds)')
    plt.xlabel('Tool')
    plt.ylabel('Runtime (seconds)')
    plt.show()

if __name__ == "__main__":
    input_mpileup = "NA12878_child.mpileup"
    output_vcf_varscan = "varscan_output.vcf"
    output_vcf_snvcaller = "snvcaller_output.vcf"

    # Run VarScan
    run_varscan(input_mpileup, output_vcf_varscan)

    # Run SNV Caller
    run_snv_caller(input_mpileup, output_vcf_snvcaller)

    # Compare tools
    results = compare_tools(input_mpileup, output_vcf_varscan, output_vcf_snvcaller)

    # Visualize results
    visualize_results(results)
