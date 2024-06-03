import subprocess
import time
import memory_profiler
import matplotlib.pyplot as plt

def benchmark_tool(tool_cmd):
    """Benchmark a given command line tool for runtime and memory usage."""
    start_time = time.time()
    mem_usage = memory_profiler.memory_usage((subprocess.run, [tool_cmd]), interval=0.1)
    runtime = time.time() - start_time
    max_memory = max(mem_usage)
    return runtime, max_memory

def run_varscan(input_mpileup, output_vcf):
    """Run VarScan tool."""
    varscan_cmd = f"varscan mpileup2snp {input_mpileup} --output-vcf {output_vcf}"
    subprocess.run(varscan_cmd, shell=True)

def run_snv_caller(input_mpileup, output_vcf):
    """Run SNV Caller tool."""
    snv_caller_cmd = f"python snv_caller.py {input_mpileup} {output_vcf}"
    subprocess.run(snv_caller_cmd, shell=True)

def compare_tools(input_mpileup, output_vcf_varscan, output_vcf_snvcaller):
    """Compare VarScan and SNV Caller in terms of runtime, memory usage, and number of variants."""
    # Benchmark VarScan
    varscan_runtime, varscan_memory = benchmark_tool(f"varscan mpileup2snp {input_mpileup} --output-vcf {output_vcf_varscan}")

    # Benchmark SNV Caller
    snvcaller_runtime, snvcaller_memory = benchmark_tool(f"python snv_caller.py {input_mpileup} {output_vcf_snvcaller}")

    return {
        'VarScan': {'runtime': varscan_runtime, 'memory': varscan_memory},
        'SNVCaller': {'runtime': snvcaller_runtime, 'memory': snvcaller_memory}
    }

def visualize_results(results):
    """Visualize benchmarking results."""
    tools = list(results.keys())
    runtimes = [results[tool]['runtime'] for tool in tools]
    memories = [results[tool]['memory'] for tool in tools]

    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    axs[0].bar(tools, runtimes, color=['blue', 'green'])
    axs[0].set_title('Runtime (seconds)')
    axs[1].bar(tools, memories, color=['blue', 'green'])
    axs[1].set_title('Memory Usage (MB)')

    for ax in axs.flat:
        ax.set_ylabel('Value')
        ax.set_xlabel('Tool')

    plt.tight_layout()
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
