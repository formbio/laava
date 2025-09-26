#!/usr/bin/env python3
"""Memory usage comparison script for get_flipflop_config.py versions."""

import psutil
import subprocess
import sys
import time
import argparse
from pathlib import Path


def monitor_memory_usage(command, description):
    """Monitor memory usage of a command and return peak memory usage."""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Command: {' '.join(command)}")
    print(f"{'='*60}")
    
    # Start the process
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Monitor memory usage
    peak_memory = 0
    memory_samples = []
    start_time = time.time()
    
    try:
        while process.poll() is None:
            try:
                # Get process and all its children
                parent = psutil.Process(process.pid)
                children = parent.children(recursive=True)
                
                # Calculate total memory usage (RSS - Resident Set Size)
                total_memory = parent.memory_info().rss
                for child in children:
                    try:
                        total_memory += child.memory_info().rss
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                
                # Convert to MB
                memory_mb = total_memory / (1024 * 1024)
                memory_samples.append(memory_mb)
                
                if memory_mb > peak_memory:
                    peak_memory = memory_mb
                
                # Print progress every 30 seconds
                elapsed = time.time() - start_time
                if len(memory_samples) % 30 == 0:  # Assuming 1 sample per second
                    print(f"  Progress: {elapsed:.0f}s elapsed, current memory: {memory_mb:.1f} MB, peak: {peak_memory:.1f} MB")
                
                time.sleep(1)  # Sample every second
                
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break
    
    except KeyboardInterrupt:
        print("\nInterrupted by user")
        process.terminate()
        return None, None, None
    
    # Wait for process to complete
    stdout, stderr = process.communicate()
    end_time = time.time()
    
    # Calculate statistics
    total_time = end_time - start_time
    avg_memory = sum(memory_samples) / len(memory_samples) if memory_samples else 0
    
    print(f"\nResults for {description}:")
    print(f"  Exit code: {process.returncode}")
    print(f"  Total time: {total_time:.1f} seconds")
    print(f"  Peak memory: {peak_memory:.1f} MB")
    print(f"  Average memory: {avg_memory:.1f} MB")
    print(f"  Memory samples: {len(memory_samples)}")
    
    if process.returncode != 0:
        print(f"  STDERR: {stderr.decode()}")
    
    return peak_memory, avg_memory, total_time


def main():
    parser = argparse.ArgumentParser(description="Compare memory usage between original and optimized versions")
    parser.add_argument("per_read_tsv", help="Per read TSV file")
    parser.add_argument("tagged_bam", help="Tagged BAM file")
    parser.add_argument("vector_type", help="Vector type (sc or ss)")
    parser.add_argument("orientation", help="Orientation (left or right)")
    parser.add_argument("-o", "--output-prefix", help="Output prefix", required=True)
    parser.add_argument("--flipflop-fasta", help="Optional flip-flop FASTA file")
    parser.add_argument("--original-only", action="store_true", help="Run only original version")
    parser.add_argument("--optimized-only", action="store_true", help="Run only optimized version")
    
    args = parser.parse_args()
    
    # Check if files exist
    script_dir = Path(__file__).parent
    original_script = script_dir / "get_flipflop_config.py"
    optimized_script = script_dir / "get_flipflop_config_optimized.py"
    
    if not original_script.exists():
        print(f"Error: Original script not found at {original_script}")
        sys.exit(1)
    
    if not optimized_script.exists():
        print(f"Error: Optimized script not found at {optimized_script}")
        sys.exit(1)
    
    # Prepare base command arguments
    base_args = [
        args.per_read_tsv,
        args.tagged_bam,
        args.vector_type,
        args.orientation,
        "-o", args.output_prefix
    ]
    
    if args.flipflop_fasta:
        base_args.extend(["--flipflop-fasta", args.flipflop_fasta])
    
    results = {}
    
    # Run original version
    if not args.optimized_only:
        original_cmd = ["python3", str(original_script)] + base_args
        original_output_prefix = args.output_prefix + "_original"
        original_cmd[-1] = original_output_prefix  # Replace output prefix
        
        peak_mem, avg_mem, total_time = monitor_memory_usage(original_cmd, "Original Version")
        if peak_mem is not None:
            results['original'] = {
                'peak_memory': peak_mem,
                'avg_memory': avg_mem,
                'total_time': total_time
            }
    
    # Run optimized version
    if not args.original_only:
        optimized_cmd = ["python3", str(optimized_script)] + base_args
        optimized_output_prefix = args.output_prefix + "_optimized"
        optimized_cmd[-1] = optimized_output_prefix  # Replace output prefix
        
        peak_mem, avg_mem, total_time = monitor_memory_usage(optimized_cmd, "Optimized Version")
        if peak_mem is not None:
            results['optimized'] = {
                'peak_memory': peak_mem,
                'avg_memory': avg_mem,
                'total_time': total_time
            }
    
    # Print comparison
    if len(results) == 2:
        print(f"\n{'='*60}")
        print("COMPARISON RESULTS")
        print(f"{'='*60}")
        
        orig = results['original']
        opt = results['optimized']
        
        print(f"Peak Memory Usage:")
        print(f"  Original:  {orig['peak_memory']:.1f} MB")
        print(f"  Optimized: {opt['peak_memory']:.1f} MB")
        print(f"  Reduction: {orig['peak_memory'] - opt['peak_memory']:.1f} MB ({((orig['peak_memory'] - opt['peak_memory']) / orig['peak_memory'] * 100):.1f}%)")
        
        print(f"\nAverage Memory Usage:")
        print(f"  Original:  {orig['avg_memory']:.1f} MB")
        print(f"  Optimized: {opt['avg_memory']:.1f} MB")
        print(f"  Reduction: {orig['avg_memory'] - opt['avg_memory']:.1f} MB ({((orig['avg_memory'] - opt['avg_memory']) / orig['avg_memory'] * 100):.1f}%)")
        
        print(f"\nExecution Time:")
        print(f"  Original:  {orig['total_time']:.1f} seconds")
        print(f"  Optimized: {opt['total_time']:.1f} seconds")
        time_diff = opt['total_time'] - orig['total_time']
        time_pct = (time_diff / orig['total_time'] * 100)
        if time_diff > 0:
            print(f"  Difference: +{time_diff:.1f} seconds (+{time_pct:.1f}% slower)")
        else:
            print(f"  Difference: {time_diff:.1f} seconds ({time_pct:.1f}% faster)")


if __name__ == "__main__":
    main()
