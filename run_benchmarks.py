#!/usr/bin/env python
"""Master benchmark runner.

Reads benchmark_config.yaml and translates each run entry into a shell command.

Usage:
  python run_benchmarks.py                          # all models, all runs
  python run_benchmarks.py Hubbard_Atom             # one model, all its runs
  python run_benchmarks.py Hubbard_Atom ctint       # one model, one solver
  python run_benchmarks.py --solver ctint           # one solver across all models
  python run_benchmarks.py --dry-run                # print commands without running
  python run_benchmarks.py --slurm                  # generate SLURM job scripts
"""

import argparse
import json
import os
import subprocess
import sys
import time

import yaml


def load_config(config_path='benchmark_config.yaml'):
    with open(config_path) as f:
        return yaml.safe_load(f)


def parse_run_entry(entry):
    """Parse a run entry from the config into (solver_name, measure_list)."""
    if isinstance(entry, str):
        return entry, []
    elif isinstance(entry, dict):
        solver = entry['solver']
        measure = entry.get('measure', [])
        if isinstance(measure, str):
            measure = [measure]
        return solver, measure
    else:
        raise ValueError(f"Invalid run entry: {entry}")


def build_command(model, solver, measure, n_mpi_ranks):
    """Build the shell command for a single solver run."""
    script_dir = os.path.join(model, 'scripts')
    script_path = os.path.join(script_dir, solver)

    if not os.path.exists(script_path):
        return None, f"Script not found: {script_path}"

    cmd_parts = []
    if n_mpi_ranks > 1:
        cmd_parts.extend(['mpirun', '-np', str(n_mpi_ranks)])
    cmd_parts.extend(['python', solver])

    if measure:
        cmd_parts.extend(['--measure'] + measure)

    cmd = ' '.join(cmd_parts)
    return cmd, None


def run_single(model, solver, measure, n_mpi_ranks, dry_run=False):
    """Run a single benchmark. Returns dict with status info."""
    cmd, err = build_command(model, solver, measure, n_mpi_ranks)
    result = {
        'model': model,
        'solver': solver,
        'measure': measure,
        'command': cmd,
    }

    if err:
        result['status'] = 'skipped'
        result['error'] = err
        return result

    script_dir = os.path.join(model, 'scripts')
    log_path = os.path.join(model, 'results', f'{solver}.log')

    if dry_run:
        print(f"  [DRY RUN] cd {script_dir} && {cmd}")
        result['status'] = 'dry_run'
        return result

    print(f"  Running: cd {script_dir} && {cmd}")

    # Ensure results directory exists
    os.makedirs(os.path.join(model, 'results'), exist_ok=True)

    start = time.time()
    try:
        with open(log_path, 'w') as logf:
            proc = subprocess.run(
                cmd, shell=True, cwd=script_dir,
                stdout=logf, stderr=subprocess.STDOUT,
                timeout=3600  # 1 hour timeout
            )
        elapsed = time.time() - start
        result['run_time'] = elapsed
        result['status'] = 'pass' if proc.returncode == 0 else 'fail'
        result['returncode'] = proc.returncode
        if proc.returncode != 0:
            print(f"    FAILED (exit {proc.returncode}, {elapsed:.1f}s) -- see {log_path}")
        else:
            print(f"    OK ({elapsed:.1f}s)")
    except subprocess.TimeoutExpired:
        result['status'] = 'timeout'
        print(f"    TIMEOUT after 3600s -- see {log_path}")
    except Exception as e:
        result['status'] = 'error'
        result['error'] = str(e)
        print(f"    ERROR: {e}")

    return result


def generate_slurm_script(model, solver, measure, n_mpi_ranks):
    """Generate a SLURM job script for a single run."""
    cmd, err = build_command(model, solver, measure, n_mpi_ranks)
    if err:
        return None

    script_dir = os.path.join(model, 'scripts')
    job_name = f"bench_{model}_{solver}"

    slurm = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --ntasks={n_mpi_ranks}
#SBATCH --time=01:00:00
#SBATCH --output={model}/results/{solver}.log

cd {script_dir}
{cmd}
"""
    return slurm


def main():
    parser = argparse.ArgumentParser(description="Run solver benchmarks")
    parser.add_argument('model', nargs='?', help="Model to run (default: all)")
    parser.add_argument('solver', nargs='?', help="Solver to run (default: all for model)")
    parser.add_argument('--solver', dest='solver_flag',
                        help="Run this solver across all models that have it")
    parser.add_argument('--dry-run', action='store_true',
                        help="Print commands without running")
    parser.add_argument('--slurm', action='store_true',
                        help="Generate SLURM job scripts instead of running")
    parser.add_argument('--config', default='benchmark_config.yaml',
                        help="Config file path")
    args = parser.parse_args()

    config = load_config(args.config)
    defaults = config.get('defaults', {})
    default_ranks = defaults.get('n_mpi_ranks', 4)

    models_config = config.get('models', {})

    # Filter models
    if args.model:
        if args.model not in models_config:
            print(f"Error: model '{args.model}' not in config")
            sys.exit(1)
        models_config = {args.model: models_config[args.model]}

    all_results = []

    for model, model_cfg in models_config.items():
        n_mpi_ranks = model_cfg.get('n_mpi_ranks', default_ranks)
        runs = model_cfg.get('runs', [])

        print(f"\n=== {model} ===")

        for entry in runs:
            solver, measure = parse_run_entry(entry)

            # Filter by solver if specified
            if args.solver and solver != args.solver:
                continue
            if args.solver_flag and solver != args.solver_flag:
                continue

            if args.slurm:
                script = generate_slurm_script(model, solver, measure, n_mpi_ranks)
                if script:
                    slurm_path = f"slurm_{model}_{solver}.sh"
                    with open(slurm_path, 'w') as f:
                        f.write(script)
                    print(f"  Generated {slurm_path}")
            else:
                result = run_single(model, solver, measure, n_mpi_ranks,
                                    dry_run=args.dry_run)
                all_results.append(result)

    # Summary
    if all_results and not args.dry_run:
        print("\n=== Summary ===")
        passed = sum(1 for r in all_results if r['status'] == 'pass')
        failed = sum(1 for r in all_results if r['status'] == 'fail')
        skipped = sum(1 for r in all_results if r['status'] == 'skipped')
        errors = sum(1 for r in all_results if r['status'] in ('error', 'timeout'))
        total = len(all_results)
        print(f"  {passed}/{total} passed, {failed} failed, "
              f"{skipped} skipped, {errors} errors")

        # Write summary JSON
        summary_path = 'benchmark_results.json'
        with open(summary_path, 'w') as f:
            json.dump(all_results, f, indent=2, default=str)
        print(f"  Details written to {summary_path}")

        if failed or errors:
            sys.exit(1)


if __name__ == '__main__':
    main()
