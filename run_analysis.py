#!/usr/bin/env python
"""Run all analysis notebooks and produce a consolidated cross-model summary.

Usage:
  python run_analysis.py                        # run notebooks + print summary
  python run_analysis.py --summary-only         # skip notebooks, just print summary
  python run_analysis.py Hubbard_Atom Dimer     # only these models
  python run_analysis.py --no-execute           # summary only (alias for --summary-only)
"""

import argparse
import os
import subprocess
import sys

import numpy as np
from glob import glob
from os.path import basename, isdir, isfile, join

from h5 import HDFArchive
from triqs.gf import inverse


# =====================================================================
# Helpers (self-contained to avoid model.py imports at top level)
# =====================================================================

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))


def discover_models(selected=None):
    """Return list of model directory names that have results/*.h5 files."""
    models = []
    for d in sorted(os.listdir(REPO_ROOT)):
        full = join(REPO_ROOT, d)
        if not isdir(full) or not isdir(join(full, 'results')):
            continue
        if not glob(join(full, 'results', '*.h5')):
            continue
        if selected and d not in selected:
            continue
        models.append(d)
    return models


def load_results(model_dir):
    """Load all solver results from a model's results/ directory."""
    data = {}
    for fpath in sorted(glob(join(model_dir, 'results', '*.h5'))):
        solver = basename(fpath).replace('.h5', '')
        with HDFArchive(fpath, 'r') as ar:
            data[solver] = {}
            for key in ar:
                data[solver][key] = ar[key]
    return data


def get_block_list(data):
    """Extract block names from the first solver that has G."""
    for solver, obs in data.items():
        if 'G' in obs:
            return list(obs['G'].indices)
    return []


def compute_sigma(data, model_dir):
    """Compute self-energy for all solvers. Imports G0_iw from the model."""
    orig_path = sys.path.copy()
    orig_cwd = os.getcwd()
    # model.py uses os.getcwd() + '/../common' to find util.py,
    # so cwd must be the model directory (matching notebook behavior)
    sys.path.insert(0, model_dir)
    try:
        import importlib
        os.chdir(model_dir)
        spec = importlib.util.spec_from_file_location('model', join(model_dir, 'model.py'))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        G0_iw = mod.G0_iw
    except Exception:
        return {}
    finally:
        os.chdir(orig_cwd)
        sys.path[:] = orig_path

    sigma = {}
    for solver, obs in data.items():
        if 'G' not in obs:
            continue
        G = obs['G']
        S = G0_iw.copy()
        S << inverse(G0_iw) - inverse(G)
        sigma[solver] = S
    return sigma


def pairwise_max_dev(obs_dict, block_lst):
    """Compute pairwise max-norm deviations. Returns {(s1,s2): float}."""
    solvers = sorted(obs_dict.keys())
    devs = {}
    for i, s1 in enumerate(solvers):
        for j in range(i + 1, len(solvers)):
            s2 = solvers[j]
            dev = 0.0
            for bl in block_lst:
                diff = obs_dict[s1][bl].data - obs_dict[s2][bl].data
                dev = max(dev, np.max(np.abs(diff)))
            devs[(s1, s2)] = dev
    return devs


def collect_static_obs(data, obs_name='density'):
    """Collect static observable values across solvers."""
    result = {}
    for solver, obs in sorted(data.items()):
        if 'static_obs' not in obs or obs_name not in obs['static_obs']:
            continue
        val = obs['static_obs'][obs_name]
        if hasattr(val, 'items'):
            result[solver] = dict(val)
        else:
            result[solver] = val
    return result


# =====================================================================
# Notebook execution
# =====================================================================

def execute_notebooks(models):
    """Run jupyter nbconvert --execute on each model's notebook."""
    results = {}
    for model in models:
        nb_path = join(REPO_ROOT, model, 'notebook.ipynb')
        if not isfile(nb_path):
            print(f"  {model}: no notebook.ipynb, skipping")
            results[model] = 'no_notebook'
            continue

        print(f"  {model}: executing notebook ...", end=' ', flush=True)
        try:
            proc = subprocess.run(
                ['jupyter', 'nbconvert', '--execute', '--inplace', nb_path],
                capture_output=True, text=True, timeout=600,
                cwd=join(REPO_ROOT, model)
            )
            if proc.returncode == 0:
                print("OK")
                results[model] = 'ok'
            else:
                print("FAILED")
                # Print last few lines of stderr for diagnosis
                err_lines = proc.stderr.strip().splitlines()
                for line in err_lines[-5:]:
                    print(f"    {line}")
                results[model] = 'failed'
        except subprocess.TimeoutExpired:
            print("TIMEOUT")
            results[model] = 'timeout'

    return results


# =====================================================================
# Summary
# =====================================================================

def print_summary(models):
    """Load results for each model and print consolidated summary."""

    sep = "=" * 78
    print(f"\n{sep}")
    print("  BENCHMARK ANALYSIS SUMMARY")
    print(sep)

    all_model_metrics = {}

    for model in models:
        model_dir = join(REPO_ROOT, model)
        data = load_results(model_dir)
        if not data:
            continue

        solvers = sorted(data.keys())
        block_lst = get_block_list(data)

        # -- G deviations --
        G_dict = {s: d['G'] for s, d in data.items() if 'G' in d}
        G_devs = pairwise_max_dev(G_dict, block_lst) if len(G_dict) >= 2 else {}

        # -- Sigma deviations --
        Sigma_dict = compute_sigma(data, model_dir)
        Sigma_devs = pairwise_max_dev(Sigma_dict, block_lst) if len(Sigma_dict) >= 2 else {}

        # -- Static observables --
        density = collect_static_obs(data, 'density')

        # -- Observable availability --
        obs_keys = set()
        for s, d in data.items():
            obs_keys.update(k for k in d if k not in ('static_obs', 'Solver_Info'))

        metrics = {
            'solvers': solvers,
            'observables': sorted(obs_keys),
            'G_max_dev': max(G_devs.values()) if G_devs else None,
            'G_devs': G_devs,
            'Sigma_max_dev': max(Sigma_devs.values()) if Sigma_devs else None,
            'Sigma_devs': Sigma_devs,
            'density': density,
        }
        all_model_metrics[model] = metrics

    # -- Print per-model summaries --
    for model, m in all_model_metrics.items():
        print(f"\n--- {model} ---")
        print(f"  Solvers ({len(m['solvers'])}): {', '.join(m['solvers'])}")
        print(f"  Observables: {', '.join(m['observables'])}")

        if m['G_max_dev'] is not None:
            print(f"  G  max deviation: {m['G_max_dev']:.2e}")
            worst = max(m['G_devs'], key=m['G_devs'].get)
            print(f"     worst pair:    {worst[0]} vs {worst[1]}")
        if m['Sigma_max_dev'] is not None:
            print(f"  Sigma max dev:    {m['Sigma_max_dev']:.2e}")
            worst = max(m['Sigma_devs'], key=m['Sigma_devs'].get)
            print(f"     worst pair:    {worst[0]} vs {worst[1]}")

        if m['density']:
            print(f"  Density:")
            for solver, val in sorted(m['density'].items()):
                if isinstance(val, dict):
                    items = ", ".join(f"{k}: {v:.6f}" for k, v in val.items())
                    print(f"    {solver:20s}  {items}")
                else:
                    print(f"    {solver:20s}  {val:.6f}")

    # -- Cross-model overview table --
    print(f"\n{sep}")
    print("  CROSS-MODEL OVERVIEW")
    print(sep)
    hdr_model = "Model".ljust(28)
    hdr_solvers = "#Solv"
    hdr_g = "G max dev"
    hdr_s = "Sig max dev"
    print(f"  {hdr_model} {hdr_solvers}  {hdr_g:>11s}  {hdr_s:>11s}")
    print(f"  {'-'*28} {'-'*5}  {'-'*11}  {'-'*11}")
    for model, m in all_model_metrics.items():
        n = len(m['solvers'])
        g = f"{m['G_max_dev']:.2e}" if m['G_max_dev'] is not None else "n/a"
        s = f"{m['Sigma_max_dev']:.2e}" if m['Sigma_max_dev'] is not None else "n/a"
        print(f"  {model:28s} {n:5d}  {g:>11s}  {s:>11s}")

    print()
    return all_model_metrics


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(description="Run analysis notebooks and print summary")
    parser.add_argument('models', nargs='*', help="Models to analyze (default: all with results)")
    parser.add_argument('--summary-only', '--no-execute', action='store_true',
                        help="Skip notebook execution, only compute and print summary")
    args = parser.parse_args()

    selected = args.models or None
    models = discover_models(selected)

    if not models:
        print("No models with results found.")
        sys.exit(1)

    print(f"Models: {', '.join(models)}\n")

    # Step 1: Execute notebooks
    if not args.summary_only:
        print("=== Executing notebooks ===")
        nb_results = execute_notebooks(models)
        ok = sum(1 for v in nb_results.values() if v == 'ok')
        total = sum(1 for v in nb_results.values() if v != 'no_notebook')
        print(f"\nNotebooks: {ok}/{total} executed successfully")

    # Step 2: Consolidated summary
    print_summary(models)


if __name__ == '__main__':
    main()
