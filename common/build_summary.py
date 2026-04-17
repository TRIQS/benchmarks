#!/usr/bin/env python
"""Cross-model / cross-solver correctness summary.

Walks every model registered in `benchmark_config.yaml`, computes the
max-norm deviation of each solver's observables vs. a per-model reference
(typically an ED solver), and writes three markdown files to `summary/`:

  summary/correctness.md  -- per-(model, observable, solver) deviation table
  summary/coverage.md     -- which solvers ran on which models, plus wall time
  summary/README.md       -- landing page with links to each MODEL/results/report.md

Usage (from repo root):
  python common/build_summary.py
"""

import importlib.util
import os
import sys
from os.path import isdir, isfile, join

import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, join(REPO_ROOT, "common"))

from analysis import (  # noqa: E402
    load_all_results,
    compute_sigma,
    deviation_vs_reference,
)

SUMMARY_DIR = join(REPO_ROOT, "summary")

# Observable ordering for the correctness table. "Sigma" is derived from G
# and the model's G0_iw — handled specially below.
OBSERVABLE_ORDER = [
    "G", "Sigma", "G_w",
    "chi2_d", "chi2_m", "chi2_s", "chi2_t",
    "chi3_d", "chi3_m", "chi3_s", "chi3_t",
    "chi4_d", "chi4_m", "chi4_s", "chi4_t",
]

# Fallback when benchmark_config.yaml does not set reference_solver.
DEFAULT_PRIORITY = ["exact", "pyed", "pomerol", "edipack", "atomdiag"]


# =====================================================================
# Config + model import
# =====================================================================

def load_config():
    with open(join(REPO_ROOT, "benchmark_config.yaml")) as f:
        return yaml.safe_load(f)


def import_model(model_dir):
    sys.modules.pop("model", None)
    sys.modules.pop("util", None)
    spec = importlib.util.spec_from_file_location(
        "model", join(model_dir, "model.py"))
    if spec is None or spec.loader is None:
        raise ImportError(f"bad spec for {model_dir}/model.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["model"] = mod
    spec.loader.exec_module(mod)
    return mod


def pick_reference(available, model_cfg, priority):
    """Choose the reference solver for a model."""
    override = model_cfg.get("reference_solver") if isinstance(model_cfg, dict) else None
    if override and override in available:
        return override
    for s in priority:
        if s in available:
            return s
    return None


# =====================================================================
# Per-model data collection
# =====================================================================

def _info_get(info, key, default=None):
    if info is None:
        return default
    try:
        return info[key] if key in info else default
    except Exception:
        return default


def collect_model(model_name):
    """Load results + derive Sigma for one model.

    Returns a dict with keys:
      solvers, observables ({name: {solver: gf}}), run_times ({solver: float|None})
    or None on failure.
    """
    model_dir = join(REPO_ROOT, model_name)
    if not isfile(join(model_dir, "model.py")):
        return None
    if not any(f.endswith(".h5") for f in os.listdir(join(model_dir, "results"))
               if isfile(join(model_dir, "results", f))):
        return None

    original_cwd = os.getcwd()
    original_path = sys.path[:]
    try:
        os.chdir(model_dir)
        sys.path.insert(0, model_dir)
        model = import_model(model_dir)
        data = load_all_results("results")
        if not data:
            return None

        # Pull Sigma in addition to the stored observables.
        sigma = compute_sigma(data, model.G0_iw) if hasattr(model, "G0_iw") else {}

        observables = {"Sigma": sigma} if sigma else {}
        for solver, obs in data.items():
            for key, val in obs.items():
                if key in ("static_obs", "Solver_Info"):
                    continue
                observables.setdefault(key, {})[solver] = val

        run_times = {}
        for solver in data:
            info = data[solver].get("Solver_Info") if hasattr(data[solver], "get") \
                else (data[solver]["Solver_Info"] if "Solver_Info" in data[solver] else None)
            run_times[solver] = _info_get(info, "run_time")

        return {
            "solvers": sorted(data.keys()),
            "observables": observables,
            "run_times": run_times,
        }
    except Exception as e:
        print(f"  {model_name}: FAILED ({type(e).__name__}: {e})")
        return None
    finally:
        os.chdir(original_cwd)
        sys.path[:] = original_path


# =====================================================================
# Markdown assembly
# =====================================================================

def _fmt_dev(x):
    return f"{x:.2e}" if x is not None else "—"


def build_correctness(per_model, refs):
    """One row per (model, observable) with the reference solver's column
    replaced by '(ref)'. Solvers without this observable show '—'.
    """
    all_solvers = sorted({s for info in per_model.values() for s in info["solvers"]})

    lines = [
        "# Correctness summary",
        "",
        "Max-norm deviation of each solver from the model's reference solver.",
        "Column `ref` lists which solver was chosen as reference. `—` means the",
        "solver did not compute this observable; `(ref)` marks the reference",
        "itself in its own column.",
        "",
        "| Model | Observable | ref | " + " | ".join(all_solvers) + " |",
        "|---|---|---|" + "---|" * len(all_solvers),
    ]

    for model in sorted(per_model):
        info = per_model[model]
        ref = refs[model]
        if ref is None:
            # No reference — skip deviation computation; emit informational row.
            lines.append(
                f"| {model} | _all_ | — | "
                + " | ".join("" for _ in all_solvers) + " |"
            )
            continue

        ordered = [o for o in OBSERVABLE_ORDER if o in info["observables"]]
        extras = [o for o in info["observables"] if o not in OBSERVABLE_ORDER]
        for obs_name in ordered + sorted(extras):
            obs_dict = info["observables"][obs_name]
            if ref not in obs_dict:
                continue
            try:
                devs = deviation_vs_reference(obs_dict, ref)
            except Exception as e:
                lines.append(f"| {model} | {obs_name} | {ref} | "
                             + " | ".join("_err_" for _ in all_solvers) + " |")
                print(f"  {model}/{obs_name}: deviation failed ({e})")
                continue

            cells = []
            for s in all_solvers:
                if s == ref:
                    cells.append("(ref)")
                elif s in devs:
                    cells.append(_fmt_dev(devs[s]))
                elif s in obs_dict:
                    cells.append(_fmt_dev(0.0))
                else:
                    cells.append("—")
            lines.append(f"| {model} | `{obs_name}` | `{ref}` | "
                         + " | ".join(cells) + " |")

    return "\n".join(lines) + "\n"


def build_coverage(per_model):
    """Matrix of wall times per (model, solver). Cell shows run_time or '✓'/'—'."""
    all_solvers = sorted({s for info in per_model.values() for s in info["solvers"]})

    lines = [
        "# Solver coverage",
        "",
        "Wall time (seconds) per solver/model. `✓` means the solver ran but",
        "did not report `run_time` in its `Solver_Info` group; `—` means the",
        "solver has no result for this model.",
        "",
        "| Model | " + " | ".join(all_solvers) + " |",
        "|---|" + "---|" * len(all_solvers),
    ]
    for model in sorted(per_model):
        info = per_model[model]
        row = []
        for s in all_solvers:
            if s not in info["solvers"]:
                row.append("—")
                continue
            t = info["run_times"].get(s)
            row.append(f"{t:.0f}" if isinstance(t, (int, float)) else "✓")
        lines.append(f"| {model} | " + " | ".join(row) + " |")
    return "\n".join(lines) + "\n"


def build_index(per_model, refs):
    lines = [
        "# Benchmark summary",
        "",
        "Cross-model correctness and coverage overview. See the individual",
        "`MODEL/results/report.md` files for full plots and deviation tables.",
        "",
        "- [correctness.md](correctness.md) -- max-norm deviation vs. reference solver",
        "- [coverage.md](coverage.md) -- solver/model wall-time matrix",
        "",
        "## Per-model reports",
        "",
    ]
    for model in sorted(per_model):
        ref = refs[model] or "_none_"
        n = len(per_model[model]["solvers"])
        lines.append(
            f"- [{model}](../{model}/results/report.md) "
            f"— {n} solvers, reference: `{ref}`"
        )
    lines.append("")
    return "\n".join(lines)


# =====================================================================
# Main
# =====================================================================

def main():
    cfg = load_config()
    defaults = cfg.get("defaults", {}) or {}
    priority = defaults.get("reference_solver_priority") or DEFAULT_PRIORITY
    models_cfg = cfg.get("models", {}) or {}

    # Accept extra model directories that have results but aren't in the config.
    discovered = [
        d for d in sorted(os.listdir(REPO_ROOT))
        if isdir(join(REPO_ROOT, d, "results"))
        and isfile(join(REPO_ROOT, d, "model.py"))
    ]
    model_list = [m for m in discovered if m in models_cfg] + \
                 [m for m in discovered if m not in models_cfg]

    per_model = {}
    refs = {}
    for m in model_list:
        print(f"  {m}: collecting ...", flush=True)
        info = collect_model(m)
        if info is None:
            continue
        per_model[m] = info
        refs[m] = pick_reference(info["solvers"], models_cfg.get(m, {}), priority)

    if not per_model:
        print("No models collected.")
        sys.exit(1)

    os.makedirs(SUMMARY_DIR, exist_ok=True)
    with open(join(SUMMARY_DIR, "correctness.md"), "w") as f:
        f.write(build_correctness(per_model, refs))
    with open(join(SUMMARY_DIR, "coverage.md"), "w") as f:
        f.write(build_coverage(per_model))
    with open(join(SUMMARY_DIR, "README.md"), "w") as f:
        f.write(build_index(per_model, refs))

    print(f"\nWrote summary/{{correctness,coverage,README}}.md "
          f"for {len(per_model)} models.")


if __name__ == "__main__":
    main()
