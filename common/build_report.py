#!/usr/bin/env python
"""Generate a static per-model comparison report.

Usage (from repo root):
  python common/build_report.py                     # all models with results
  python common/build_report.py Hubbard_Atom Dimer  # selected models

For each model the script writes:
  MODEL/results/figures/*.png   -- one figure per observable
  MODEL/results/tables/*.md     -- deviation tables, static obs, solver info
  MODEL/results/report.md       -- assembled report (entry point)

The model description + LaTeX Hamiltonian is pulled from model.py's
module docstring. Parameter values are pulled from a whitelist of
module-level numeric attributes.
"""

import argparse
import importlib.util
import os
import sys
import traceback
from glob import glob
from os.path import isdir, isfile, join

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COMMON_DIR = join(REPO_ROOT, "common")
sys.path.insert(0, COMMON_DIR)

from analysis import (  # noqa: E402
    load_all_results,
    compute_sigma,
    deviation_table,
    plot_iw_comparison,
    plot_iw_scalar_comparison,
    plot_w_comparison,
    plot_chi_contour,
    plot_static_obs_table,
)


# --- Parameters rendered in the report's "Parameters" section ---------
PARAMETER_WHITELIST = [
    "beta", "n_iw", "n_w", "broadening",
    "U", "J", "mu", "h",
    "t", "t_perp", "t_prime",
    "w0", "D_coupling",
    "n_orb", "n_orb_bath",
    "block_names",
]

CHANNELS = ["d", "m", "s", "t"]


# =====================================================================
# Model import
# =====================================================================

def import_model(model_dir):
    """Import a model's model.py with cwd set to MODEL/ so its
    `sys.path.append(os.getcwd()+'/../common')` trick works.

    Returns the imported module object. Caller is responsible for
    restoring cwd and sys.path.
    """
    sys.modules.pop("model", None)
    sys.modules.pop("util", None)
    spec = importlib.util.spec_from_file_location(
        "model", join(model_dir, "model.py")
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"could not build spec for {model_dir}/model.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["model"] = mod
    spec.loader.exec_module(mod)
    return mod


# =====================================================================
# Parameter extraction
# =====================================================================

def _format_value(v):
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return f"{float(v):g}"
    if isinstance(v, (list, tuple)):
        return repr(list(v))
    if isinstance(v, np.ndarray):
        return f"array(shape={v.shape}, dtype={v.dtype})"
    return str(v)


def render_parameters(model):
    """Return a markdown parameter table pulled from model.py attributes."""
    rows = []
    for name in PARAMETER_WHITELIST:
        if not hasattr(model, name):
            continue
        val = getattr(model, name)
        if callable(val):
            continue
        rows.append((name, _format_value(val)))
    if not rows:
        return "_No parameters in whitelist found on model.py._\n"
    lines = ["| Name | Value |", "|------|-------|"]
    lines.extend(f"| `{n}` | {v} |" for n, v in rows)
    return "\n".join(lines) + "\n"


# =====================================================================
# Solver info table
# =====================================================================

def _info_get(info, key, default=None):
    """Read a field from a Solver_Info HDFArchiveGroup or dict."""
    if info is None:
        return default
    try:
        return info[key] if key in info else default
    except Exception:
        return default


def render_solver_info(data):
    """Markdown table of per-solver metadata from Solver_Info groups."""
    rows = []
    for solver in sorted(data.keys()):
        info = data[solver].get("Solver_Info") if hasattr(data[solver], "get") \
            else (data[solver]["Solver_Info"] if "Solver_Info" in data[solver] else None)
        git_hash = _info_get(info, "solver_git_hash", "") or ""
        rows.append({
            "solver": solver,
            "solver_name": _info_get(info, "solver_name", "?"),
            "solver_version": _info_get(info, "solver_version", "?"),
            "solver_git_hash": git_hash[:8],
            "triqs_version": _info_get(info, "triqs_version", "?"),
            "num_threads": _info_get(info, "num_threads", "?"),
            "run_time": _info_get(info, "run_time"),
        })
    lines = [
        "| solver | name | version | git | TRIQS | MPI | run time (s) |",
        "|--------|------|---------|-----|-------|-----|--------------|",
    ]
    for r in rows:
        t = r["run_time"]
        t_str = f"{t:.1f}" if isinstance(t, (int, float)) else "—"
        lines.append(
            f"| `{r['solver']}` | {r['solver_name']} | {r['solver_version']} | "
            f"`{r['solver_git_hash']}` | {r['triqs_version']} | "
            f"{r['num_threads']} | {t_str} |"
        )
    return "\n".join(lines) + "\n"


# =====================================================================
# Per-observable section builders
# =====================================================================

def _write(path, text):
    with open(path, "w") as f:
        f.write(text)


def _save_fig(fig, path):
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def section_g(data, block_lst, fig_dir, tab_dir):
    G = {s: d["G"] for s, d in data.items() if "G" in d}
    if not G:
        return None
    fig = plot_iw_comparison(G, "G", block_lst)
    _save_fig(fig, join(fig_dir, "G_iw.png"))
    md = deviation_table(G, block_lst, label="G", return_markdown=True)
    _write(join(tab_dir, "deviation_G_iw.md"), md)
    return f"## $G(i\\omega_n)$\n\n![](figures/G_iw.png)\n\n{md}\n"


def section_sigma(data, model, block_lst, fig_dir, tab_dir):
    if not hasattr(model, "G0_iw"):
        return None
    Sigma = compute_sigma(data, model.G0_iw)
    if not Sigma:
        return None
    fig = plot_iw_comparison(Sigma, r"$\Sigma$", block_lst)
    _save_fig(fig, join(fig_dir, "Sigma_iw.png"))
    md = deviation_table(Sigma, block_lst, label="Sigma", return_markdown=True)
    _write(join(tab_dir, "deviation_Sigma_iw.md"), md)
    return f"## $\\Sigma(i\\omega_n)$\n\n![](figures/Sigma_iw.png)\n\n{md}\n"


def section_g_w(data, block_lst, fig_dir, tab_dir):
    G_w = {s: d["G_w"] for s, d in data.items() if "G_w" in d}
    if not G_w:
        return None
    fig = plot_w_comparison(G_w, "G_w", block_lst, spectral=True)
    _save_fig(fig, join(fig_dir, "G_w.png"))
    md = deviation_table(G_w, block_lst, label="G(w)", return_markdown=True)
    _write(join(tab_dir, "deviation_G_w.md"), md)
    return f"## Spectral function $A(\\omega) = -\\tfrac{{1}}{{\\pi}} \\mathrm{{Im}}\\, G^R(\\omega)$\n\n![](figures/G_w.png)\n\n{md}\n"


def section_static_obs(data, tab_dir):
    pieces = []
    for obs in ("density", "nn_ab"):
        md = plot_static_obs_table(data, obs, return_markdown=True)
        if md is None:
            continue
        pieces.append(md)
    if not pieces:
        return None
    combined = "\n".join(pieces)
    _write(join(tab_dir, "static_obs.md"), combined)
    return "## Static observables\n\n" + combined + "\n"


def section_chi(data, prefix, fig_dir, tab_dir, plotter, latex_name):
    """Generic handler for chi2/chi3/chi4 — loops over physical channels."""
    any_section = False
    pieces = [f"## {latex_name}\n"]
    for ch in CHANNELS:
        key = f"{prefix}_{ch}"
        chi_dict = {s: d[key] for s, d in data.items() if key in d}
        if not chi_dict:
            continue
        any_section = True
        try:
            fig = plotter(chi_dict, ch)
        except Exception as e:
            pieces.append(f"### {key}\n\n_Plot failed: {e}_\n")
            continue
        img_name = f"{key}.png"
        _save_fig(fig, join(fig_dir, img_name))
        try:
            md = deviation_table(chi_dict, label=key, return_markdown=True)
        except Exception as e:
            md = f"_Deviation table failed: {e}_\n"
        _write(join(tab_dir, f"deviation_{key}.md"), md)
        pieces.append(f"### {key}\n\n![](figures/{img_name})\n\n{md}\n")
    return "".join(pieces) if any_section else None


def _chi2_plotter(chi_dict, channel):
    return plot_iw_scalar_comparison(chi_dict, f"chi2_{channel}")


def _chi3_plotter(chi_dict, channel):
    return plot_chi_contour(chi_dict, channel)


def _chi4_plotter(chi_dict, channel):
    return plot_chi_contour(chi_dict, channel)


# =====================================================================
# Report assembly
# =====================================================================

def _block_lst_from_data(data):
    for d in data.values():
        if "G" in d and hasattr(d["G"], "indices"):
            return list(d["G"].indices)
    return []


def assemble_report(model_name, model, data):
    block_lst = _block_lst_from_data(data)
    fig_dir = "results/figures"
    tab_dir = "results/tables"
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(tab_dir, exist_ok=True)

    # Header: prefer the docstring's own H1 if present, else use the folder name.
    doc = (model.__doc__ or "").strip()
    title = model_name
    body_doc = doc
    if doc.startswith("# "):
        first_line, _, rest = doc.partition("\n")
        title = first_line[2:].strip()
        body_doc = rest.lstrip()
    if not doc:
        body_doc = f"_(no module docstring in `{model_name}/model.py`)_"

    parameters_md = render_parameters(model)
    solver_info_md = render_solver_info(data)
    _write(join(tab_dir, "solver_info.md"), solver_info_md)

    sections = [
        section_g(data, block_lst, fig_dir, tab_dir),
        section_sigma(data, model, block_lst, fig_dir, tab_dir),
        section_static_obs(data, tab_dir),
        section_g_w(data, block_lst, fig_dir, tab_dir),
        section_chi(data, "chi2", fig_dir, tab_dir, _chi2_plotter, r"Two-point susceptibilities $\chi_2$"),
        section_chi(data, "chi3", fig_dir, tab_dir, _chi3_plotter, r"Three-point correlators $\chi_3$"),
        section_chi(data, "chi4", fig_dir, tab_dir, _chi4_plotter, r"Four-point correlators $\chi_4 / G_{2c}$"),
    ]
    body = "\n\n".join(s for s in sections if s)

    report = (
        f"# {title}\n\n"
        f"{body_doc}\n\n"
        f"## Parameters\n\n{parameters_md}\n"
        f"## Solvers\n\n{solver_info_md}\n"
        f"{body}\n"
    )
    _write("results/report.md", report)


# =====================================================================
# Per-model driver
# =====================================================================

def build_one(model_name):
    model_dir = join(REPO_ROOT, model_name)
    if not isdir(model_dir) or not isfile(join(model_dir, "model.py")):
        print(f"  {model_name}: not a model directory, skipping")
        return False
    if not glob(join(model_dir, "results", "*.h5")):
        print(f"  {model_name}: no results/*.h5, skipping")
        return False

    print(f"  {model_name}: building report ...", flush=True)
    original_cwd = os.getcwd()
    original_path = sys.path[:]
    try:
        os.chdir(model_dir)
        sys.path.insert(0, model_dir)
        model = import_model(model_dir)
        data = load_all_results("results")
        assemble_report(model_name, model, data)
        print(f"    -> {model_name}/results/report.md")
        return True
    except Exception as e:
        print(f"  {model_name}: FAILED ({type(e).__name__}: {e})")
        traceback.print_exc()
        return False
    finally:
        os.chdir(original_cwd)
        sys.path[:] = original_path


def discover_models():
    return sorted(
        d for d in os.listdir(REPO_ROOT)
        if isdir(join(REPO_ROOT, d))
        and isfile(join(REPO_ROOT, d, "model.py"))
        and glob(join(REPO_ROOT, d, "results", "*.h5"))
    )


def main():
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("models", nargs="*",
                        help="Model directories (default: all with results)")
    args = parser.parse_args()

    models = args.models or discover_models()
    if not models:
        print("No models with results/*.h5 found.")
        sys.exit(1)

    ok = 0
    for m in models:
        if build_one(m):
            ok += 1
    print(f"\nBuilt {ok}/{len(models)} reports.")


if __name__ == "__main__":
    main()
