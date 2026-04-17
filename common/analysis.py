"""Shared analysis and plotting library for the static report generator.

Loads per-solver HDF5 results, derives quantities like the self-energy,
computes max-norm deviations (against a reference or pairwise), and
produces either matplotlib figures or markdown table strings for the
`common/build_report.py` / `common/build_summary.py` pipeline.
"""

import numpy as np
from glob import glob
from os.path import basename

from h5 import HDFArchive
from triqs.gf import BlockGf, inverse


# =====================================================================
# Data loading
# =====================================================================

def load_all_results(results_dir='results'):
    """Auto-discover and load all solver results from HDF5 archives.

    Returns
    -------
    dict : {solver_name: {key: value, ...}, ...}
        e.g. {'cthyb': {'G': BlockGf, ...}, 'ctint': {'G': ..., 'chi3_d': ...}}
    """
    data = {}
    for fpath in sorted(glob(f'{results_dir}/*.h5')):
        solver = basename(fpath).replace('.h5', '')
        with HDFArchive(fpath, 'r') as ar:
            data[solver] = {}
            for key in ar:
                data[solver][key] = ar[key]
    return data


# =====================================================================
# Derived quantities
# =====================================================================

def compute_sigma(data, G0_iw):
    """Derive self-energy Sigma(iw) = G0^{-1} - G^{-1} for all solvers that have G.

    Parameters
    ----------
    data : dict
        Output of load_all_results().
    G0_iw : BlockGf
        Non-interacting Green function.

    Returns
    -------
    dict : {solver_name: BlockGf Sigma_iw}
    """
    sigma = {}
    for solver, obs in data.items():
        if 'G' not in obs:
            continue
        G = obs['G']
        S = G0_iw.copy()
        S << inverse(G0_iw) - inverse(G)
        sigma[solver] = S
    return sigma


# =====================================================================
# Deviation tables
# =====================================================================

def _max_abs_diff(a, b, block_lst=None):
    """Max-norm |a - b| for BlockGf (iterate blocks) or single Gf (direct .data).

    block_lst restricts BlockGf iteration; for Gf it is ignored.
    """
    if isinstance(a, BlockGf):
        blocks = block_lst if block_lst is not None else [bl for bl, _ in a]
        dev = 0.0
        for bl in blocks:
            dev = max(dev, float(np.max(np.abs(a[bl].data - b[bl].data))))
        return dev
    return float(np.max(np.abs(a.data - b.data)))


def compute_deviations(obs_dict, block_lst=None):
    """Compute pairwise max-norm deviations for an observable.

    Returns
    -------
    dict : {(solver_a, solver_b): float}
        Upper-triangle pairwise max-norm deviations.
    """
    solvers = sorted(obs_dict.keys())
    devs = {}
    for i, s1 in enumerate(solvers):
        for s2 in solvers[i + 1:]:
            devs[(s1, s2)] = _max_abs_diff(obs_dict[s1], obs_dict[s2], block_lst)
    return devs


def deviation_vs_reference(obs_dict, ref_solver, block_lst=None):
    """Max-norm deviation of each solver vs a chosen reference.

    Returns
    -------
    dict : {solver_name: float}
        Empty if `ref_solver` is not in `obs_dict`.
    """
    if ref_solver not in obs_dict:
        return {}
    ref = obs_dict[ref_solver]
    return {
        solver: _max_abs_diff(val, ref, block_lst)
        for solver, val in obs_dict.items()
        if solver != ref_solver
    }


def deviation_table(obs_dict, block_lst=None, label='G', return_markdown=False):
    """Pairwise max-norm deviation table. Prints by default, or returns markdown.

    Parameters
    ----------
    obs_dict : dict
        {solver_name: BlockGf or Gf}
    block_lst : list of str, optional
        Block names to compare. Required for BlockGf; ignored for single Gf.
    label : str
        Label for the table header.
    return_markdown : bool
        If True, return a markdown upper-triangle deviation matrix as a string
        (diagonal marked '—', below-diagonal left blank).
    """
    solvers = sorted(obs_dict.keys())
    n = len(solvers)
    if n < 2:
        msg = f"Need at least 2 solvers for deviation table (have {n})"
        if return_markdown:
            return f"_{msg}._\n"
        print(f"  {msg}")
        return

    if return_markdown:
        return _deviation_table_markdown(obs_dict, solvers, block_lst, label)

    width = max(len(s) for s in solvers) + 2
    header = " " * width + "".join(s.rjust(width) for s in solvers)
    print(f"\n  {label} max-norm deviations:")
    print(f"  {header}")
    for i, s1 in enumerate(solvers):
        row = s1.rjust(width)
        for j, s2 in enumerate(solvers):
            if j <= i:
                row += " " * width
            else:
                d = _max_abs_diff(obs_dict[s1], obs_dict[s2], block_lst)
                row += f"{d:{width}.2e}"
        print(f"  {row}")


def _deviation_table_markdown(obs_dict, solvers, block_lst, label):
    """Upper-triangle deviation matrix rendered as GitHub-flavored markdown."""
    lines = [f"### {label} max-norm deviations", ""]
    lines.append("| | " + " | ".join(solvers) + " |")
    lines.append("|---|" + "---|" * len(solvers))
    for i, s1 in enumerate(solvers):
        cells = []
        for j, s2 in enumerate(solvers):
            if j < i:
                cells.append("")
            elif j == i:
                cells.append("—")
            else:
                d = _max_abs_diff(obs_dict[s1], obs_dict[s2], block_lst)
                cells.append(f"{d:.2e}")
        lines.append(f"| **{s1}** | " + " | ".join(cells) + " |")
    lines.append("")
    return "\n".join(lines)


# =====================================================================
# Plotting helpers
# =====================================================================

MARKERS = ['-x', '-+', '-^', '-v', '-<', '->', '-*', '-p', '-s', '-d',
           '-o', '-1', '-2', '-3', '-4']


def plot_iw_comparison(obs_dict, name, block_lst, n_iw_plot=None, component=(0, 0)):
    """Plot Matsubara-frequency observable comparison across solvers.

    Parameters
    ----------
    obs_dict : dict
        {solver_name: BlockGf}
    name : str
        Observable name for labels (e.g. 'G', 'Sigma').
    block_lst : list of str
        Block names to plot.
    n_iw_plot : int, optional
        Number of positive Matsubara frequencies to plot. None = all.
    component : tuple
        Orbital indices (i, j) to plot.
    """
    from triqs.plot.mpl_interface import oplot, plt

    solvers = sorted(obs_dict.keys())
    n_blocks = len(block_lst)
    fig, axes = plt.subplots(n_blocks, 1, figsize=(10, 6 * n_blocks), squeeze=False)

    i_orb, j_orb = component
    for idx, block in enumerate(block_lst):
        ax = axes[idx, 0]
        ax.set_title(f"{name}[{block}]")
        for k, solver in enumerate(solvers):
            marker = MARKERS[k % len(MARKERS)]
            g = obs_dict[solver][block][i_orb, j_orb]
            if n_iw_plot is not None:
                oplot(g, marker, name=f"{solver}", x_window=(0, n_iw_plot), axes=ax)
            else:
                oplot(g, marker, name=f"{solver}", axes=ax)
        ax.set_xlabel(r"$\omega_n$")
        ax.set_ylabel(f"{name}[{block}]" + r"$(i\omega_n)$")
        ax.legend(fontsize='small')

    plt.tight_layout()
    return fig


def plot_iw_scalar_comparison(obs_dict, name):
    """Plot a scalar Matsubara quantity (single Gf, no block structure) across solvers.

    Used for chi2_{d,m,s,t} which are single Gf objects on a bosonic mesh.
    """
    from triqs.plot.mpl_interface import oplot, plt

    solvers = sorted(obs_dict.keys())
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_title(name)
    for k, solver in enumerate(solvers):
        oplot(obs_dict[solver], MARKERS[k % len(MARKERS)], name=solver, axes=ax)
    ax.set_xlabel(r"$\omega_n$")
    ax.set_ylabel(name)
    ax.legend(fontsize='small')
    plt.tight_layout()
    return fig


def plot_w_comparison(obs_dict, name, block_lst, component=(0, 0), spectral=False):
    """Plot real-frequency observable comparison across solvers.

    Parameters
    ----------
    obs_dict : dict
        {solver_name: BlockGf with MeshReFreq}
    name : str
        Observable name for labels.
    block_lst : list of str
        Block names to plot.
    component : tuple
        Orbital indices.
    spectral : bool
        If True, plot A(w) = -1/pi * Im G^R(w).
    """
    from triqs.plot.mpl_interface import oplot, plt

    solvers = sorted(obs_dict.keys())
    n_blocks = len(block_lst)
    fig, axes = plt.subplots(n_blocks, 1, figsize=(10, 6 * n_blocks), squeeze=False)

    i_orb, j_orb = component
    for idx, block in enumerate(block_lst):
        ax = axes[idx, 0]
        ylabel = f"A[{block}](w)" if spectral else f"{name}[{block}](w)"
        ax.set_title(ylabel)
        for k, solver in enumerate(solvers):
            g = obs_dict[solver][block][i_orb, j_orb]
            w = np.array([float(x) for x in g.mesh])
            if spectral:
                y = -1.0 / np.pi * g.data.imag
            else:
                y = g.data.real
            ax.plot(w, y, MARKERS[k % len(MARKERS)], label=solver, markersize=3)
        ax.set_xlabel(r"$\omega$")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize='small')

    plt.tight_layout()
    return fig


def plot_static_obs_table(data, obs_name='density', return_markdown=False):
    """Static-observable comparison across solvers. Prints by default, or returns markdown.

    Parameters
    ----------
    data : dict
        Output of load_all_results().
    obs_name : str
        Name of static observable ('density' or 'nn_ab').
    return_markdown : bool
        If True, return a markdown table string (solvers as rows, sub-keys as columns).
    """
    solvers = sorted(s for s, d in data.items()
                     if 'static_obs' in d and obs_name in d['static_obs'])
    if not solvers:
        msg = f"No solver has static_obs/{obs_name}"
        if return_markdown:
            return f"_{msg}._\n"
        print(f"  {msg}")
        return

    if return_markdown:
        return _static_obs_markdown(data, solvers, obs_name)

    print(f"\n  {obs_name}:")
    for solver in solvers:
        val = data[solver]['static_obs'][obs_name]
        if hasattr(val, 'items'):
            items = ", ".join(f"{k}: {v:.6f}" for k, v in val.items())
            print(f"    {solver:20s}  {items}")
        else:
            print(f"    {solver:20s}  {val}")


def _static_obs_markdown(data, solvers, obs_name):
    """Render a static observable as a markdown table.

    Handles both scalar and dict-valued observables (e.g. density per orbital,
    nn_ab per orbital pair). Column order is taken from the first solver that
    has sub-keys; missing cells are left blank.
    """
    subkeys = None
    for solver in solvers:
        val = data[solver]['static_obs'][obs_name]
        if hasattr(val, 'items'):
            subkeys = list(val.keys())
            break

    lines = [f"### {obs_name}", ""]
    if subkeys is None:
        lines.append("| solver | value |")
        lines.append("|---|---|")
        for solver in solvers:
            lines.append(f"| {solver} | {data[solver]['static_obs'][obs_name]} |")
    else:
        header_keys = [str(k) for k in subkeys]
        lines.append("| solver | " + " | ".join(header_keys) + " |")
        lines.append("|---|" + "---|" * len(subkeys))
        for solver in solvers:
            val = data[solver]['static_obs'][obs_name]
            cells = []
            for k in subkeys:
                v = val.get(k) if hasattr(val, 'get') else None
                cells.append(f"{v:.6f}" if isinstance(v, (int, float, np.floating)) else "")
            lines.append(f"| {solver} | " + " | ".join(cells) + " |")
    lines.append("")
    return "\n".join(lines)


def plot_chi_contour(chi_dict, channel, omega_idx=0):
    """2D contour plot for chi3 or chi4 comparison.

    Parameters
    ----------
    chi_dict : dict
        {solver_name: Gf or array}
    channel : str
        Physical channel label ('d', 'm', 's', 't').
    omega_idx : int
        Bosonic frequency index for the slice (default Omega=0).
    """
    from triqs.plot.mpl_interface import plt

    solvers = sorted(chi_dict.keys())
    n = len(solvers)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5), squeeze=False)

    for k, solver in enumerate(solvers):
        ax = axes[0, k]
        chi = chi_dict[solver]
        if hasattr(chi, 'data'):
            # data layout: (freq1, freq2, ..., orb1, orb2, ...)
            # Number of frequency dims = number of mesh components
            n_freq = len(chi.mesh.components) if hasattr(chi.mesh, 'components') else 1
            sliced = chi.data[omega_idx]
            # Select first orbital component for all trailing dims
            for _ in range(sliced.ndim - (n_freq - 1)):
                sliced = sliced[..., 0]
            data = sliced.real
        else:
            data = np.array(chi)

        ax.set_title(f"chi_{channel} - {solver}")
        if data.ndim >= 2:
            im = ax.imshow(data, origin='lower', aspect='auto', cmap='RdBu_r')
            plt.colorbar(im, ax=ax)
        else:
            ax.plot(data)
            ax.set_xlabel(r"$\nu_n$")

    plt.tight_layout()
    return fig
