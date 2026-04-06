"""Shared analysis and plotting library for solver benchmark notebooks.

Extends the functionality of the original plot.py with support for
multi-observable comparison, static observables, real-frequency data,
and two-particle correlators.
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
        ar = HDFArchive(fpath, 'r')
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

def deviation_table(obs_dict, block_lst, label='G'):
    """Print pairwise max-norm deviation table for an observable.

    Parameters
    ----------
    obs_dict : dict
        {solver_name: BlockGf or Gf}
    block_lst : list of str
        Block names to compare.
    label : str
        Label for the table header.
    """
    solvers = sorted(obs_dict.keys())
    n = len(solvers)
    if n < 2:
        print(f"  Need at least 2 solvers for deviation table (have {n})")
        return

    # Compute max deviation across all blocks
    def max_dev(a, b):
        dev = 0.0
        for bl in block_lst:
            diff = a[bl].data - b[bl].data
            dev = max(dev, np.max(np.abs(diff)))
        return dev

    # Header
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
                d = max_dev(obs_dict[s1], obs_dict[s2])
                row += f"{d:{width}.2e}"
        print(f"  {row}")


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
                y = -1.0 / np.pi * g.data[:, 0, 0].imag
            else:
                y = g.data[:, 0, 0].real
            ax.plot(w, y, MARKERS[k % len(MARKERS)], label=solver, markersize=3)
        ax.set_xlabel(r"$\omega$")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize='small')

    plt.tight_layout()
    return fig


def plot_static_obs_table(data, obs_name='density'):
    """Print static observable comparison table across solvers.

    Parameters
    ----------
    data : dict
        Output of load_all_results().
    obs_name : str
        Name of static observable ('density' or 'nn_ab').
    """
    solvers = sorted(s for s, d in data.items() if 'static_obs' in d and obs_name in d['static_obs'])
    if not solvers:
        print(f"  No solver has static_obs/{obs_name}")
        return

    print(f"\n  {obs_name}:")
    for solver in solvers:
        val = data[solver]['static_obs'][obs_name]
        if hasattr(val, 'items'):
            items = ", ".join(f"{k}: {v:.6f}" for k, v in val.items())
            print(f"    {solver:20s}  {items}")
        else:
            print(f"    {solver:20s}  {val}")


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
            # Assume 3-index object chi(Omega, nu, nu') -- take Omega slice
            data_2d = chi.data[omega_idx, :, :].real
        else:
            data_2d = np.array(chi)
        im = ax.imshow(data_2d, origin='lower', aspect='auto', cmap='RdBu_r')
        ax.set_title(f"chi_{channel} - {solver}")
        plt.colorbar(im, ax=ax)

    plt.tight_layout()
    return fig
