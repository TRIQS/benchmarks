# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Sr2RuO4
#
# Three-band effective model for the layered perovskite Sr$_2$RuO$_4$.

# %%
import sys
sys.path.append('../common')
from model import *
from analysis import load_all_results, compute_sigma, deviation_table, \
    plot_iw_comparison, plot_w_comparison, plot_static_obs_table, plot_chi_contour

# %%
data = load_all_results()
block_lst = list(data[next(iter(data))]['G'].indices) if data else []

# %% [markdown]
# ## 1. Single-Particle Green Function $G(i\omega_n)$

# %%
G_dict = {s: d['G'] for s, d in data.items() if 'G' in d}
if G_dict:
    plot_iw_comparison(G_dict, 'G', block_lst)
    deviation_table(G_dict, block_lst, label='G')

# %% [markdown]
# ## 2. Self-Energy $\Sigma(i\omega_n)$

# %%
Sigma_dict = compute_sigma(data, G0_iw)
if Sigma_dict:
    plot_iw_comparison(Sigma_dict, r'$\Sigma$', block_lst)
    deviation_table(Sigma_dict, block_lst, label='Sigma')

# %% [markdown]
# ## 3. Static Observables

# %%
plot_static_obs_table(data, 'density')
plot_static_obs_table(data, 'nn_ab')

# %% [markdown]
# ## 4. Real-Frequency Spectral Function $A(\omega)$

# %%
G_w_dict = {s: d['G_w'] for s, d in data.items() if 'G_w' in d}
if G_w_dict:
    plot_w_comparison(G_w_dict, 'G_w', block_lst, spectral=True)
else:
    print("No solver has real-frequency data.")

# %% [markdown]
# ## 5. Three-Point Correlator $\chi_3$

# %%
for ch in ['d', 'm', 's', 't']:
    key = f'chi3_{ch}'
    chi3 = {s: d[key] for s, d in data.items() if key in d}
    if chi3:
        plot_chi_contour(chi3, ch)

# %% [markdown]
# ## 6. Two-Particle Green Function $\chi_4$ / $G_{2c}$

# %%
for ch in ['d', 'm', 's', 't']:
    key = f'chi4_{ch}'
    chi4 = {s: d[key] for s, d in data.items() if key in d}
    if chi4:
        plot_chi_contour(chi4, ch)
