import sys, os
sys.path.append(os.getcwd() + "/..")
sys.path.append(os.getcwd() + "/../../common")
from model import *

from pytriqs.archive import HDFArchive
from pytriqs.plot.mpl_interface import oplot, plt
from glob import glob
from os.path import basename

# === Load chi2pp and chi2ph for every solver

solver_lst = [ basename(f).strip('.h5') for f in glob('results/*.h5') ]
chi2pp_tau, chi2ph_tau = {}, {}

for solver in solver_lst:
    dat = HDFArchive('results/' + solver + '.h5','r')
    chi2pp_tau[solver] = dat['chi2pp_tau']
    chi2ph_tau[solver] = dat['chi2ph_tau']
    
# === For every block and solver, plot chi2pp and chi2ph

block_lst = [bl for bl in chi2pp_tau[solver_lst[0]].indices]
n_blocks = len(block_lst)

for chi2, name in [[chi2pp_tau, r'$\chi_2^{\rm PP}$'], [chi2ph_tau, r'$\chi_2^{\rm PH}$']]:

    for block in block_lst:

        if max(norm_inf(chi2[solver][block][0,0,0,0]) for solver in solver_lst) < 1e-14: continue
        
        plt.figure(figsize=(10,6))
        for solver in solver_lst:
            oplot(chi2[solver][block][0,0,0,0], name = name + "_%s" % solver)

        plt.suptitle(name + "[" + str(block) + r"][0,0,0,0]", fontsize=14, y=1.05)
        plt.xlabel(r"$\tau$")

        plt.tight_layout()
        plt.show()
