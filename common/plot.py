import sys, os
sys.path.append(os.getcwd() + "/..")
sys.path.append(os.getcwd() + "/../../common")
from model import *

from h5 import HDFArchive
from triqs.plot.mpl_interface import oplot, plt
from glob import glob
from os.path import basename
import numpy as np

# === Load chi3pp and chi3ph for every solver

solver_lst = [ basename(f).strip('.h5') for f in glob('results/*.h5') ]
chi3pp_iw, chi3ph_iw = {}, {}

for solver in solver_lst:
    dat = HDFArchive('results/' + solver + '.h5','r')
    chi3pp_iw[solver] = dat['chi3pp_iw']
    chi3ph_iw[solver] = dat['chi3ph_iw']
    
# === For every block and solver, plot chi3pp and chi3ph

block_lst = [bl for bl in chi3pp_iw[solver_lst[0]].indices]
n_blocks = len(block_lst)
n_solvers = len(solver_lst)

for chi3, name in [[chi3pp_iw, r'$\chi_3^{\rm PP}$'], [chi3ph_iw, r'$\chi_3^{\rm PH}$']]:

    for block in block_lst:
        
        vmin = min([np.min(chi3[solver][block].data.real) for solver in solver_lst])
        vmax = max([np.max(chi3[solver][block].data.real) for solver in solver_lst])
        
        if vmin == vmax: continue
            
        plt.subplots(n_solvers, figsize=(3*n_solvers + 1,3))
        plt.suptitle("Re " + name + "[" + str(block) + r"][0,0,0,0]", fontsize=14, y=1.1)
        count = 1
        
        for solver in solver_lst:
            fig = plt.subplot(1, n_solvers, count); 
            plt.pcolormesh(chi3[solver][block].data[:,:,0,0,0,0].real, vmin=vmin, vmax=vmax)
            fig.set_title(solver)
            fig.axes.set_aspect('equal')
            plt.xlabel(r"$\Omega$")
            plt.ylabel(r"$\omega$")
            count += 1
            
        plt.colorbar() 
        plt.tight_layout()

    plt.show()
