from pytriqs.archive import HDFArchive
from pytriqs.gf import Gf
from pytriqs.plot.mpl_interface import oplot, plt
from glob import glob
from os.path import basename

# === Load Green function for every solver

file_lst = glob('results/*.h5')
solver_lst, G = [], {}

for f in file_lst:
    solver = basename(f).strip('.h5')
    solver_lst.append(solver)
    G[solver] = HDFArchive(f,'r')['G']

# === For every block, plot Green functions of all solvers

block_lst = list(G[solver_lst[0]].indices)
n_blocks = len(block_lst)

plt.subplots(n_blocks,1,figsize=(10,6*n_blocks))

for i, block in enumerate(block_lst,1):
    fig = plt.subplot(n_blocks,1,i)
    fig.set_title("G[" + block + "]")
    for solver in solver_lst:
        oplot(G[solver][block], name = "G_%s" % solver)
    plt.xlabel("$\omega_n$")
    plt.ylabel("$G[" + block + "](i\omega_n)$")

plt.tight_layout()
plt.savefig("plot.png", dpi = 200)
plt.show()
