import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from triqs.utility import mpi

# print on master node
def mpi_print(arg):
    if mpi.is_master_node():
        print(arg)
