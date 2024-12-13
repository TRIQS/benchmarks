import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from nrgljubljana_interface import Flat, MeshReFreqPts
import numpy as np

# ==== System Parameters ====
model = 'SIAM'
symtype = 'ISO'

# N. B.: Kondo temperature for these parameters is T_K = 1.35e-8  
TK = 1.3533e-8
T = 1.3533e-11
D = 1.0
U = 2e-3
V = 0.04 * U

mesh_max = 2.0
mesh_min = 1e-12
mesh_ratio = 1.01

mesh_array = log_mesh(mesh_max, mesh_min, mesh_ratio)

mesh = MeshReFreqPts(mesh_array)

gf = Gf(mesh=mesh, target_shape=(1, 1))
Delta_w = BlockGf(name_list=['imp'], block_list=[gf])

Delta_w['imp'] << (2 * V / np.pi) * Flat(D)