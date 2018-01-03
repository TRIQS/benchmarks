execfile('../system.py')

from pytriqs.archive import HDFArchive
from math import exp

# Exact expressions for parition function, occupation numbers and Green function
# cf. Eq. (B.3) and (B.4) of S. Pairault, et al., EPJ B 16, 85-105 (2000)

Z = 1. + exp((mu+h)*beta) + exp((mu-h)*beta) + exp((2.*mu-U)*beta)
n = {}

G_iw = BlockGf_from_struct(mesh=iw_mesh, struct=gf_struct)

n['up'] = ( exp((mu+h)*beta) + exp((2.*mu-U)*beta) ) / Z
n['dn'] = ( exp((mu-h)*beta) + exp((2.*mu-U)*beta) ) / Z

G_iw['up'] = Gf(mesh=iw_mesh, target_shape=(1,1))
G_iw['dn'] = Gf(mesh=iw_mesh, target_shape=(1,1))

G_iw['up'] << ( 1 - n['dn'] ) * inverse( iOmega_n + mu + h ) + \
                    n['dn']   * inverse( iOmega_n + mu + h  - U )
G_iw['dn'] << ( 1 - n['up'] ) * inverse( iOmega_n + mu - h ) + \
                    n['up']   * inverse( iOmega_n + mu - h  - U )


# -------- Save in archive ---------
with HDFArchive("../results/exact.h5",'w') as res:
    res["G"] = G_iw
