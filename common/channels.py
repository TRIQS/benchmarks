"""Channel translation utilities for two-particle quantities.

Conventions follow Appendix C of Wentzell et al., PRB 93, 2016 (arXiv:1610.06520).

Spin decomposition to physical channels:
  Particle-hole:      chi_d = chi_uu + chi_ud,  chi_m = chi_uu - chi_ud
  Particle-particle:  chi_s = chi_ud - chi_udc,  chi_t = chi_ud + chi_udc

These apply uniformly to chi2, chi3, and chi4.
"""


# =====================================================================
# Spin decomposition to physical channels
# =====================================================================

def to_physical_channels_ph(chi_uu, chi_ud):
    """Convert ph-channel spin components to density (d) and magnetic (m).

    chi_d = chi_uu + chi_ud
    chi_m = chi_uu - chi_ud
    """
    chi_d = chi_uu.copy()
    chi_d << chi_uu + chi_ud
    chi_m = chi_uu.copy()
    chi_m << chi_uu - chi_ud
    return chi_d, chi_m


def to_physical_channels_pp(chi_ud, chi_udc):
    """Convert pp-channel spin components to singlet (s) and triplet (t).

    chi_s = chi_ud - chi_udc  (antisymmetric, S=0)
    chi_t = chi_ud + chi_udc  (symmetric, S=1)
    """
    chi_s = chi_ud.copy()
    chi_s << chi_ud - chi_udc
    chi_t = chi_ud.copy()
    chi_t << chi_ud + chi_udc
    return chi_s, chi_t


# =====================================================================
# Frequency parameterization translations for chi3
# =====================================================================
#
# Fermionic convention: chi3(nu1, nu2) with both fermionic
#
# Channel parameterizations (nu1, nu2) <-> (Omega, nu):
#   ph:  nu1=nu, nu2=Omega+nu    =>  Omega=nu2-nu1 (bosonic), nu=nu1
#   pp:  nu1=nu, nu2=Omega-nu    =>  Omega=nu1+nu2 (bosonic), nu=nu1


def chi3_ph_to_ferm(Omega, nu):
    """(Omega, nu) in ph channel -> fermionic (nu1, nu2).

    nu1 = nu,  nu2 = Omega + nu
    """
    return nu, Omega + nu


def chi3_ferm_to_ph(nu1, nu2):
    """Fermionic (nu1, nu2) -> (Omega, nu) in ph channel.

    Omega = nu2 - nu1 (bosonic),  nu = nu1
    """
    return nu2 - nu1, nu1


def chi3_pp_to_ferm(Omega, nu):
    """(Omega, nu) in pp channel -> fermionic (nu1, nu2).

    nu1 = nu,  nu2 = Omega - nu
    """
    return nu, Omega - nu


def chi3_ferm_to_pp(nu1, nu2):
    """Fermionic (nu1, nu2) -> (Omega, nu) in pp channel.

    Omega = nu1 + nu2 (bosonic),  nu = nu1
    """
    return nu1 + nu2, nu1


# =====================================================================
# Frequency parameterization translations for chi4
# =====================================================================
#
# Fermionic convention: chi4(w1, w2, w3) with w4 = w1 - w2 + w3
#
# Channel parameterizations (w1, w2, w3, w4) <-> (nu, nu', Omega):
#   pp:  w1=nu, w2=Omega-nu', w3=Omega-nu, w4=nu'
#   ph:  w1=nu, w2=nu+Omega, w3=nu'+Omega, w4=nu'
#   xph: w1=nu, w2=nu', w3=Omega+nu', w4=nu+Omega


def chi4_ferm_to_pp(chi4_ferm):
    """Translate chi4(w1, w2, w3) -> chi4^pp(nu, nu', Omega).

    nu=w1, nu'=w4=w1-w2+w3, Omega=w1+w3.
    Inverse: w1=nu, w2=Omega-nu', w3=Omega-nu.
    """
    raise NotImplementedError("Requires mesh index remapping -- implement when needed")


def chi4_pp_to_ferm(chi4_pp):
    """Translate chi4^pp(nu, nu', Omega) -> chi4(w1, w2, w3)."""
    raise NotImplementedError("Requires mesh index remapping -- implement when needed")


def chi4_ferm_to_ph(chi4_ferm):
    """Translate chi4(w1, w2, w3) -> chi4^ph(nu, nu', Omega).

    nu=w1, nu'=w4=w1-w2+w3, Omega=w2-w1.
    Inverse: w1=nu, w2=nu+Omega, w3=nu'+Omega.
    """
    raise NotImplementedError("Requires mesh index remapping -- implement when needed")


def chi4_ph_to_ferm(chi4_ph):
    """Translate chi4^ph(nu, nu', Omega) -> chi4(w1, w2, w3)."""
    raise NotImplementedError("Requires mesh index remapping -- implement when needed")


def chi4_ferm_to_xph(chi4_ferm):
    """Translate chi4(w1, w2, w3) -> chi4^xph(nu, nu', Omega).

    Inverse: w1=nu, w2=nu', w3=Omega+nu'.
    """
    raise NotImplementedError("Requires mesh index remapping -- implement when needed")


def chi4_xph_to_ferm(chi4_xph):
    """Translate chi4^xph(nu, nu', Omega) -> chi4(w1, w2, w3)."""
    raise NotImplementedError("Requires mesh index remapping -- implement when needed")


# =====================================================================
# Inter-channel translations (go through fermionic as intermediate)
# =====================================================================

def chi4_pp_to_ph(chi4_pp):
    """Translate chi4 from pp to ph frequency parameterization."""
    return chi4_ferm_to_ph(chi4_pp_to_ferm(chi4_pp))


def chi4_ph_to_pp(chi4_ph):
    """Translate chi4 from ph to pp frequency parameterization."""
    return chi4_ferm_to_pp(chi4_ph_to_ferm(chi4_ph))


def chi4_ph_to_xph(chi4_ph):
    """Translate chi4 from ph to xph frequency parameterization."""
    return chi4_ferm_to_xph(chi4_ph_to_ferm(chi4_ph))
