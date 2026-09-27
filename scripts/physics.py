"""SI blackbody model. I_dot counts irreversibly erased unbiased bits/s.

Baseline: eta=1, T_bath=T_rad, zero incident background, isotropic total
luminosity. A is TOTAL emitting area, not projected area. No facility
sensitivity is silently built into this module.
"""
import numpy as np

K_B = 1.380649e-23
H = 6.62607015e-34
C = 299792458.0
SIGMA = 5.670374419e-8
PC_TO_M = 3.085677581491367e16
AU = 149597870700.0
JY_TO_W_M2_HZ = 1e-26
L_SUN = 3.828e26
N_STARS_MW = 1e11
L_TIR_MW = 1e10 * L_SUN  # illustrative normalization, not a measured limit


def positive(value, name):
    value = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(value)) or np.any(value <= 0):
        raise ValueError(f"{name} must be positive and finite")
    return value


def planck_b_nu(nu, T):
    """W m^-2 Hz^-1 sr^-1; preserves the Wien exponential and RJ limit."""
    nu, T = np.broadcast_arrays(positive(nu, 'nu'), positive(T, 'T'))
    x = H * nu / (K_B * T)
    with np.errstate(over='ignore', under='ignore', divide='ignore'):
        log_denominator = np.where(x < 50, np.log(np.expm1(x)),
                                   x + np.log1p(-np.exp(-x)))
        return np.exp(np.log(2 * H / C**2) + 3 * np.log(nu) - log_denominator)


def waste_power(I_dot=1e40, T_bath=10., eta=1.):
    if np.any(np.asarray(eta) < 1):
        raise ValueError('eta must be >= 1 for this Landauer-normalized model')
    return positive(I_dot, 'I_dot') * K_B * positive(T_bath, 'T_bath') * np.log(2) * positive(eta, 'eta')


def radiator_area(T_rad=10., I_dot=1e40, T_bath=None, eta=1., T_bg=0.):
    T_rad = positive(T_rad, 'T_rad')
    if not np.isfinite(T_bg) or T_bg < 0 or np.any(T_rad <= T_bg):
        raise ValueError('require 0 <= T_bg < T_rad')
    T_bath = T_rad if T_bath is None else T_bath
    return waste_power(I_dot, T_bath, eta) / (SIGMA * (T_rad**4 - T_bg**4))


def flux_ujy(wavelength_um, T_rad=10., I_dot=1e40, d_pc=10.,
             T_bath=None, eta=1., T_bg=0.):
    """Unresolved isotropic flux; if T_bg>0, returns background contrast."""
    nu = C / (positive(wavelength_um, 'wavelength_um') * 1e-6)
    A = radiator_area(T_rad, I_dot, T_bath, eta, T_bg)
    B = planck_b_nu(nu, T_rad)
    if T_bg:
        B = B - planck_b_nu(nu, T_bg)
    return A * B / (4 * (positive(d_pc, 'd_pc') * PC_TO_M)**2) / 1e-32


def spectral_index(T, short_um=450., long_um=850., beta=0.):
    if short_um >= long_um:
        raise ValueError('require short_um < long_um')
    nu1, nu2 = C / (positive(np.array([short_um, long_um]), 'wavelength') * 1e-6)
    return float(np.log(planck_b_nu(nu1, T) / planck_b_nu(nu2, T)) / np.log(nu1 / nu2) + beta)


def detection_horizon(wavelength_um, T_rad, rms_ujy, snr=5., **kwargs):
    """Conditional distance in pc, given TOTAL effective 1-sigma noise."""
    return 10 * np.sqrt(flux_ujy(wavelength_um, T_rad, d_pc=10., **kwargs)
                        / (positive(rms_ujy, 'rms') * positive(snr, 'snr')))


def illustrative_fraction_limit(I_dot, T=15., epsilon=1e-3):
    """Unclipped conditional luminosity budget. Values >1 mean no bound."""
    return positive(epsilon, 'epsilon') * L_TIR_MW / (N_STARS_MW * waste_power(I_dot, T))
