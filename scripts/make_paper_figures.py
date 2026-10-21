"""Reproducible spectra, numerical tables and CONDITIONAL budget diagram.

No hard-coded PRIMA/ALMA detection promises. Figure 3 uses a declared,
illustrative total RMS of 10 microJy and an assumed luminosity budget.
Neither is an observed exclusion. Data accompany each figure.
"""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from physics import (flux_ujy, waste_power, radiator_area, spectral_index,
                     detection_horizon, illustrative_fraction_limit, AU,
                     L_TIR_MW, N_STARS_MW)

ROOT = Path(__file__).resolve().parents[1]


def save(fig, name):
    for extension in ('png', 'pdf'):
        fig.savefig(ROOT / 'figures' / f'{name}.{extension}', dpi=180, bbox_inches='tight')
    plt.close(fig)


def main():
    for folder in ('figures', 'results'):
        (ROOT / folder).mkdir(exist_ok=True)
    waves = [28, 92, 100, 126, 183, 235, 450, 850, 1300, 3000]
    with (ROOT / 'results' / 'flux_table.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['wavelength_um', '10K_ujy', '15K_ujy', '30K_ujy'])
        for w in waves:
            writer.writerow([w] + [float(flux_ujy(w, t)) for t in (10, 15, 30)])
    numbers = {
        'baseline': {'I_dot_bits_per_s': 1e40, 'd_pc': 10, 'eta': 1, 'T_bg_K': 0},
        'temperatures': {str(t): dict(power_W=float(waste_power(T_bath=t)),
                                      area_m2=float(radiator_area(t)),
                                      disk_one_face_radius_AU=float(np.sqrt(radiator_area(t)/np.pi)/AU),
                                      spherical_radius_AU=float(np.sqrt(radiator_area(t)/(4*np.pi))/AU),
                                      alpha_450_850=spectral_index(t),
                                      peak_Snu_wavelength_um=5099.4428295/t)
                         for t in (10, 15, 30)},
        'dust_15K_beta1.7_alpha_450_850': spectral_index(15, beta=1.7),
        'conditional_budget': {'epsilon': .001, 'f': .14, 'T_K': 15,
                               'L_excess_W': float(.14*N_STARS_MW*waste_power(T_bath=15)),
                               'I_dot_max': float(.001*L_TIR_MW/(.14*N_STARS_MW*waste_power(1.,15)))},
        'equal_SNR_for_sigma_alpha_0.15': float(np.sqrt(2)/(.15*np.log(850/450))),
        'conditional_dmax_pc_rms10_snr5': {f'{t}K_{w}um': float(detection_horizon(w,t,10))
                                         for t in (10,15,30) for w in (100,183,235,450,850)}
    }
    (ROOT / 'results' / 'numerical_summary.json').write_text(json.dumps(numbers, indent=2)+'\n')
    wavelengths = np.geomspace(20, 3500, 800)
    fig, ax = plt.subplots(figsize=(8,5), layout='constrained')
    for t, color in [(10,'#2874a6'), (15,'#218c74'), (30,'#8e44ad')]:
        ax.loglog(wavelengths, flux_ujy(wavelengths,t), label=f'{t} K blackbody', color=color)
    # Normalize the illustrative dust template at a stated wavelength, not its sampled maximum.
    dust = flux_ujy(wavelengths,15) * (850/wavelengths)**1.7
    dust *= flux_ujy(850,10)/flux_ujy(850,15)
    ax.loglog(wavelengths,dust,ls='--',color='#777777',label='15 K dust, beta=1.7; matched at 850 um')
    ax.scatter([92,126,183,235], flux_ujy([92,126,183,235],10), color='#2874a6',s=25,label='PRIMAger nominal band centers (2025)')
    ax.scatter([250,350,500], [5800,6300,6800], marker='x',color='#c0392b',label='SPIRE 1-sigma confusion (Nguyen et al. 2010)')
    ax.set(xlim=(20,3500),ylim=(.01,2e4),xlabel='Wavelength (micrometers)',ylabel='Flux density (microJy)',title='Model fluxes: 10^40 erased bits/s, 10 pc, eta=1')
    ax.legend(fontsize=8,loc='upper left',bbox_to_anchor=(1.02,1))
    ax.grid(alpha=.2,which='both')
    save(fig, 'fig2_spectra')
    rates = np.geomspace(1e38,1e47,500)
    limit = illustrative_fraction_limit(rates)
    fig, ax = plt.subplots(figsize=(8,5),layout='constrained')
    ax.loglog(rates,limit,color='#c0392b',label='Assumed budget: 0.001 x 10^10 solar luminosities')
    ax.fill_between(rates, limit, 1, where=limit<1, color='#c0392b',alpha=.1,label='Above assumed budget; NOT an empirical exclusion')
    for w,color in [(235,'#218c74'),(850,'#8e44ad')]:
        threshold = 1e40 * 5 * 10 / float(flux_ujy(w,15))
        ax.axvline(threshold,color=color,ls='--',label=f'{w} um: conditional 5-sigma, RMS=10 microJy')
    ax.scatter([1e40],[.14],color='#2874a6',label='Illustrative f=0.14 (not calibrated to stars)')
    ax.set(xlim=(1e38,1e47),ylim=(1e-4,1),xlabel='Irreversible erasure rate (bits/s)',ylabel='Assumed occupied stellar fraction',title='Conditional luminosity budget and detection thresholds: 15 K')
    ax.legend(fontsize=8,loc='lower left')
    ax.grid(alpha=.2,which='both')
    save(fig,'fig3_conditional_budget')
    print(json.dumps(numbers,indent=2))


if __name__ == '__main__':
    main()
