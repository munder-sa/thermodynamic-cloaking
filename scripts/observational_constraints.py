"""Conditional benchmarks, NOT an end-to-end survey completeness calculation."""
import json
from pathlib import Path
import numpy as np
from physics import flux_ujy, waste_power, detection_horizon, L_SUN
ROOT=Path(__file__).resolve().parents[1]

def main():
    budget=3.828e33
    f,N=.14,1e11
    per_node=budget/(f*N)
    temperatures={str(t):dict(energy_per_bit_J=float(waste_power(1,t)),
        per_node_erasure_limit_bits_s=float(per_node/waste_power(1,t)),
        dex_above_benchmark=float(np.log10(per_node/waste_power(1e40,t)))) for t in (15,20)}
    d=float(detection_horizon(235,10,46000/5))
    result=dict(prima_235=dict(classical_5sigma_confusion_ujy=46000,
        model_flux_at_10pc_ujy=float(flux_ujy(235,10)),conditional_distance_pc=d,
        conditional_distance_AU=d*206264.806247,
        beam_fwhm_arcsec_2025_design=27.9,
        source='https://academic.oup.com/mnras/article/532/2/1966/7696741'),
        budget=dict(galactic_W=budget,assumed_f=f,assumed_N_stars=N,
                    per_node_W=per_node,per_node_solar_luminosities=per_node/L_SUN,
                    temperatures=temperatures),
        parallax_arcsec={str(au):206264.806247/au for au in (1000,10000,52600)})
    (ROOT/'results'/'observational_constraints.json').write_text(json.dumps(result,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,4.5),layout='constrained')
    distances=np.geomspace(.05,20,300)
    ax.loglog(distances,flux_ujy(235,10,d_pc=distances)/1000,label='10 K benchmark at 235 micrometers')
    ax.axhline(46,color='#c0392b',ls='--',label='Classical confusion: 5 sigma = 46 mJy')
    ax.axhline(.05,color='#218c74',ls=':',label='Hypothetical 5 sigma = 50 microJy')
    ax.scatter([10],[float(flux_ujy(235,10))/1000],color='black',label='10 pc: 0.0298 mJy')
    ax.set(xlabel='Distance (pc)',ylabel='Monochromatic flux density (mJy)',
           title='Source model versus two distinct threshold assumptions')
    ax.grid(alpha=.2,which='both');ax.legend(fontsize=8)
    for ext in ('pdf','png'):fig.savefig(ROOT/'figures'/f'fig5_confusion_comparison.{ext}',dpi=180)
    plt.close(fig)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
