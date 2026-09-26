import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8.5,
    "mathtext.fontset": "cm",
    "pdf.fonttype": 42,
})

K_B = 1.380649e-23
H = 6.62607015e-34
C = 2.99792458e8
SIGMA = 5.670374419e-8
PC_TO_M = 3.085677581e16
JY_TO_W_M2_HZ = 1e-26
L_SUN = 3.828e26
N_STARS_MW = 1.0e11
L_TIR_MW = 1.0e10 * L_SUN

def planck_b_nu(nu, T):
    x = (H * nu) / (K_B * T)
    x = np.clip(x, 1e-10, 700.0)
    return (2.0 * H * (nu**3) / (C**2)) / (np.expm1(x))

def make_figure_spectra():
    wavelengths_um = np.logspace(0.5, 3.5, 600)
    nu = C / (wavelengths_um * 1e-6)
    d_m = 10.0 * PC_TO_M

    def get_flux(T):
        P_waste = 1e40 * K_B * T * np.log(2)
        A_rad = P_waste / (SIGMA * (T**4))
        B_nu = planck_b_nu(nu, T)
        return ((A_rad / (4.0 * d_m**2)) * B_nu / JY_TO_W_M2_HZ) * 1e6

    flux_10k = get_flux(10.0)
    flux_30k = get_flux(30.0)
    raw_dust = (nu**1.7) * planck_b_nu(nu, 15.0)
    flux_dust = (raw_dust / np.max(raw_dust)) * np.max(flux_10k)

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    ax.axvspan(5.0, 28.0, color="#f0883e", alpha=0.08, label="JWST/MIRI")
    ax.axvspan(25.0, 250.0, color="#3fb950", alpha=0.08, label="Far-IR Gap (PRIMA)")
    ax.axvspan(350.0, 3000.0, color="#d29922", alpha=0.08, label="ALMA Band")

    ax.plot(wavelengths_um, flux_10k, color="#1f77b4", lw=2.0, label=r"Inward Node ($10\ \mathrm{K},\ \beta=0$)")
    ax.plot(wavelengths_um, flux_30k, color="#9467bd", lw=1.8, ls="--", label=r"Inward Node ($30\ \mathrm{K},\ \beta=0$)")
    ax.plot(wavelengths_um, flux_dust, color="#2ca02c", lw=1.8, ls="-.", label=r"Cold Dust ($15\ \mathrm{K},\ \beta=1.7$)")

    ax.plot([5.6, 7.7, 10.0, 15.0, 21.0, 25.5], [1.5, 2.5, 4.0, 10.0, 30.0, 70.0], "o-", color="#e6550d", lw=1.8, markersize=4, label=r"JWST/MIRI (10k-sec, $10\sigma$)")
    ax.plot([70.0, 100.0, 160.0, 250.0, 350.0, 500.0], [4000.0, 4500.0, 7000.0, 6000.0, 7000.0, 8000.0], "x--", color="#d62728", lw=1.5, markersize=5, label="Herschel (PACS/SPIRE Archive)")
    ax.plot([25.0, 40.0, 70.0, 100.0, 150.0, 200.0], [25.0, 15.0, 10.0, 12.0, 18.0, 25.0], "*-", color="#2ca02c", lw=1.8, markersize=6, label="PRIMA (Future Probe Goal)")
    ax.plot([450.0, 850.0, 1300.0, 3000.0], [1000.0, 200.0, 100.0, 50.0], "s-", color="#bc6b00", lw=1.8, markersize=4, label=r"ALMA (1-hr, $10\sigma$)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(3.0, 3200.0)
    ax.set_ylim(1e-1, 5e5)
    ax.set_xlabel(r"Wavelength $\lambda\ [\mu\mathrm{m}]$")
    ax.set_ylabel(r"Flux Density $S_\nu\ [\mu\mathrm{Jy}]\ (d = 10\ \mathrm{pc})$")
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(loc="upper right", frameon=True, framealpha=0.9)
    plt.tight_layout()
    plt.savefig("figures/fig2_spectra.pdf", bbox_inches="tight")
    print("Saved: figures/fig2_spectra.pdf")

def make_figure_exclusion():
    log_I_dot = np.linspace(25, 45, 400)
    I_dot = 10.0**log_I_dot
    P_waste = I_dot * K_B * 15.0 * np.log(2)
    f_cobe = np.clip((0.001 * L_TIR_MW) / (N_STARS_MW * P_waste), 1e-6, 1.0)

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.fill_between(log_I_dot, f_cobe, 1.0, color="#d62728", alpha=0.2, label="Excluded by COBE/FIRAS (Diffuse IR)")
    ax.plot(log_I_dot, f_cobe, color="#d62728", lw=2.0)
    ax.axvline(x=np.log10(2.5e38), color="#2ca02c", lw=1.8, ls="--", label=r"PRIMA Discovery Limit ($d \leq 10\ \mathrm{pc}$)")
    ax.axvline(x=np.log10(1.2e40), color="#bc6b00", lw=1.8, ls=":", label=r"ALMA Sensitivity Limit ($d \leq 10\ \mathrm{pc}$)")
    ax.axhline(y=0.14, color="#1f77b4", lw=2.0, ls="-.", label=r"Steady-State Occupancy ($f_{\mathrm{intro}} \approx 14\%$)")
    ax.scatter([40], [0.14], color="#1f77b4", s=90, edgecolors="black", zorder=5)
    ax.text(37.0, 0.18, r"Simulated Node ($\dot{I}=10^{40},\ f=14\%$)", color="#1f77b4", weight="bold")
    ax.scatter([25], [1.0 / N_STARS_MW], color="#2ca02c", s=70, edgecolors="black", zorder=5)
    ax.text(25.3, 1.5e-6, "Earth (2026)", color="#2ca02c")

    ax.set_yscale("log")
    ax.set_xlim(25, 45)
    ax.set_ylim(1e-6, 1.0)
    ax.set_xlabel(r"Computational Capacity per Node $\log_{10}(\dot{I})\ [\mathrm{bits\ s^{-1}}]$")
    ax.set_ylabel(r"Galactic Inward Occupancy Fraction $f_{\mathrm{civ}}$")
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(loc="lower left", frameon=True, framealpha=0.9)
    plt.tight_layout()
    plt.savefig("figures/fig3_exclusion.pdf", bbox_inches="tight")
    print("Saved: figures/fig3_exclusion.pdf")

if __name__ == "__main__":
    make_figure_spectra()
    make_figure_exclusion()
