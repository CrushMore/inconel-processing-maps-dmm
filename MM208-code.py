import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import UnivariateSpline, RectBivariateSpline

# ─────────────────────────────────────────────────────────────────────────────
# IN718 Processing Map  —  Prasad DMM approach
# ─────────────────────────────────────────────────────────────────────────────

temperatures = np.array([900, 950, 975, 1000, 1025, 1050, 1100])
strain_rates  = np.array([0.001, 0.01, 0.1, 1.0])
log_edot      = np.log10(strain_rates)

stress_data = {
    0.1: np.array([
        [279.9,158.2,124.7,103.5, 83.0, 61.2, 46.9],
        [359.4,227.0,184.0,151.3,134.2,122.1,116.0],
        [496.2,324.3,288.2,244.0,217.2,190.7,168.6],
        [599.8,466.1,418.1,394.5,333.9,289.9,230.7],
    ]),
    0.2: np.array([
        [299.2,165.5,129.4,109.8, 90.5, 66.1, 43.5],
        [376.4,237.5,191.8,158.4,140.7,126.5,114.7],
        [512.8,331.8,290.3,250.4,220.8,197.1,169.3],
        [672.1,487.3,423.5,390.9,348.5,294.7,232.3],
    ]),
    0.3: np.array([
        [307.0,166.8,128.0,109.2, 85.3, 61.8, 43.4],
        [389.3,247.2,193.7,160.0,139.4,124.7,109.3],
        [524.0,337.3,289.6,254.6,224.1,194.6,162.3],
        [685.0,491.5,421.9,379.7,353.5,292.6,228.3],
    ]),
    0.4: np.array([
        [306.7,164.7,124.0,104.8, 82.1, 60.3, 43.0],
        [398.9,253.7,193.5,157.8,136.3,121.1,105.4],
        [537.6,341.8,288.0,253.4,222.5,191.5,157.1],
        [687.4,489.4,415.9,360.3,352.5,287.9,222.4],
    ]),
    0.5: np.array([
        [306.8,159.3,119.8,100.5, 81.9, 60.5, 44.8],
        [408.4,258.5,193.0,153.8,133.3,118.2,102.9],
        [555.9,348.0,286.9,252.1,220.4,189.3,153.4],
        [686.7,486.3,412.5,342.8,346.1,282.5,215.7],
    ]),
    0.6: np.array([
        [306.8,157.8,117.5, 99.8, 84.7, 61.9, 46.4],
        [418.9,263.0,194.8,153.5,133.2,117.4,102.5],
        [582.3,358.1,288.4,250.1,219.1,189.0,152.3],
        [701.7,488.9,412.9,329.3,338.7,279.3,212.4],
    ]),
}


def make_processing_map(target_strain=0.4):
    sigma     = stress_data[target_strain]
    log_sigma = np.log10(sigma)

    # Step 1: UnivariateSpline k=3 per temperature column
    le_eval = np.linspace(-3, 0, 60)
    m_grid  = np.zeros((len(le_eval), len(temperatures)))

    for j in range(len(temperatures)):
        spl  = UnivariateSpline(log_edot, log_sigma[:, j], k=3, s=0)
        dspl = spl.derivative()
        m_grid[:, j] = np.clip(dspl(le_eval), 1e-4, 0.99)

    eta_grid = 2 * m_grid / (m_grid + 1) * 100

    # Step 2: 2-D spline on eta for smooth contours
    T_fine  = np.linspace(900, 1100, 300)
    le_fine = np.linspace(-3, 0, 300)
    TT, LL  = np.meshgrid(T_fine, le_fine)

    spl2d   = RectBivariateSpline(le_eval, temperatures, eta_grid, kx=3, ky=3)
    eta_map = spl2d(le_fine, T_fine)

    # Step 3: Instability xi by finite difference on smooth eta
    m_map  = (eta_map / 100) / np.clip(2 - eta_map / 100, 1e-6, None)
    dm_dle = np.gradient(m_map, le_fine, axis=0)
    xi_map = (1 / np.clip(m_map, 1e-6, None) - 1 / (m_map + 1)) * dm_dle + m_map

    # Plot
    fig, ax = plt.subplots(figsize=(7, 6))

    levels = list(range(14, 64, 2))
    cs = ax.contour(TT, LL, eta_map, levels=levels,
                    colors='black', linewidths=0.8)
    ax.clabel(cs, inline=True, fontsize=7, fmt='%d', inline_spacing=4)

    ax.contourf(TT, LL, xi_map, levels=[-1e6, 0],
                colors=['lightgrey'], alpha=0.8)
    ax.contour(TT, LL, xi_map, levels=[0],
               colors='black', linewidths=1.0, linestyles='--')

    dom1 = plt.Rectangle((900,-2.6), 110, 1.6,
                          lw=1.5, edgecolor='black', facecolor='none')
    dom2 = plt.Rectangle((1040,-2.6), 65, 1.2,
                          lw=1.5, edgecolor='black', facecolor='none')
    ax.add_patch(dom1); ax.add_patch(dom2)
    ax.text(912, -2.4, 'Domain 1', fontsize=9, fontweight='bold')
    ax.text(1048,-2.4, 'Domain 2', fontsize=9, fontweight='bold')
    ax.text(902, -0.25,'Instability', fontsize=8, rotation=90, style='italic')
    ax.text(960,  0.35, 'Instability', fontsize=8, style='italic')

    ax.set_xlim(900, 1100); ax.set_ylim(-3, 0)
    ax.set_xlabel('Temperature, °C', fontsize=11)
    ax.set_ylabel(r'Log (Strain rate, s$^{-1}$)', fontsize=11)
    ax.set_title(f'IN718 (25 microns)          Strain = {target_strain}',
                 fontsize=11, pad=8)
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax.tick_params(which='both', direction='in', top=True, right=True)
    plt.tight_layout()

    fname = f'processing_map_IN718_strain{target_strain}.png'
    plt.savefig(fname, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fname}")


for s in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]:
    make_processing_map(s)