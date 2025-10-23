from __future__ import annotations
import math
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl

from ion_simulation_plan.plan_sims import build_plan
from conc2time import _load_records, _density_stats, _performance_stats, _load_microstate_config

ION_VALENCE = {
    'NA': 1,
    'K': 1,
    'LI': 1,
    'CA': 2,
    'MG': 2,
}

def main():
    # Set typographic style
    mpl.rcParams.update({
        'font.family': 'serif',
        'font.size': 10,
        'mathtext.fontset': 'stix',
    })
    process_root = Path('ip6-run/process')
    microstate = 'IP_010101'
    sim_ns = 100.0

    # Compute constants once from historical dataset (median density and k)
    records = _load_records(process_root)
    if not records:
        raise SystemExit(f"No historical simulations found under {process_root}")
    dens = _density_stats(records)
    perf = _performance_stats(records)
    rho_med = dens.median         # atoms per nm^3
    k_med = perf.median_k         # (ns/day) * atoms

    # Microstate parameters (positive ion count, box type, charge unused here)
    cfg = _load_microstate_config(process_root, microstate)
    Nplus_default = cfg.positive_count
    box_type = cfg.box_type  # assumed cubic in these datasets

    # We plot only water-model dependence across a shared concentration range.
    cmin, cmax = 0.5, 300.0  # mM

    def wall_time_days(c_mM: float, multiplier: float, Nplus: int = Nplus_default) -> tuple[float, float]:
        # Convert concentration and compute volume in nm^3
        c = c_mM / 1000.0  # mol/L
        V_nm3 = (Nplus / 6.02214076e23) / c * 1e24
        # Estimate atoms and throughput
        Natoms = rho_med * V_nm3
        ns_per_day = k_med / Natoms
        hours = 24.0 * sim_ns / ns_per_day
        return hours / 24.0 * multiplier, V_nm3 ** (1.0 / 3.0)

    water_models = ['TIP3P', 'TIP4P/2005']
    colors = {'TIP3P': '#1f77b4', 'TIP4P/2005': '#d62728'}
    markers = {'TIP3P': 'o', 'TIP4P/2005': 's'}

    # 8 cm wide figure (height ~5.5 cm)
    cm = 1 / 2.54
    plt.figure(figsize=(8.0 * cm, 5.5 * cm))
    import numpy as np
    xs = np.logspace(math.log10(cmin), math.log10(cmax), 400)
    for wm in water_models:
        mult = 1.0 if wm == 'TIP3P' else (4.0 / 3.0)
        ys_days = [wall_time_days(x, mult)[0] for x in xs]
        plt.plot(xs, ys_days,
                 #marker=markers.get(wm, 'o'),
                 #markevery=[-1],
                 color=colors.get(wm, 'k'),
                 linestyle='-', label=wm)

    plt.xscale('log')  # spans 0.5 mM to 300 mM nicely
    plt.xlabel('Concentration (mM)')
    plt.ylabel('t (days) per 100 ns')
    plt.axhline(10, color='gray', linestyle='--', alpha=0.5, label='10 days')
    #plt.title('Wall time vs. concentration (TIP3P vs TIP4P/2005, 100 ns, 310 K)')
    plt.grid(True, which='both', linestyle=':', alpha=0.4)
    plt.legend(ncol=1, fontsize=8)
    plt.tight_layout()
    out = Path('ion_simulation_plan/walltime_vs_conc.pdf')
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out)
    print(f'Wrote {out}')

if __name__ == '__main__':
    main()
