"""Robustness experiments for the synchronous automaton.

  size           lattice-size dependence, L = 15 ... 240
  sweep          one-at-a-time parameter sweep and a (p, delta_X) grid
  vulnerability  relaxing the protection of inward sites (v_intro, ruin/capture)
  conflict       conflict lethality s_purge x launch probability p at delta_X=0.025

Usage: python robustness.py [size|sweep|vulnerability|conflict|figures|all] [--quick]
Each experiment writes results/robustness_<name>.json; `figures` draws
figures/fig6-8 from those files. Occupancy protocol matches the main runs:
20,000 steps from an empty lattice, means over steps 10,001-20,000, seeds
listed in each file. Lineage trials start from one expanding site with
emergence and mutation disabled, right-censored at 1000 steps.
"""
from pathlib import Path
from multiprocessing import Pool
import argparse
import json
import os
import numpy as np
from run_simulation import GalacticEcosystem

ROOT = Path(__file__).resolve().parents[1]
STATES = ['empty', 'introvert', 'extrovert', 'ruin']
BASE = dict(p_birth=1e-4, p_mutate=5e-4, p_expand=.05, s_travel_cost=.3, s_purge=.99,
            decay_introvert=1e-4, decay_extrovert=.05, decay_ruin=.005)
CASES = {'baseline': {}, 'higher_expansion': dict(p_expand=.1, s_travel_cost=.05)}
QUICK = False


def occupancy(job):
    params, L, seed, steps, burn = job
    sim = GalacticEcosystem(grid_size=L, seed=seed, **params).run(steps)
    h = np.array([sim.history[s] for s in STATES])[:, burn:]
    return h.mean(axis=1).tolist()


def lineage(job):
    params, L, seed, cap = job
    params = {**params, 'p_birth': 0., 'p_mutate': 0.}
    sim = GalacticEcosystem(grid_size=L, seed=seed, **params)
    sim.grid[L // 2, L // 2] = 2
    for t in range(1, cap + 1):
        sim.step()
        if not np.any(sim.grid == 2):
            return t, False
    return cap, True


def run_jobs(fn, jobs):
    with Pool(os.cpu_count()) as pool:
        return pool.map(fn, jobs, chunksize=1)


def summarize_occ(values):
    v = np.asarray(values)
    return dict(per_seed=v.tolist(), mean=v.mean(axis=0).tolist(),
                between_seed_sd=v.std(axis=0, ddof=1).tolist())


def summarize_lin(values, cap):
    t = np.array([x[0] for x in values]); c = np.array([x[1] for x in values])
    return dict(n=len(t), censored=int(c.sum()), median=float(np.median(t)),
                max=int(t.max()), survival_at={str(k): float(np.mean(t > k)) for k in (10, 100, 500)},
                surviving_at_cap=float(c.mean()), cap=cap)


def protocol(seeds, steps, burn, lineage_trials, cap):
    return dict(seeds=list(seeds), steps=steps, burn_in=burn, lineage_trials=lineage_trials,
                lineage_cap=cap, lineage_seeds_start=10000, base_parameters=BASE,
                quick=QUICK, update='synchronous')


def settings(seeds):
    if QUICK:
        return range(2), 2000, 1000, 20, 300
    return range(seeds), 20000, 10000, 200, 1000


def evaluate(params, L, seeds, steps, burn, trials, cap):
    occ = run_jobs(occupancy, [(params, L, s, steps, burn) for s in seeds])
    lin = run_jobs(lineage, [(params, L, 10000 + s, cap) for s in range(trials)]) if trials else None
    out = dict(parameters=params, grid_size=L, occupancy=summarize_occ(occ))
    if lin is not None:
        out['lineage'] = summarize_lin(lin, cap)
    q = params['p_expand'] * (1 - params['s_travel_cost'])
    out['R_upper'] = q / params['decay_extrovert'] if params['decay_extrovert'] else None
    loss = params['p_mutate'] + (1 - params['p_mutate']) * params['decay_introvert']
    out['inward_envelope'] = params['p_birth'] / (params['p_birth'] + loss)
    return out


def save(name, payload):
    (ROOT / 'results').mkdir(exist_ok=True)
    (ROOT / 'results' / f'robustness_{name}.json').write_text(json.dumps(payload, indent=1) + '\n')


def size_experiment():
    seeds, steps, burn, trials, cap = settings(10)
    sizes = [15, 30, 60, 120] if QUICK else [15, 30, 60, 120, 240]
    rows = []
    for case, change in CASES.items():
        for L in sizes:
            r = evaluate({**BASE, **change}, L, seeds, steps, burn, trials, cap)
            r['case'] = case
            rows.append(r)
            print('size', case, L, np.round(np.array(r['occupancy']['mean']) * 100, 3), r['lineage']['median'], flush=True)
    save('size', dict(protocol=protocol(seeds, steps, burn, trials, cap), rows=rows))


SWEEP = {
    'p_birth': [2.5e-5, 5e-5, 1e-4, 2e-4, 4e-4],
    'p_mutate': [1.25e-4, 2.5e-4, 5e-4, 1e-3, 2e-3],
    'decay_introvert': [2.5e-5, 5e-5, 1e-4, 2e-4, 4e-4],
    'p_expand': [.0125, .025, .05, .1, .2],
    's_travel_cost': [0., .15, .3, .5, .7],
    's_purge': [.1, .3, .6, .9, .99],
    'decay_extrovert': [.0125, .025, .05, .1, .2],
    'decay_ruin': [.00125, .0025, .005, .01, .02],
}
GRID_P = [.025, .05, .1, .2, .4, .8]
GRID_DX = [.0125, .025, .05, .1, .2]


def sweep_experiment():
    seeds, steps, burn, trials, cap = settings(5)
    trials = trials // 2
    oat = []
    for name, values in SWEEP.items():
        for value in values:
            r = evaluate({**BASE, name: value}, 60, seeds, steps, burn, trials, cap)
            r.update(swept=name, value=value)
            oat.append(r)
            print('sweep', name, value, np.round(np.array(r['occupancy']['mean']) * 100, 3), flush=True)
    grid = []
    for p in GRID_P:
        for dx in GRID_DX:
            r = evaluate({**BASE, 'p_expand': p, 'decay_extrovert': dx}, 60, seeds, steps, burn, trials, cap)
            grid.append(r)
            print('grid', p, dx, round(r['R_upper'], 2), np.round(np.array(r['occupancy']['mean']) * 100, 3),
                  r['lineage']['surviving_at_cap'], flush=True)
    save('sweep', dict(protocol=protocol(seeds, steps, burn, trials, cap), one_at_a_time=oat,
                       grid=dict(p_expand=GRID_P, decay_extrovert=GRID_DX, rows=grid)))


V_VALUES = [0., .001, .01, .03, .1, .3, 1.]


def vulnerability_experiment():
    seeds, steps, burn, trials, cap = settings(10)
    rows = []
    for case, change in CASES.items():
        for invasion in ('ruin', 'capture'):
            for v in V_VALUES:
                if v == 0 and invasion == 'capture':
                    continue  # identical to v=0 'ruin' draw-for-draw
                params = {**BASE, **change, 'v_intro': v, 'invasion': invasion}
                r = evaluate(params, 60, seeds, steps, burn, 0, cap)
                r.update(case=case, invasion=invasion, v_intro=v)
                rows.append(r)
                print('vuln', case, invasion, v, np.round(np.array(r['occupancy']['mean']) * 100, 3), flush=True)
    save('vulnerability', dict(protocol=protocol(seeds, steps, burn, 0, cap), rows=rows))


CONFLICT_U = [0., .1, .3, .5, .7, .9, .99]
CONFLICT_P = [.05, .1, .2, .4]


def conflict_experiment():
    # Conflict applies to ANY occupied target, including the lineage's own
    # expanding sites and inert remnants, so s_purge also sets self-conflict.
    seeds, steps, burn, trials, cap = settings(5)
    trials = trials // 2
    rows = []
    for p in CONFLICT_P:
        for u in CONFLICT_U:
            r = evaluate({**BASE, 'p_expand': p, 'decay_extrovert': .025, 's_purge': u},
                         60, seeds, steps, burn, trials, cap)
            rows.append(r)
            print('conflict', p, u, round(r['R_upper'], 2), np.round(np.array(r['occupancy']['mean']) * 100, 3),
                  r['lineage']['surviving_at_cap'], flush=True)
    save('conflict', dict(protocol=protocol(seeds, steps, burn, trials, cap), s_purge=CONFLICT_U,
                          p_expand=CONFLICT_P, decay_extrovert=.025, rows=rows))


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = {'baseline': '#2874a6', 'higher_expansion': '#c0392b'}
    load = lambda n: json.loads((ROOT / 'results' / f'robustness_{n}.json').read_text())

    d = load('size')
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6), layout='constrained')
    for case, color in colors.items():
        rows = [r for r in d['rows'] if r['case'] == case]
        L = [r['grid_size'] for r in rows]
        for k, (state, axis) in enumerate([(1, ax[0]), (2, ax[1])]):
            m = np.array([r['occupancy']['mean'][state] for r in rows]) * 100
            s = np.array([r['occupancy']['between_seed_sd'][state] for r in rows]) * 100
            axis.errorbar(L, m, yerr=s, fmt='o-', color=color, label=case, capsize=3)
        ax[2].plot(L, [r['lineage']['median'] for r in rows], 'o-', color=color, label=case)
    ax[0].set(xscale='log', xlabel='Lattice side L', ylabel='Inward occupancy (%)', title='Inward fraction vs lattice size')
    ax[1].set(xscale='log', xlabel='Lattice side L', ylabel='Expanding occupancy (%)', title='Expanding fraction vs lattice size')
    ax[2].set(xscale='log', xlabel='Lattice side L', ylabel='Median lineage lifetime (steps)', title='Single-founder lineages')
    for a in ax:
        a.minorticks_off()
        a.set_xticks(L, labels=[str(v) for v in L]); a.grid(alpha=.2)
    ax[0].legend(fontsize=8)
    for ext in ('pdf', 'png'):
        fig.savefig(ROOT / 'figures' / f'fig6_lattice_size.{ext}', dpi=180)
    plt.close(fig)

    d = load('sweep')
    fig = plt.figure(figsize=(12, 7.5), layout='constrained')
    sub = fig.subfigures(1, 2, width_ratios=[1.45, 1])
    axs = sub[0].subplots(4, 2, sharey=True)
    for a, name in zip(axs.ravel(), SWEEP):
        rows = [r for r in d['one_at_a_time'] if r['swept'] == name]
        x = [r['value'] for r in rows]
        for state, color, label in [(1, '#2874a6', 'inward'), (2, '#c0392b', 'expanding x10')]:
            scale = 100 * (10 if state == 2 else 1)
            m = np.array([r['occupancy']['mean'][state] for r in rows]) * scale
            s = np.array([r['occupancy']['between_seed_sd'][state] for r in rows]) * scale
            a.errorbar(x, m, yerr=s, fmt='o-', ms=3, color=color, label=label)
        a.plot(x, np.array([r['inward_envelope'] for r in rows]) * 100, ':', color='#2874a6', lw=1, label='b/(b+l) envelope')
        a.axvline(BASE[name], color='grey', lw=.6)
        if name not in ('s_travel_cost', 's_purge'):
            a.set_xscale('log')
            a.minorticks_off()
            a.set_xticks(x, labels=[f'{v:g}' for v in x])
        a.set_title(name, fontsize=9); a.grid(alpha=.2); a.tick_params(labelsize=7)
    axs[0, 0].legend(fontsize=6)
    sub[0].supylabel('Occupancy (%)', fontsize=9)
    ax = sub[1].subplots(2, 1)
    rows = d['grid']['rows']
    P, DX = d['grid']['p_expand'], d['grid']['decay_extrovert']
    X = np.array([[next(r for r in rows if r['parameters']['p_expand'] == p and r['parameters']['decay_extrovert'] == dx)
                   for p in P] for dx in DX], dtype=object)
    for a, key, title, fmt in [(ax[0], lambda r: r['occupancy']['mean'][2] * 100, 'Expanding occupancy (%)', '{:.2g}'),
                               (ax[1], lambda r: r['lineage']['median'], 'Median lineage lifetime (none survived to cap)', '{:.0f}')]:
        Z = np.vectorize(key)(X).astype(float)
        im = a.imshow(Z, origin='lower', cmap='magma', aspect='auto')
        for i in range(len(DX)):
            for j in range(len(P)):
                a.text(j, i, fmt.format(Z[i, j]), ha='center', va='center', fontsize=6,
                       color='white' if Z[i, j] < (Z.max() + Z.min()) / 2 else 'black')
        a.set_xticks(range(len(P)), labels=[str(p) for p in P], fontsize=7)
        a.set_yticks(range(len(DX)), labels=[str(v) for v in DX], fontsize=7)
        a.set(xlabel='p (launch probability)', ylabel='delta_X', title=title)
        Ru = np.array([[p * .7 / dx for p in P] for dx in DX])
        a.contour(Ru, levels=[1], colors='cyan', linewidths=1)
    sub[1].suptitle('s=0.3; cyan: R_upper = 1', fontsize=9)
    for ext in ('pdf', 'png'):
        fig.savefig(ROOT / 'figures' / f'fig7_parameter_sweep.{ext}', dpi=180)
    plt.close(fig)

    d = load('vulnerability')
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6), layout='constrained')
    for case, color in colors.items():
        for invasion, ls in [('ruin', '-'), ('capture', '--')]:
            rows = [r for r in d['rows'] if r['case'] == case and (r['invasion'] == invasion or r['v_intro'] == 0)]
            rows.sort(key=lambda r: r['v_intro'])
            x = np.array([max(r['v_intro'], 3e-4) for r in rows])
            for state, a in [(1, ax[0]), (2, ax[1]), (3, ax[2])]:
                m = np.array([r['occupancy']['mean'][state] for r in rows]) * 100
                s = np.array([r['occupancy']['between_seed_sd'][state] for r in rows]) * 100
                a.errorbar(x, m, yerr=s, fmt='o' + ls, ms=3, color=color, label=f'{case}, {invasion}', capsize=2)
    for a, t in zip(ax, ['Inward', 'Expanding', 'Remnant']):
        a.set(xscale='log', xlabel='v_intro (v=0 plotted at 3e-4)', ylabel=f'{t} occupancy (%)', title=f'{t} fraction')
        a.grid(alpha=.2)
    ax[0].legend(fontsize=7)
    for ext in ('pdf', 'png'):
        fig.savefig(ROOT / 'figures' / f'fig8_inward_vulnerability.{ext}', dpi=180)
    plt.close(fig)

    d = load('conflict')
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6), layout='constrained')
    shades = plt.cm.viridis(np.linspace(0, .85, len(d['p_expand'])))
    for p, color in zip(d['p_expand'], shades):
        rows = sorted([r for r in d['rows'] if r['parameters']['p_expand'] == p], key=lambda r: r['parameters']['s_purge'])
        u = [r['parameters']['s_purge'] for r in rows]
        label = f'p={p} (R_upper={p * .7 / .025:.1f})'
        for state, a in [(1, ax[0]), (2, ax[1])]:
            m = np.array([r['occupancy']['mean'][state] for r in rows]) * 100
            sd = np.array([r['occupancy']['between_seed_sd'][state] for r in rows]) * 100
            a.errorbar(u, m, yerr=sd, fmt='o-', ms=3, color=color, label=label, capsize=2)
        ax[2].plot(u, [r['lineage']['surviving_at_cap'] for r in rows], 'o-', ms=3, color=color, label=label)
    ax[0].set(xlabel='Conflict probability s_purge', ylabel='Inward occupancy (%)', title='Inward fraction')
    ax[1].set(xlabel='Conflict probability s_purge', ylabel='Expanding occupancy (%)', title='Expanding fraction', yscale='log')
    ax[2].set(xlabel='Conflict probability s_purge', ylabel='Fraction alive at cap', title='Single-founder lineages (cap 1000)', ylim=(-.03, 1.03))
    for a in ax:
        a.grid(alpha=.2); a.axvline(.99, color='grey', lw=.6)
    ax[2].legend(fontsize=7)
    for ext in ('pdf', 'png'):
        fig.savefig(ROOT / 'figures' / f'fig9_conflict.{ext}', dpi=180)
    plt.close(fig)


def main():
    global QUICK
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('what', choices=['size', 'sweep', 'vulnerability', 'conflict', 'figures', 'all'])
    parser.add_argument('--quick', action='store_true', help='tiny smoke-test protocol; do not commit its output')
    args = parser.parse_args()
    QUICK = args.quick
    todo = ['size', 'sweep', 'vulnerability', 'conflict', 'figures'] if args.what == 'all' else [args.what]
    for name in todo:
        dict(size=size_experiment, sweep=sweep_experiment,
             vulnerability=vulnerability_experiment, conflict=conflict_experiment,
             figures=figures)[name]()


if __name__ == '__main__':
    main()
