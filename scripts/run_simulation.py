"""Synchronous four-state stochastic cellular automaton (not Moran).

All proposals use the old grid. Colonization of an old empty site wins
over simultaneous spontaneous birth; any ruin proposal wins over survival.
An old expansionist may launch once in its final step. Incoming colonies
act only next step. The chosen rules are model assumptions, not empirical
astrophysical rates. Run this file from any working directory.

Optional inward vulnerability (default off, reproducing earlier results
draw-for-draw): when a conflict occurs at an old inward target, the target
is also lost with probability v_intro. With invasion='ruin' it becomes a
remnant and the source is also ruined (mutual destruction, as for two
expansionists). With invasion='capture' it becomes an expanding site and
the source survives the encounter; a failed capture still ruins the source.
"""
from pathlib import Path
import argparse
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OFFSETS = np.array([(y, x) for y in (-1, 0, 1) for x in (-1, 0, 1) if (y, x) != (0, 0)])


def resolve_proposals(local, colony_targets, ruin_targets):
    """Order-independent union resolution, including duplicate proposals."""
    result = local.copy().ravel()
    result[np.asarray(colony_targets, dtype=int)] = 2
    result[np.asarray(ruin_targets, dtype=int)] = 3
    return result.reshape(local.shape)


class GalacticEcosystem:
    EMPTY, INTROVERT, EXTROVERT, RUIN = range(4)
    def __init__(self, grid_size=60, p_birth=0.0001, p_mutate=0.0005,
                 p_expand=0.05, s_travel_cost=0.3, s_purge=0.99,
                 decay_introvert=0.0001, decay_extrovert=0.05,
                 decay_ruin=0.005, v_intro=0.0, invasion='ruin', seed=0):
        if not isinstance(grid_size, int) or grid_size < 3:
            raise ValueError('grid_size must be an integer >=3 (8 distinct neighbors)')
        self.L = grid_size
        self.params = dict(p_birth=p_birth, p_mutate=p_mutate, p_expand=p_expand,
                           s_travel_cost=s_travel_cost, s_purge=s_purge,
                           decay_introvert=decay_introvert, decay_extrovert=decay_extrovert,
                           decay_ruin=decay_ruin, v_intro=v_intro)
        for name, value in self.params.items():
            if not np.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f'{name} must lie in [0,1]')
            setattr(self, name, value)
        if invasion not in ('ruin', 'capture'):
            raise ValueError("invasion must be 'ruin' or 'capture'")
        self.invasion = invasion
        if v_intro > 0:
            self.params['invasion'] = invasion
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.grid = np.zeros((self.L, self.L), dtype=np.int8)
        self.history = {name: [] for name in ('empty', 'introvert', 'extrovert', 'ruin')}

    def step(self):
        old = self.grid.ravel()
        local = old.copy()
        empty, intro, extro, ruin = (np.flatnonzero(old == s) for s in range(4))
        rand = self.rng.random
        local[empty[rand(len(empty)) < self.p_birth]] = self.INTROVERT
        # Mutation gets priority. Therefore unconditional inward loss is
        # p_mutate + (1-p_mutate)*decay_introvert, not the simple sum.
        mutate = rand(len(intro)) < self.p_mutate
        local[intro[mutate]] = self.EXTROVERT
        local[intro[(~mutate) & (rand(len(intro)) < self.decay_introvert)]] = self.EMPTY
        local[ruin[rand(len(ruin)) < self.decay_ruin]] = self.EMPTY
        ruin_targets = list(extro[rand(len(extro)) < self.decay_extrovert])
        sources = extro[rand(len(extro)) < self.p_expand]
        directions = OFFSETS[self.rng.integers(8, size=len(sources))]
        targets = ((sources // self.L + directions[:, 0]) % self.L) * self.L + ((sources % self.L + directions[:, 1]) % self.L)
        arrives = rand(len(sources)) < (1 - self.s_travel_cost)
        target_states = old[targets]
        colonies = targets[arrives & (target_states == self.EMPTY)]
        conflicts = arrives & (target_states != self.EMPTY) & (rand(len(sources)) < self.s_purge)
        ruin_targets.extend(targets[conflicts & (target_states == self.EXTROVERT)])
        if self.v_intro > 0:
            # Drawn only when enabled, so v_intro=0 keeps the original stream.
            lost = conflicts & (target_states == self.INTROVERT) & (rand(len(sources)) < self.v_intro)
            if self.invasion == 'capture':
                colonies = np.concatenate([colonies, targets[lost]])
                conflicts = conflicts & ~lost
            else:
                ruin_targets.extend(targets[lost])
        ruin_targets.extend(sources[conflicts])
        self.grid = resolve_proposals(local.reshape(self.grid.shape), colonies, ruin_targets)
        for state, name in enumerate(self.history):
            self.history[name].append(float(np.count_nonzero(self.grid == state) / self.grid.size))

    def run(self, steps=20000):
        for _ in range(steps):
            self.step()
        return self


def ensemble(seeds, steps, burn_in, **params):
    if not 0 <= burn_in < steps:
        raise ValueError('require 0 <= burn_in < steps')
    runs, first = [], None
    for seed in seeds:
        sim = GalacticEcosystem(seed=seed, **params).run(steps)
        first = sim if first is None else first
        runs.append([np.mean(values[burn_in:]) for values in sim.history.values()])
    values = np.asarray(runs)
    return first, dict(parameters=first.params, grid_size=first.L, seeds=list(seeds),
                      steps=steps, burn_in=burn_in, states=list(first.history),
                      per_seed_time_means=values.tolist(), mean=values.mean(axis=0).tolist(),
                      between_seed_sd=values.std(axis=0, ddof=1).tolist() if len(seeds)>1 else None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steps', type=int, default=20000)
    parser.add_argument('--burn-in', type=int, default=10000)
    parser.add_argument('--seeds', type=int, default=10)
    args = parser.parse_args()
    if args.seeds < 2:
        parser.error('--seeds must be >=2')
    results = {}
    first = None
    for name, changes in [('baseline', {}), ('higher_expansion', {'p_expand': .1, 's_travel_cost': .05})]:
        sim, result = ensemble(range(args.seeds), args.steps, args.burn_in, **changes)
        results[name] = result
        first = sim if first is None else first
        print(name, result['mean'], flush=True)
    out = ROOT / 'results'
    out.mkdir(exist_ok=True)
    (out / 'simulation_summary.json').write_text(json.dumps(results, indent=2) + '\n')
    for name, history in first.history.items():
        np.save(out / f'baseline_seed0_{name}.npy', history)
    np.save(out / 'baseline_seed0_grid.npy', first.grid)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
    axes[0].imshow(first.grid, cmap=ListedColormap(['#f0f1f3', '#2874a6', '#c23b22', '#737373']), vmin=0, vmax=3)
    axes[0].set(title=f'Baseline, seed 0, step {args.steps}', xticks=[], yticks=[])
    for name, color in [('introvert', '#2874a6'), ('extrovert', '#c23b22'), ('ruin', '#737373')]:
        axes[1].plot(np.arange(1, args.steps+1), first.history[name], label=name, color=color, lw=1)
    axes[1].axvline(args.burn_in, ls=':', color='black', label='averaging starts')
    axes[1].set(xlabel='Dimensionless step', ylabel='Site fraction', title='Illustrative cellular automaton')
    axes[1].legend()
    (ROOT / 'figures').mkdir(exist_ok=True)
    for suffix in ('pdf', 'png'):
        fig.savefig(ROOT / 'figures' / f'fig1_dynamics.{suffix}', dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
