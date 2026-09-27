"""Exact random-sequential comparison, with discrete null-selection skipping.

One MCS is N independent uniform site selections with replacement. Local
probabilities are unchanged. This is NOT Gillespie or a continuous-time model.
Geometric skipping and state-weighted selection reproduce the discrete chain
exactly while avoiding millions of empty selections. X arrivals and own death
are sampled jointly, preserving their independent Bernoulli marginals.
"""
from pathlib import Path
import json
import numpy as np
from run_simulation import GalacticEcosystem, OFFSETS

ROOT = Path(__file__).resolve().parents[1]


class AsynchronousEcosystem(GalacticEcosystem):
    def rebuild(self):
        self.flat = self.grid.ravel()
        self.pools = [list(np.flatnonzero(self.flat == s)) for s in range(4)]
        self.positions = np.empty(self.grid.size, dtype=int)
        for pool in self.pools:
            for j, cell in enumerate(pool):
                self.positions[cell] = j
        q = self.p_expand * (1-self.s_travel_cost)
        self.q = q
        self.weights = np.array([self.p_birth,
            self.p_mutate+(1-self.p_mutate)*self.decay_introvert,
            q+self.decay_extrovert-q*self.decay_extrovert, self.decay_ruin])
        self.selection = 0
        self.boundary = 0
        self.pending = None

    def change(self, cell, state):
        old = int(self.flat[cell])
        if old == state:
            return
        pool = self.pools[old]
        pos = self.positions[cell]
        last = pool[-1]
        pool[pos] = last
        self.positions[last] = pos
        pool.pop()
        self.positions[cell] = len(self.pools[state])
        self.pools[state].append(cell)
        self.flat[cell] = state

    def event(self):
        masses = self.weights * [len(p) for p in self.pools]
        state = int(np.searchsorted(np.cumsum(masses), self.rng.random()*sum(masses), side='right'))
        cell = self.pools[state][self.rng.integers(len(self.pools[state]))]
        if state == 0:
            self.change(cell, 1)
        elif state == 1:
            self.change(cell, 2 if self.rng.random() < self.p_mutate/self.weights[1] else 0)
        elif state == 3:
            self.change(cell, 0)
        else:
            # Conditional on (arrival OR own death); source can launch then die.
            u = self.rng.random()*self.weights[2]
            both = self.q*self.decay_extrovert
            arrival = u < self.q
            dies = u < both or u >= self.q
            if arrival:
                dy, dx = OFFSETS[self.rng.integers(8)]
                target = ((cell//self.L+dy)%self.L)*self.L+(cell%self.L+dx)%self.L
                target_state = self.flat[target]
                if target_state == 0:
                    self.change(target, 2)
                elif self.rng.random() < self.s_purge:
                    self.change(cell, 3)
                    if target_state == 2:
                        self.change(target, 3)
            if dies:
                self.change(cell, 3)

    def step(self):
        if not hasattr(self, 'pools'):
            self.rebuild()
        self.boundary += self.grid.size
        while True:
            if self.pending is None:
                event_probability = float(self.weights @ np.array([len(p) for p in self.pools])) / self.grid.size
                if event_probability == 0:
                    break
                self.pending = self.selection + int(self.rng.geometric(min(event_probability,1.)))
            if self.pending > self.boundary:
                break
            self.selection = self.pending
            self.pending = None
            self.event()
        for s, name in enumerate(self.history):
            self.history[name].append(len(self.pools[s])/self.grid.size)


def spatial(grid):
    """Axial connected indicator correlations, r=1,2,4,8; population mean removed."""
    rows=[]
    for state in (1,2):
        v=(grid==state).astype(float)
        mu=v.mean()
        rows.append([float(((v*np.roll(v,r,0)).mean()+(v*np.roll(v,r,1)).mean())/2-mu**2)
                     for r in (1,2,4,8)])
    return rows


def measure(cls, params, seed, steps=20000, burn=10000):
    sim=cls(seed=seed,**params)
    correlations=[]
    for step in range(steps):
        sim.step()
        if step>=burn and (step+1)%100==0:
            correlations.append(spatial(sim.grid))
    h=np.array(list(sim.history.values()))[:,burn:]
    ac=[]
    for row in h[1:3]:
        ac.append([float(np.corrcoef(row[:-lag],row[lag:])[0,1]) for lag in (1,10,100)])
    return dict(mean=h.mean(axis=1).tolist(),within_run_variance=h.var(axis=1).tolist(),
                spatial_covariance=np.mean(correlations,axis=0).tolist(),
                temporal_correlation=ac)


def lineage(cls, params, seed, cap=1000):
    # Every X descends from the single founder: emergence and mutation disabled.
    sim=cls(seed=seed, p_birth=0,p_mutate=0,**params)
    sim.grid[sim.L//2,sim.L//2]=2
    for t in range(1,cap+1):
        sim.step()
        if not np.any(sim.grid==2):
            return t,False
    return cap,True


def main():
    results={}
    for name,params in [('baseline',{}),('higher_expansion',dict(p_expand=.1,s_travel_cost=.05))]:
        results[name]={}
        for mode,cls in [('synchronous',GalacticEcosystem),('asynchronous',AsynchronousEcosystem)]:
            runs=[measure(cls,params,seed) for seed in range(10)]
            means=np.array([r['mean'] for r in runs])
            trials=[lineage(cls,params,10000+s) for s in range(200)]
            times=np.array([v[0] for v in trials])
            censored=np.array([v[1] for v in trials])
            results[name][mode]=dict(runs=runs,mean=means.mean(axis=0).tolist(),
                between_seed_sd=means.std(axis=0,ddof=1).tolist(),
                lineage_times=times.tolist(),lineage_right_censored=censored.tolist(),
                lineage_survival_at={str(t):float(np.mean(times>t)) for t in (10,50,100,500)},
                lineage_median_observed_time=float(np.median(times)))
            print(name,mode,results[name][mode]['mean'],'censored',int(censored.sum()),flush=True)
    payload=dict(protocol=dict(grid_size=60,steps=20000,burn_in=10000,seeds=list(range(10)),
        spatial_sampling_every=100,spatial_r=[1,2,4,8],temporal_lags=[1,10,100],
        lineage_seeds=list(range(10000,10200)),lineage_cap=1000,
        lineage_initial='single X, otherwise empty; emergence and mutation disabled',
        async_clock='N uniform site selections with replacement per MCS; exact geometric null skipping'),results=results)
    (ROOT/'results'/'update_comparison.json').write_text(json.dumps(payload,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(12,3.5),layout='constrained')
    for mode,color in [('synchronous','#2874a6'),('asynchronous','#c0392b')]:
        r=results['baseline'][mode]
        offset=-.12 if mode=='synchronous' else .12
        axes[0].errorbar(np.arange(3)+offset,np.array(r['mean'])[1:]*100,
            yerr=np.array(r['between_seed_sd'])[1:]*100,fmt='o',label=mode,color=color)
        cov=np.mean([v['spatial_covariance'] for v in r['runs']],axis=0)
        axes[1].plot([1,2,4,8],cov[1],'-o',color=color,label=mode)
        t=np.arange(0,501)
        axes[2].step(t,[np.mean(np.array(r['lineage_times'])>v) for v in t],where='post',color=color,label=mode)
    axes[0].set(xticks=[0,1,2],xticklabels=['Inward','Expanding','Remnant'],ylabel='Occupancy (%)',title='Baseline; mean and between-seed SD')
    axes[1].set(xlabel='Axial lattice separation',ylabel='Connected X-indicator covariance',title='Late-time spatial correlation')
    axes[2].set(xlabel='Step / MCS (integer observation)',ylabel='Fraction surviving',title='Single-founder experiment; 200 trials',ylim=(0,1))
    axes[0].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=.2)
    for suffix in ('pdf','png'): fig.savefig(ROOT/'figures'/f'fig4_update_comparison.{suffix}',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    main()
