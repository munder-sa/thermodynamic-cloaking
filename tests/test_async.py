import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from compare_updates import AsynchronousEcosystem, spatial


class AsyncTests(unittest.TestCase):
    def test_exact_finite_selection_two_state_expectation(self):
        # Independent analytic expectation for 18 individual selections.
        b,d,N,steps=.2,.3,9,2
        expected=b/(b+d)*(1-(1-(b+d)/N)**(N*steps))
        values=[]
        for seed in range(1000):
            sim=AsynchronousEcosystem(grid_size=3,p_birth=b,p_mutate=0,
                                      decay_introvert=d,seed=seed).run(steps)
            values.append(np.mean(sim.grid==1))
        se=np.std(values,ddof=1)/np.sqrt(len(values))
        self.assertLess(abs(np.mean(values)-expected),4*se)

    def test_compare_naive_random_selection_distribution(self):
        values=[]
        offsets=[(y,x) for y in (-1,0,1) for x in (-1,0,1) if (y,x)!=(0,0)]
        params=dict(grid_size=3,p_birth=.2,p_mutate=.15,p_expand=.6,
                    s_travel_cost=.2,s_purge=.7,decay_introvert=.1,
                    decay_extrovert=.25,decay_ruin=.2)
        exact=[]
        for seed in range(600):
            sim=AsynchronousEcosystem(seed=seed,**params).run(8)
            exact.append(np.bincount(sim.grid.ravel(),minlength=4)/9)
            rng=np.random.default_rng(seed+20000)
            g=np.zeros((3,3),dtype=int)
            for _ in range(72):
                y,x=divmod(int(rng.integers(9)),3)
                state=g[y,x]
                if state==0:
                    if rng.random()<.2:g[y,x]=1
                elif state==1:
                    if rng.random()<.15:g[y,x]=2
                    elif rng.random()<.1:g[y,x]=0
                elif state==3:
                    if rng.random()<.2:g[y,x]=0
                else:
                    if rng.random()<.6:
                        dy,dx=offsets[int(rng.integers(8))]
                        ty,tx=(y+dy)%3,(x+dx)%3
                        if rng.random()<.8:
                            target=g[ty,tx]
                            if target==0:g[ty,tx]=2
                            elif rng.random()<.7:
                                g[y,x]=3
                                if target==2:g[ty,tx]=3
                    if rng.random()<.25:g[y,x]=3
            values.append(np.bincount(g.ravel(),minlength=4)/9)
        a,b=np.array(exact),np.array(values)
        se=np.sqrt(a.var(axis=0,ddof=1)/600+b.var(axis=0,ddof=1)/600)
        self.assertTrue(np.all(np.abs(a.mean(axis=0)-b.mean(axis=0))<4*se))

    def test_joint_launch_and_death(self):
        sim=AsynchronousEcosystem(grid_size=3,p_birth=0,p_mutate=0,p_expand=1,
                                  s_travel_cost=0,decay_extrovert=1)
        sim.grid[1,1]=2
        sim.rebuild()
        sim.event()
        self.assertEqual(sim.grid[1,1],3)
        self.assertEqual(np.count_nonzero(sim.grid==2),1)

    def test_pool_integrity_and_reproducibility(self):
        a=AsynchronousEcosystem(grid_size=8,p_birth=.2,p_mutate=.1,seed=3).run(100)
        b=AsynchronousEcosystem(grid_size=8,p_birth=.2,p_mutate=.1,seed=3).run(100)
        np.testing.assert_array_equal(a.grid,b.grid)
        for state,pool in enumerate(a.pools):
            self.assertEqual(set(pool),set(np.flatnonzero(a.grid.ravel()==state)))
            for j,cell in enumerate(pool):self.assertEqual(a.positions[cell],j)

    def test_absorbing_empty_and_constant_spatial_covariance(self):
        sim=AsynchronousEcosystem(grid_size=3,p_birth=0).run(10)
        self.assertFalse(np.any(sim.grid))
        np.testing.assert_allclose(spatial(np.ones((8,8))),0)

if __name__=='__main__':unittest.main()
