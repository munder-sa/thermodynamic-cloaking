"""Independent physical identities and adversarial transition checks."""
import sys
from pathlib import Path
import unittest
import numpy as np
from scipy.integrate import quad
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import physics as p
from run_simulation import GalacticEcosystem, resolve_proposals


class PhysicsTests(unittest.TestCase):
    def test_landauer_and_area_regressions(self):
        self.assertAlmostEqual(float(p.waste_power()) / 9.56992961692908e17, 1., places=12)
        self.assertAlmostEqual(float(p.radiator_area()) / 1.68770683e21, 1., places=7)

    def test_planck_integrates_to_stefan_boltzmann(self):
        # Integrate over dimensionless frequency, independently of flux code.
        for T in (3., 10., 30., 300.):
            scale = p.K_B*T/p.H
            integral = quad(lambda x: np.pi*float(p.planck_b_nu(scale*x,T))*scale,
                            1e-9, 100, epsabs=1e-16)[0]
            self.assertAlmostEqual(integral/(p.SIGMA*T**4),1.,places=8)

    def test_bolometric_flux_and_background_contrast(self):
        for bg in (0., 2.725):
            scale = p.K_B*10/p.H
            integral = quad(lambda x: float(p.flux_ujy(p.C/(scale*x)*1e6,T_bg=bg))*1e-32*scale,
                            1e-7,100,epsabs=1e-25)[0]
            expected = p.waste_power()/(4*np.pi*(10*p.PC_TO_M)**2)
            self.assertAlmostEqual(integral/expected,1.,places=8)

    def test_flux_regressions_and_scaling(self):
        for wavelength, expected in [(100,.0993377485),(183,11.06522225),(235,29.81564538),(450,82.33171356),(850,64.65297312)]:
            self.assertLess(abs(float(p.flux_ujy(wavelength))/expected-1),1e-8)
        self.assertAlmostEqual(float(p.flux_ujy(850,I_dot=2e40,d_pc=20)/p.flux_ujy(850)),.5)

    def test_wien_and_rj_limits(self):
        nu,T = 1e3,300.
        self.assertAlmostEqual(float(p.planck_b_nu(nu,T))/(2*p.K_B*T*nu**2/p.C**2),1.,places=7)
        self.assertEqual(float(p.planck_b_nu(1e20,1.)),0.)

    def test_spectral_slope(self):
        self.assertAlmostEqual(p.spectral_index(10),.378,places=2)
        self.assertAlmostEqual(p.spectral_index(10,beta=1.7)-p.spectral_index(10),1.7)

    def test_conditional_limits(self):
        horizon = p.detection_horizon(235,10,10)
        self.assertAlmostEqual(float(p.flux_ujy(235,10,d_pc=horizon)),50.)
        self.assertGreater(float(p.illustrative_fraction_limit(1e40)),1.)

    def test_invalid_input(self):
        for args in [(0,10),(100,0)]:
            with self.assertRaises(ValueError):
                p.flux_ujy(*args)
        with self.assertRaises(ValueError):
            p.radiator_area(2.,T_bg=2.725)
        with self.assertRaises(ValueError):
            p.waste_power(eta=.5)


class DynamicsTests(unittest.TestCase):
    def test_conflict_resolution_independent_of_order(self):
        local = np.zeros((3,3),dtype=np.int8)
        local.flat[1] = 1  # spontaneous birth competing with colonization
        a = resolve_proposals(local,[1,2,1],[2,3,2])
        b = resolve_proposals(local,[1,2,1][::-1],[2,3,2][::-1])
        np.testing.assert_array_equal(a,b)
        self.assertEqual(a.flat[1],2)
        self.assertEqual(a.flat[2],3)

    def test_birth_before_parent_death_and_no_same_step_cascade(self):
        sim = GalacticEcosystem(grid_size=5,p_birth=0,p_expand=1,s_travel_cost=0,
                                decay_extrovert=1,decay_ruin=0,seed=1)
        sim.grid[2,2]=2
        sim.step()
        self.assertEqual(np.count_nonzero(sim.grid==2),1)
        self.assertEqual(sim.grid[2,2],3)

    def test_mutation_priority(self):
        sim = GalacticEcosystem(grid_size=3,p_mutate=1,decay_introvert=1)
        sim.grid[:]=1
        sim.step()
        self.assertTrue(np.all(sim.grid==2))

    def test_conflict_kills_both_expansionists(self):
        sim = GalacticEcosystem(grid_size=3,p_expand=1,s_travel_cost=0,s_purge=1,decay_extrovert=0)
        sim.grid[:]=2
        sim.step()
        self.assertTrue(np.all(sim.grid==3))

    def test_reproducibility_and_state_conservation(self):
        a = GalacticEcosystem(grid_size=8,p_birth=.1,seed=42).run(100)
        b = GalacticEcosystem(grid_size=8,p_birth=.1,seed=42).run(100)
        np.testing.assert_array_equal(a.grid,b.grid)
        np.testing.assert_allclose(np.sum(list(a.history.values()),axis=0),1)

    def test_two_state_stationary_balance(self):
        # Independent two-state Markov-chain prediction, no expansion/mutation.
        sim = GalacticEcosystem(grid_size=60,p_birth=.1,p_mutate=0,decay_introvert=.2,seed=4).run(1200)
        self.assertLess(abs(np.mean(sim.history['introvert'][200:])-1/3),.005)

    def test_invalid_probabilities(self):
        with self.assertRaises(ValueError):
            GalacticEcosystem(p_expand=1.1)
        with self.assertRaises(ValueError):
            GalacticEcosystem(grid_size=2)


if __name__ == '__main__':
    unittest.main(verbosity=2)
