"""Run optional deterministic checks against an ORIGINAL repository checkout.

Usage: python scripts/audit_upstream.py /path/to/original-checkout
This imports the two original scripts but never executes their main blocks
or edits them. The output proves the traversal-order bug and independently
compares the original Planck routine at the actual manuscript wavelengths.
"""
from pathlib import Path
from unittest.mock import patch
import argparse
import importlib.util
import json
import numpy as np
from physics import planck_b_nu, C


def load(path, name):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('upstream',type=Path)
    args = parser.parse_args()
    old = load(args.upstream/'scripts'/'run_simulation.py','old_sim')
    spectra = load(args.upstream/'scripts'/'make_paper_figures.py','old_spectra')
    cases = []
    for source,target in [((1,1),(1,2)),((1,2),(1,1))]:
        sim = old.GalacticEcosystem(grid_size=4,p_birth=1,p_expand=1,s_travel_cost=0,
                                   decay_extrovert=0,decay_introvert=0,p_mutate=0)
        sim.grid[source]=2
        neighbor_index = sim._get_neighbors(*source).index(target)
        # Every empty site certainly births; source certainly colonizes.
        # It must be irrelevant which site is visited first.
        with patch.object(old.np.random,'rand',return_value=.5), patch.object(old.np.random,'randint',return_value=neighbor_index):
            sim.step()
        cases.append(dict(source=source,target=target,target_state=int(sim.grid[target])))
    wavelengths = np.array([28,92,100,126,183,235,450,850,1300,3000.])
    differences = []
    for T in (10,15,30):
        nu = C/(wavelengths*1e-6)
        differences.append(float(np.max(np.abs(spectra.planck_b_nu(nu,T)/planck_b_nu(nu,T)-1))))
    result = dict(original_empty_target_overwrite_cases=cases,
                  original_planck_max_relative_difference_at_benchmark_wavelengths=max(differences))
    assert cases[0]['target_state']==1 and cases[1]['target_state']==2
    assert max(differences)<1e-12
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
