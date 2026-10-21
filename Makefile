PYTHON ?= python3
S = scripts

.PHONY: all test figures simulations robustness paper checksums clean

all: figures simulations robustness checksums

test:
	$(PYTHON) -m unittest discover -s tests -v

figures:
	cd $(S) && $(PYTHON) make_paper_figures.py > ../results/figure_generation.txt
	cd $(S) && $(PYTHON) observational_constraints.py > ../results/observational_constraints.txt

# Long: 40 x 20,000-step runs plus 800 lineage trials.
simulations:
	cd $(S) && $(PYTHON) run_simulation.py
	cd $(S) && $(PYTHON) compare_updates.py

# Long: lattice-size, parameter-sweep and inward-vulnerability experiments.
robustness:
	cd $(S) && $(PYTHON) robustness.py all

paper:
	cd paper && latexmk -pdf -interaction=nonstopmode main.tex

checksums:
	git ls-files | grep -v '^SHA256SUMS.txt$$' | xargs sha256sum > SHA256SUMS.txt

clean:
	cd paper && latexmk -C
