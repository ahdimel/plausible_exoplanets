# Kepler-dichotomy study reproduction (paper/dichotomy_section.md).
#
# reproduce-dichotomy re-derives every number and figure from the
# simulation caches committed in results/ (the sweep scripts are
# resumable and skip cells whose JSON already exists, so with intact
# caches the sweeps are no-ops and this takes minutes). Deleting
# results/grid* / results/robustness turns the same target into a full
# ~10 h re-simulation — but that also needs the DR25 snapshots in
# data/, which are NOT committed (see data/PROVENANCE.md):
# dichotomy-from-scratch fetches them first.

PY := .venv/bin/python

.PHONY: reproduce-dichotomy dichotomy-from-scratch dr25-data test

# Order matters from a cold start: grid_topup and robustness --run pick
# their cells from results/grid_inference.json, so the inference runs
# once at M=8 before the top-up and once more at final depth after it.
reproduce-dichotomy:
	$(PY) analysis/grid_sweep.py --fine --workers 4
	$(PY) analysis/grid_surface.py
	$(PY) analysis/grid_inference.py
	$(PY) analysis/grid_topup.py --workers 4
	$(PY) analysis/grid_inference.py
	$(PY) analysis/mixture_sweep.py --workers 4
	$(PY) analysis/mixture_inference.py
	$(PY) analysis/robustness.py --run all --workers 4
	$(PY) analysis/robustness.py --analyze
	$(PY) analysis/paper_figures.py

dr25-data:
	$(PY) -c "import sys; sys.path.insert(0, 'src'); \
	from exoverse.kepler_data import fetch_dr25; fetch_dr25()"

dichotomy-from-scratch: dr25-data reproduce-dichotomy

test:
	$(PY) -m pytest -q
