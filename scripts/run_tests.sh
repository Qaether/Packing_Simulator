#!/usr/bin/env bash
set -euo pipefail
export PYTHONPYCACHEPREFIX="${PYTHONPYCACHEPREFIX:-.pycache_tmp}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-.mplconfig}"
python3 -m unittest discover -s tests -p 'test_*.py'
python3 scripts/run_smoke_experiment.py --out-dir results_smoke --n 32 --seed 0
