#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from qaether_sim.analysis import run_smoke_pipeline


def main():
    parser = argparse.ArgumentParser(description="Run Qaether smoke experiment pipeline.")
    parser.add_argument("--out-dir", default="results_smoke", help="Output directory")
    parser.add_argument("--n", type=int, default=32, help="Number of unit spaces")
    parser.add_argument("--seed", type=int, default=0, help="Random seed")
    args = parser.parse_args()
    summary = run_smoke_pipeline(args.out_dir, n=args.n, seed=args.seed)
    print(summary)


if __name__ == "__main__":
    main()
