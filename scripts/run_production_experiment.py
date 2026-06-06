#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from qaether_sim.analysis import run_production_pipeline


def main():
    parser = argparse.ArgumentParser(description="Run Qaether experiment pipeline through Stage 7 or Stage 10.")
    parser.add_argument("--out-dir", default="results_production", help="Output directory")
    parser.add_argument("--n", type=int, default=64, help="Number of unit spaces")
    parser.add_argument("--seeds", type=str, default="0", help="Comma-separated seed integers, e.g., '0,1,2'")
    parser.add_argument("--smoke", action="store_true", help="Run in lightweight smoke mode for verification")
    parser.add_argument(
        "--max-stage",
        type=int,
        choices=(7, 10),
        default=10,
        help="Stop after the geometry atlas (7) or include dynamics pilots (10)",
    )
    args = parser.parse_args()

    seeds = [int(s.strip()) for s in args.seeds.split(",") if s.strip()]

    print(f"Starting Qaether Pipeline (Stages 0-{args.max_stage})")
    print(f"Config: N={args.n}, Seeds={seeds}, Mode={'Smoke (lightweight)' if args.smoke else 'Full Sweep'}, OutDir={args.out_dir}")

    summary = run_production_pipeline(
        args.out_dir,
        n=args.n,
        seeds=seeds,
        smoke=args.smoke,
        max_stage=args.max_stage,
    )
    print(f"Execution complete. Summary: {summary}")


if __name__ == "__main__":
    main()
