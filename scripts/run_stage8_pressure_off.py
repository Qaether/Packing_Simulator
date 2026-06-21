#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from qaether_sim.analysis import run_stage8_pressure_off


def main():
    parser = argparse.ArgumentParser(description="Run Stage 8 pressure-off dynamics from Phase A states.")
    parser.add_argument("--phase-a-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--snapshot-stride", type=int, default=50)
    args = parser.parse_args()

    summary = run_stage8_pressure_off(
        args.phase_a_dir,
        args.out_dir,
        steps=args.steps,
        snapshot_stride=args.snapshot_stride,
    )
    print(summary)


if __name__ == "__main__":
    main()
