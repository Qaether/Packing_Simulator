# Stage 8 Pressure-off Survival Report

## Protocol

- Source states: `results_phase_a_N64`, `results_phase_a_N128`
- System sizes: N=64 and N=128
- Seeds: 0, 1, 2 where a selected state exists
- Integration: overdamped exclusion-only dynamics
- Fixed periodic cell, no phase coupling, no explicit perturbation
- `dt=0.01`, 1000 steps, snapshot stride 50
- Contact tracking: one hysteretic graph per frame

## Results

- N=64: 13 runs, 273 recorded frames, 13 HDF5 trajectories
- N=128: 12 runs, 252 recorded frames, 12 HDF5 trajectories
- Every initial contact survived: final `S_E=1.0` for every run.
- Every initially present T-motif survived: final `S_T=1.0`.
- Every initially present O-motif survived: final `S_O=1.0`.
- N=128 seed 0 random gained two contact edges without losing an initial edge.
- N=128 seed 0 jammed contained two initial O-motifs and retained both.
- All other runs had final graph edit distance zero.

## Interpretation

Turning off the compression protocol while holding the periodic cell fixed does not
produce substantial rearrangement. The selected stressed states carry overlap energy,
but their residual forces are already nearly balanced. Hard random and lattice states
are force-free and therefore exactly static under pure exclusion dynamics.

This baseline establishes that spontaneous pressure-off decay is absent on the tested
time interval. A non-trivial survival response requires the explicit Stage 10
perturbation protocol or a future variable-cell pressure-release protocol.

## Data

- `results_stage8_pressure_off_N64`
- `results_stage8_pressure_off_N128`
- `stage8_pressure_off_scaling_N64_N128.csv`
