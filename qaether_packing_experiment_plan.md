# Qaether Unit-Space Compression, Graph, Motif, and Dynamics Experiment Plan

## 0. Core Objective

This experiment tests three questions in sequence.

```text
1. When unit spaces are compressed, what spatial graph appears at each packing ratio phi?

2. In jammed or frustrated states, what structure is revealed by
   primitive triangle/square cycles and T/O motif detection?

3. Starting from those states and from HCP/FCC lattice states, what motion appears
   when the external pressure/compression protocol is removed and only exclusion remains?
   How does that motion change when unit spaces oscillate with a common quaternion frequency
   but different phases, and how does it change when a small pre-pressure-off vibration is added?
```

Interpretation:

- Qaether = vertex or unit-space center
- sphere = effective exclusion proxy for a Qaether unit space
- `positions` are computational embedding coordinates for contact-graph construction, not fundamental particle coordinates in a background space.
- `velocities` are effective relaxation variables, not literal microscopic particle velocities.
- The primary objects are `G=(V,E)`, primitive cycles, T/O motifs, and their survival/transition.
- `ell_Q = 1`
- `R_Q = ell_Q / 2`
- packing fraction:

```text
phi = N * (4*pi*R_Q^3/3) / V_cell
```

Production geometry:

- Use a periodic unit cell by default.
- Do not include physical wall, boundary shell, peeling, or detached-layer experiments in this plan.
- "Pressure off" means turning off the compression or isotropic pressure-control protocol and integrating only exclusion-based dynamics inside the periodic cell.
- HCP/FCC states are constructed lattice benchmarks.
- This is interpreted as embedded graph/motif response of minimal space-units, not particle motion through pre-existing space.

Execution phases:

- Phase A, geometry-only atlas: Stage 0-7, establishing `phi -> G_phi -> C3/C4/T/O`.
- Phase B, pressure-off survival: Stage 8, 10, 11, 12, tracking `G(t)`, `T(t)`, and `O(t)` survival.
- Phase C, phase coupling pilot: Stage 9, using scalar phase as a controlled proxy.
- Phase D, scaling and topology transition: Stage 13-17.

## Stage Status Checklist

Mark a stage complete only when required outputs exist and smoke/validation checks pass.

- [x] Stage 0: freeze definitions, protocols, state/config/HDF5/metadata schema
- [x] Stage 1: implement initial condition generation
- [x] Stage 2: run packing-ratio sweep periodic compression
- [x] Stage 3: build spatial contact graph atlas and classify jamming/frustration
- [x] Stage 4: detect primitive triangle/square cycles
- [x] Stage 5: detect T/O motifs
- [x] Stage 6: build compression graph/cycle/motif atlas
- [x] Stage 7: select representative states and HCP/FCC benchmarks
- [ ] Stage 8: run pressure-off exclusion-only dynamics and time-series analysis
- [ ] Stage 9: run quaternion phase coupling and lambda/phase-seed sweep
- [ ] Stage 10: run pre-pressure-off small-vibration perturbation sweep
- [ ] Stage 11: compare baseline, phase-conditioned, perturbed, HCP, and FCC dynamics
- [ ] Stage 12: analyze dynamics sensitivity by packing ratio
- [ ] Stage 13: run large-N compression graph atlas scaling
- [ ] Stage 14: run large-N pressure-off dynamics scaling
- [ ] Stage 15: analyze O-motif rarity and abundance-like observables
- [ ] Stage 16: analyze T/O motif survival and transitions
- [ ] Stage 17: analyze graph topology transitions
- [ ] Stage 18: generate final tables, figures, and report
- [ ] Stage 19: write final interpretation and separate future boundary/finite-wall plans

Progress rules:

- Keep the Korean and English plan files synchronized.
- Do not mark a stage complete without outputs and validation.
- Store random-compressed, jammed/frustrated, HCP, FCC, perturbed, and phase-conditioned results under separate labels.

Phase A execution record:

- 2026-06-06: completed Stages 0-7 with `N=64`, seeds `0,1,2`, and 15 `phi_target` values
- Output directory: `results_phase_a_N64`
- Compression relaxation: `relax_dt=0.2` separated from the dynamics timestep, up to 3000 steps with early convergence
- No Stage 8+ dynamics, phase, or perturbation outputs were generated
- 2026-06-06 scaling extension: repeated Stages 0-7 with `N=128`, seeds `0,1,2`, and the same 15 `phi_target` values
- Scaling output directory: `results_phase_a_N128`
- Topology outputs: `structure_topology_by_phi.csv` and `structure_topology_summary_by_phi.csv` in each result directory
- Cross-size comparison: `structure_topology_scaling_N64_N128.csv` and `structure_topology_scaling_N64_N128.md`

Correctness gates:

- A zero-pressure, zero-energy, contact-free state must be classified as `flowing`, never jammed or frustrated.
- Wrapped coordinate differences must not be used directly for MSD; use minimum-image increments or an unwrapped trajectory.
- FCC/HCP benchmarks must contain complete periodic unit cells. For the current four-site orthorhombic cells, `N` must be divisible by four and the replicated cell count must equal `N / 4` exactly.
- A lattice may be labeled `toto_valid=true` only after both the hard-sphere overlap check and the same T/O detector used for experimental states pass.
- Dynamic edges, cycles, and motifs at one time point must all be evaluated from the same hysteretic contact graph.
- If the initial edge or motif set is empty, its survival quotient is undefined and must be stored as `NaN` together with the initial count. It must not be reported as zero or one.
- A smoke run verifies execution and invariants only. It does not authorize production-scale or publication-level interpretation.

## 1. Recommended Python Stack

- `numpy`, `scipy`: positions, velocities, distances, `cKDTree`, relaxation
- `numba`: accelerated force loops, periodic minimum-image distances, contact updates
- `networkx`: spatial contact graphs, cycles, cliques, motif graph analysis
- `h5py`: snapshots, trajectories, per-frame observables
- `pandas`: summary tables
- `matplotlib`, `plotly`: static and interactive figures

## 2. Proposed Code Structure

```text
qaether_sim/
    __init__.py
    config.py
    state.py
    initial_conditions.py
    compression.py
    forces.py
    contact_graph.py
    cycles.py
    motifs_T.py
    motifs_O.py
    graph_atlas.py
    jamming.py
    lattice_benchmarks.py
    bulk_dynamics.py
    phase.py
    perturbation.py
    survival.py
    analysis.py
    visualization.py
    run_compression_atlas.py
    run_pressure_off_dynamics.py
    run_phase_sweep.py
    run_perturbation_sweep.py
```

Core state:

```python
positions: np.ndarray        # computational embedding coordinates, not fundamental background space
velocities: np.ndarray       # effective relaxation variable, not literal particle velocity
radius: float
ell_Q: float
box_matrix: np.ndarray       # periodic unit cell
graph_edges: np.ndarray      # primary spatial adjacency/contact structure
q: np.ndarray | None         # optional quaternion state
theta: np.ndarray | None     # optional oscillator phase
motif_data: dict             # cycles, T/O motifs, motif incidence
metadata: dict
```

## 3. Stage 0: Definitions, Protocols, and Schema

Purpose:

- Define unit space, contact, graph, cycle, T/O motif, jamming, frustration, and pressure-off dynamics.
- Treat quaternion phase oscillation as a controlled modulation, not as arbitrary full q-dynamics.
- Ensure every stage uses the same HDF5 snapshot and metadata schema.

Definitions:

- `contact`: minimum-image distance `r_ij <= ell_Q * (1 + epsilon_contact)`
- `contact_hysteresis`: edge birth at `r_ij <= ell_Q * (1 + eps_on)`, edge death at `r_ij >= ell_Q * (1 + eps_off)`, with `eps_off > eps_on`
- `spatial graph`: contact graph at each snapshot, `G_phi = (V, E_phi)`
- `jamming`: mechanical constraint stability, measured by overlap energy, displacement, coordination, force-balance residual, and pressure response.
- `frustration`: failed further compression/rearrangement caused by incompatible local constraints, measured by graph-transition failure and motif incompatibility.
- `pressure-off`: remove compression or external pressure-control terms and integrate pure exclusion dynamics
- If a state is strictly non-overlapping and velocity-free, pure exclusion pressure-off dynamics is trivial; Stage 8 is residual-stress relaxation when overlaps exist and perturbation-response dynamics otherwise.

Required metadata:

- `N`, `seed`, `phi`, `box_matrix`, `compression_protocol`
- `phi_target`, `phi_achieved`, `max_overlap`, `mean_overlap`, `overlap_energy_per_vertex`, `fraction_overlapping_pairs`
- `force_model`, `epsilon_contact`, `epsilon_cycle`, `epsilon_motif`
- `epsilon_contact_on`, `epsilon_contact_off`, `contact_rule`
- `phase_enabled`, `lambda_phase`, `phase_seed`, `omega_Q`
- `perturbation_enabled`, `perturbation_amplitude`, `perturbation_seed`

Outputs:

- `config_schema.json`
- `metadata_schema.json`
- `protocol_definitions.md`
- HDF5 snapshot read/write smoke test

## 4. Stage 1: Initial Condition Generation

Initial condition families:

- random gas
- random gas plus annealing
- weak HCP/FCC seed
- ideal HCP
- ideal FCC
- dense HCP/FCC benchmark near `phi = 0.74`

Purpose:

- Separate graphs that arise through compression from graphs imposed by constructed lattices.

Outputs:

- `initial_state_seed_XXXX.h5`
- `hcp_state_phi_XXXX.h5`
- `fcc_state_phi_XXXX.h5`

## 5. Stage 2: Packing-Ratio Sweep Periodic Compression

Purpose:

- Compress unit spaces and generate snapshots for observing the spatial graph at each `phi`.

Target packing ratios:

```python
phi_target_list = [
    0.20, 0.30, 0.40, 0.50, 0.55, 0.58, 0.60,
    0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.735, 0.740,
]
```

Baseline force model:

```text
U_core(r) = 0.5*k*(ell_Q - r)^2, if r < ell_Q
          = 0, otherwise
```

Outputs:

- `snapshot_phi_XXXX.h5`
- `energy_curve.csv`
- `pressure_estimator_curve.csv`
- `coordination_curve.csv`
- `overlap_diagnostics_by_phi.csv`
- `phi_target_vs_achieved.csv`

Acceptance:

- `phi_target` is the packing ratio attempted by the external protocol.
- `phi_achieved` is the recorded snapshot packing ratio.
- Above `phi ~ 0.64`, accept a state as hard-sphere-like only if overlap metrics are below tolerance; otherwise label it `soft-overcompressed/frustrated`.

## 6. Stage 3: Spatial Contact Graph Atlas and Jamming/Frustration Classification

Purpose:

- Record what spatial graph appears at each packing ratio.
- Select states that are actually jammed or frustrated, not merely dense.

Graph observables:

- `E_phi`
- degree distribution
- connected component size
- clustering coefficient
- graph distance distribution
- contact persistence between neighboring `phi`
- graph edit distance between consecutive `phi`

Jamming/frustration metrics:

```text
Delta MSD < epsilon_msd
Delta Z < epsilon_Z
Delta E < epsilon_E
P_bulk > P_jam or dP/dphi sharply increases
graph_edit_rate(phi) -> plateau
Z_mean
Z_rattler_removed
force_balance_residual
bulk_modulus_proxy dP/dphi
failed_rearrangement_count
```

Outputs:

- `graph_edges_phi_XXXX.csv`
- `graph_summary_by_phi.csv`
- `graph_transition_by_phi.csv`
- `jammed_state_candidates.csv`
- `frustrated_state_candidates.csv`
- `state_selection_report.md`

## 7. Stage 4: Primitive Triangle/Square Cycle Detection

Purpose:

- Identify the basic face-like structure of the spatial graph.
- Before geometric validation, unwrap each candidate into a local minimum-image patch relative to one anchor vertex.

Triangle:

- induced `K3`
- non-collinear
- all edge lengths near `ell_Q`

Square:

- chordless `i-j-k-l-i`
- no diagonal contacts
- near-planar
- edge and diagonal lengths within tolerance

Outputs:

- `primitive_triangles_phi_XXXX.csv`
- `primitive_squares_phi_XXXX.csv`
- `cycle_summary_by_phi.csv`

## 8. Stage 5: T/O Motif Detection

T-motif conditions:

- four vertices induce `K4`
- six edges near `ell_Q`
- tetrahedron volume `> epsilon_volume`
- four primitive triangle faces
- no internal chordless square

O-motif conditions:

- six vertices
- internal contact edge count `12`
- exactly three disjoint opposite non-edges
- opposite pair midpoints share one center
- three opposite axes are nearly orthogonal
- three primitive square cycles
- eight primitive triangle cycles
- the three primitive square cycles cover all 12 O-edges exactly once
- each O-edge is incident to exactly two O-triangles

Outputs:

- `T_motifs_phi_XXXX.csv`
- `O_motifs_phi_XXXX.csv`
- `T_density_by_phi.csv`
- `O_density_by_phi.csv`
- `O_cluster_distribution.csv`
- `P_O_by_phi.csv`

## 9. Stage 6: Compression Graph/Cycle/Motif Atlas

Purpose:

- Combine spatial graphs, primitive cycles, and T/O motifs from the compression path into one atlas.

Questions:

- How does graph topology change with `phi`?
- Do triangle/square/T/O structures grow or collapse near jamming/frustration?
- Which motif signatures distinguish random-compressed graphs from HCP/FCC graphs?

Outputs:

- `compression_graph_motif_atlas.csv`
- `cycle_motif_by_phi_summary.csv`
- `graph_motif_atlas_report.md`
- summary figures

## 10. Stage 7: Representative States and HCP/FCC Benchmarks

Selected states:

- random typical state at selected `phi`
- jammed state
- frustrated state
- T-rich/O-poor state
- O-positive state
- ideal HCP close-packing benchmark
- ideal FCC close-packing benchmark
- TOTO-labeled HCP/FCC-compatible constructed benchmark

HCP/FCC benchmark rules:

- Separate ideal close-packing benchmarks from target-phi constructed benchmarks.
- Generate with the same `N` and periodic cell convention, while recording `phi_target` and `phi_achieved/contact scale` separately.
- Reject incompatible `N` rather than truncating a periodic lattice or shrinking a partially filled box.
- Run the same contact graph, primitive cycle, and T/O motif pipeline.

Outputs:

- `selected_states_manifest.csv`
- `selected_state_*.h5`
- `hcp_fcc_graph_summary.csv`
- `hcp_fcc_motif_summary.csv`

## 11. Stage 8: Pressure-Off Exclusion-Only Dynamics and Time Series

Purpose:

- Starting from selected states and HCP/FCC states, remove external pressure/compression and observe motion driven only by exclusion.
- The central observable is graph/motif response: contact survival, graph edit distance, cycle/motif birth-death, and T/O survival.
- Label a run `A. Residual-stress relaxation` only when initial overlap/stress is present.
- Label a run `B. Perturbation-response` only when the initial state is hard-sphere-like and an explicit perturbation was applied.
- Label an unforced, non-overlapping, unperturbed run `none`; scalar phase assignment alone is not a perturbation when no exclusion force is active.

Baseline dynamics:

- fixed periodic cell
- external compression/pressure term off
- pure exclusion force only
- overdamped and weakly damped dynamics kept under separate protocol labels

Primary observables:

- contact survival `S_E(t)`
- graph edit distance from `t0`
- motif/cycle birth-death events
- primitive triangle/square count over time
- `N_T(t)`, `N_O(t)`
- `S_T(t)`, `S_O(t)`
- graph transition events

Secondary observables:

- displacement MSD
- effective velocity/relaxation energy

Outputs:

- `pressure_off_trajectory_*.h5`
- `pressure_off_summary.csv`
- `dynamic_graph_timeseries.csv`
- `dynamic_motif_timeseries.csv`

## 12. Stage 9: Scalar Phase Proxy Coupling and Lambda/Phase-Seed Sweep

Purpose:

- Test how movement, graph transition, and motif survival change when all unit spaces share one intrinsic frequency `omega_Q` but have different initial phases.
- This is not full SU(2) quaternion dynamics; it is a scalar phase proxy for testing whether internal phase mismatch can bias exclusion stiffness, graph transition, or O-motif survival.

Minimal model:

```text
theta_i(t) = omega_Q * t + theta_i(0)
theta_i(0) ~ Uniform(0, 2*pi)
k_ij(t) = k_core * [1 + lambda * cos(theta_i(t) - theta_j(t))]
0 <= lambda < 1
```

Experiment:

- selected random-compressed states
- selected jammed/frustrated states
- HCP/FCC states
- `lambda = 0.0, 0.1, 0.3, 0.5`
- multiple phase seeds

Metrics:

- `P(O > 0 | phi, lambda)`
- `E[N_O/N | phi, lambda]`
- `S_E(t | lambda)`
- `S_T(t | lambda)`
- `S_O(t | lambda)`
- phase coherence of all contacts
- phase coherence of O-motif edges

Outputs:

- `phase_model_spec.md`
- `phase_sweep_summary.csv`
- `phase_sweep_by_lambda.csv`
- `phase_conditioned_o_rarity.csv`
- `phase_conditioned_survival.csv`

Follow-up pilot:

```text
Stage 9b: SU(2) relative-frame coupling pilot
h_ij = q_i^{-1} q_j
k_ij = k_core * [1 + lambda * f(h_ij)]
Start with f(h_ij) = 0.5 * Re Tr(h_ij) or another simple similarity function.
```

## 13. Stage 10: Pre-Pressure-Off Small-Vibration Perturbation Sweep

Purpose:

- Add small vibration or position perturbation before pressure-off dynamics and repeat the experiment.

Perturbation:

```text
x_i <- x_i + A * u_i
u_i: random unit vector or selected lattice vibration mode
A / ell_Q = 1e-4, 1e-3, 1e-2
```

Protocols:

- random micro-jitter
- acoustic-like long-wavelength mode
- HCP/FCC phonon-like mode
- phase off / phase on, stored under separate labels

Outputs:

- `perturbed_initial_state_*.h5`
- `perturbation_metadata.csv`
- `perturbed_pressure_off_summary.csv`
- `perturbed_dynamic_motif_timeseries.csv`

## 14. Stage 11: Baseline/Phase/Perturbed/HCP/FCC Dynamics Comparison

Purpose:

- Summarize differences across dynamics and initial-structure conditions.

Groups:

- baseline pressure-off
- phase-conditioned pressure-off
- pre-vibrated pressure-off
- phase-conditioned plus pre-vibrated pressure-off
- random-compressed
- HCP
- FCC

Metrics:

- degree distribution distance
- cycle distribution distance
- T/O motif density difference
- graph edit distance to HCP/FCC contact graph
- dynamics survival difference

Outputs:

- `dynamics_comparison_summary.csv`
- `random_vs_lattice_graph_comparison.csv`
- `random_vs_lattice_dynamics_comparison.csv`
- `dynamics_comparison_report.md`

## 15. Stage 12: Dynamics Sensitivity by Packing Ratio

Purpose:

- Test how the same dynamics protocol changes with `phi`.

Metrics:

- `MSD_final(phi)`
- `S_E_final(phi)`
- `S_T_final(phi)`
- `S_O_final(phi)`
- graph transition count versus `phi`

Outputs:

- `phi_dynamics_sensitivity.csv`
- `phi_dynamics_sensitivity_report.md`

## 16. Stage 13: Large-N Compression Graph Atlas Scaling

Purpose:

- Test whether the compression graph atlas is stable with system size.

Smoke matrix:

```text
N = 64, 128
seeds = 3
phi_targets = [0.30, 0.50, 0.64]
states = [random_compressed, HCP, FCC]
dynamics = [baseline, phase, perturbed]
```

Scaling matrix:

```text
Pilot: N = 256, seeds = 10
Production: N = 512, seeds = 20..50
After optimized O-detector: N = 1024, seeds = 10..30
```

Do not use brute-force `N^6` O-motif search. Generate candidates from local contact neighborhoods or 2-hop candidate pools.

Outputs:

- `scaling_graph_atlas_summary.csv`
- `scaling_graph_transition_summary.csv`

## 17. Stage 14: Large-N Pressure-Off Dynamics Scaling

Purpose:

- Test whether pressure-off motion and motif survival persist with `N`.

Outputs:

- `scaling_pressure_off_summary.csv`
- `scaling_motif_survival_summary.csv`

## 18. Stage 15: O-Motif Rarity and Abundance-Like Observables

Purpose:

- Quantify whether O-motifs are rare and how occurrence changes under phase, perturbation, HCP/FCC, and random-compressed protocols.

Metrics:

- `P(O > 0 | phi, N, protocol)`
- `E[N_O/N | phi, N, protocol]`
- `Var(N_O/N | phi, N, protocol)`
- O cluster size distribution
- O survival after pressure-off dynamics

Outputs:

- `o_rarity_by_protocol.csv`
- `o_abundance_observable_report.md`

## 19. Stage 16: T/O Motif Survival and Transitions

Purpose:

- Record how T/O motifs are preserved, destroyed, born, or transformed during pressure-off dynamics.

Metrics:

- motif birth/death events
- T-to-O adjacency changes
- motif cluster split/merge
- motif lifetime distribution

Outputs:

- `motif_lifetime_distribution.csv`
- `motif_transition_events.csv`

## 20. Stage 17: Graph Topology Transition Analysis

Purpose:

- Understand how unit-space motion appears as graph topology change.

Metrics:

- edge birth/death rate
- graph edit distance
- component stability
- cycle birth/death rate
- motif transition correlated with graph transition

Outputs:

- `graph_topology_transition_summary.csv`
- `graph_transition_events.csv`

## 21. Stage 18: Final Tables, Figures, and Report

Final tables:

- `compression_graph_atlas_summary.csv`
- `cycle_motif_by_phi_summary.csv`
- `selected_state_summary.csv`
- `pressure_off_dynamics_summary.csv`
- `phase_perturbation_comparison_summary.csv`
- `random_vs_lattice_summary.csv`
- `scaling_summary.csv`

Final figures:

- graph topology versus `phi`
- triangle/square cycle count versus `phi`
- T/O density versus `phi`
- pressure-off MSD/energy/contact survival
- phase sweep comparison
- perturbation comparison
- HCP/FCC versus random-compressed comparison

## 22. Stage 19: Final Interpretation and Follow-Up Separation

Interpretation questions:

- What spatial graph family appears at each `phi` under compression?
- What are the cycle/T/O motif signatures of jammed or frustrated states?
- Are HCP/FCC and random-compressed graphs in the same structural family or different families?
- After pressure is removed, how much of the graph and motif structure survives under exclusion-only dynamics?
- Do quaternion phase and pre-pressure-off small vibration significantly change motion and motif survival?

Follow-up:

- Boundary, finite-wall, and peeling experiments remain separate from this plan.
