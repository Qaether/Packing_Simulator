# Qaether 단위공간 압축/그래프/motif/dynamics 실험 계획안

## 0. 전체 목적

이번 실험은 세 질문을 순서대로 검증한다.

```text
1. 단위공간을 압축할 때 각 packing ratio phi에서 어떤 공간 그래프가 형성되는가?

2. jamming 또는 frustration 상태의 공간 그래프에서
   삼각/사각 primitive cycle, T-motif, O-motif가 어떤 구조를 보이는가?

3. random-compressed 상태와 HCP/FCC 격자 상태에서 외부 압력/압축을 제거하면,
   배제력만 가진 단위공간들이 어떻게 움직이며,
   quaternion phase oscillation 또는 pressure-off 전 작은 진동 perturbation이
   그 움직임과 graph/motif survival을 어떻게 바꾸는가?
```

해석 기준:

- Qaether = vertex 또는 단위공간 중심
- sphere = Qaether 단위공간의 effective exclusion proxy
- `positions`는 배경공간 속 입자 좌표가 아니라 contact graph를 만들기 위한 computational embedding proxy이다.
- `velocities`는 literal particle velocity가 아니라 relaxation protocol의 effective update variable이다.
- 따라서 근본 관측 대상은 `rho_i` 자체가 아니라 `G=(V,E)`, primitive cycle, T/O motif 및 그 survival/transition이다.
- `ell_Q = 1`
- `R_Q = ell_Q / 2`
- packing fraction:

```text
phi = N * (4*pi*R_Q^3/3) / V_cell
```

Production geometry:

- periodic unit cell을 기본으로 한다.
- 물리적 wall, boundary shell, peeling, detached-layer 실험은 포함하지 않는다.
- "압력 제거"는 외부 압축 또는 isotropic pressure-control protocol을 끄고, periodic cell 안에서 배제력 기반 dynamics만 적분하는 것으로 정의한다.
- HCP/FCC는 비교용 constructed lattice benchmark로 사용한다.
- 이 실험은 "공간 안 입자의 운동"이 아니라, 최소 단위공간의 embedded graph/motif response 실험으로 해석한다.

실행 phase:

- Phase A, geometry-only atlas: Stage 0-7. Quaternion, perturbation, dynamics를 끄고 `phi -> G_phi -> C3/C4/T/O`를 먼저 확정한다.
- Phase B, pressure-off survival: Stage 8, 10, 11, 12. `G(t)`, `T(t)`, `O(t)` survival과 transition을 본다.
- Phase C, phase coupling pilot: Stage 9. scalar phase proxy가 graph/motif survival에 주는 영향을 본다.
- Phase D, scaling/transition: Stage 13-17. `N` 증가와 graph topology transition의 안정성을 확인한다.

## Stage 진행 체크리스트

각 stage가 끝날 때마다 required output과 smoke/validation check가 있을 때만 완료 처리한다.

- [x] Stage 0: 정의, protocol, state/config/HDF5/metadata schema 고정
- [x] Stage 1: 초기 조건 생성기 구현
- [x] Stage 2: packing-ratio sweep periodic compression 실행
- [x] Stage 3: spatial contact graph atlas 및 jamming/frustration 판정
- [x] Stage 4: primitive triangle/square cycle 검출
- [x] Stage 5: T/O motif 검출
- [x] Stage 6: compression graph/cycle/motif atlas 생성
- [x] Stage 7: representative state 및 HCP/FCC benchmark 선정
- [ ] Stage 8: pressure-off exclusion-only dynamics 및 time series 분석
- [ ] Stage 9: quaternion phase coupling 및 lambda/phase-seed sweep
- [ ] Stage 10: pressure-off 전 small-vibration perturbation sweep
- [ ] Stage 11: baseline/phase/perturbed/HCP/FCC dynamics 비교
- [ ] Stage 12: packing ratio별 dynamics sensitivity 분석
- [ ] Stage 13: large-N compression graph atlas scaling
- [ ] Stage 14: large-N pressure-off dynamics scaling
- [ ] Stage 15: O-motif rarity 및 abundance-like observable 분석
- [ ] Stage 16: T/O motif survival 및 transition 분석
- [ ] Stage 17: graph topology transition 분석
- [ ] Stage 18: final tables/figures/report 생성
- [ ] Stage 19: 최종 해석과 다음 boundary/finite-wall 후속 계획 분리

진행 규칙:

- 한글/영문 plan 파일을 항상 같이 업데이트한다.
- 완료 결과 없이 stage를 완료 처리하지 않는다.
- random-compressed, jammed/frustrated, HCP, FCC, perturbed, phase-conditioned 결과는 반드시 별도 라벨로 저장한다.

Phase A 실행 기록:

- 2026-06-06: `N=64`, seeds `0,1,2`, 15개 `phi_target`으로 Stage 0-7 완료
- 결과 디렉터리: `results_phase_a_N64`
- compression relaxation: dynamics timestep과 분리된 `relax_dt=0.2`, 최대 3000 step, 수렴 시 조기 종료
- Stage 8 이후 dynamics/phase/perturbation 산출물은 생성하지 않음

정합성 통과 조건:

- 압력, 에너지, 접촉이 모두 0인 상태는 반드시 `flowing`이어야 하며 jammed/frustrated로 판정하면 안 된다.
- wrapped 좌표를 직접 빼서 MSD를 계산하지 않는다. minimum-image increment 또는 unwrapped trajectory를 사용한다.
- FCC/HCP benchmark는 완전한 periodic unit cell만 복제한다. 현재 4-site orthorhombic cell에서는 `N`이 4의 배수이고 복제 cell 수가 정확히 `N / 4`여야 한다.
- lattice의 `toto_valid=true`는 hard-sphere overlap 검사와 실험 상태에 사용하는 동일한 T/O detector를 모두 통과한 뒤에만 부여한다.
- 동일 시점의 dynamic edge, cycle, motif는 모두 동일한 hysteretic contact graph에서 계산한다.
- 초기 edge/motif 집합이 비어 있으면 survival quotient는 정의되지 않으므로 초기 개수와 함께 `NaN`으로 저장한다. 0 또는 1로 기록하지 않는다.
- smoke run은 실행 가능성과 invariant만 검증하며 production 또는 논문 수준 해석을 승인하지 않는다.

## 1. 추천 Python 라이브러리

- `numpy`, `scipy`: positions, velocities, distances, `cKDTree`, relaxation
- `numba`: force loop, periodic minimum-image distance, contact update 가속
- `networkx`: spatial contact graph, cycle, clique, motif graph 분석
- `h5py`: snapshots, trajectories, per-frame observables
- `pandas`: summary table
- `matplotlib`, `plotly`: static/interactive figures

## 2. 추천 코드 구조

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

핵심 state:

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

## 3. Stage 0: 정의, protocol, schema 고정

목적:

- 단위공간, contact, graph, cycle, T/O motif, jamming, frustration, pressure-off dynamics를 명확히 정의한다.
- quaternion phase oscillation을 arbitrary full q-dynamics가 아니라 controlled modulation으로 둔다.
- 모든 stage가 같은 HDF5 snapshot과 metadata schema를 쓰게 한다.

고정 정의:

- `contact`: minimum-image distance 기준 `r_ij <= ell_Q * (1 + epsilon_contact)`
- `contact_hysteresis`: edge 생성은 `r_ij <= ell_Q * (1 + eps_on)`, edge 제거는 `r_ij >= ell_Q * (1 + eps_off)`, `eps_off > eps_on`
- `spatial graph`: 각 snapshot의 contact graph `G_phi = (V, E_phi)`
- `jamming`: mechanical constraint stability. relaxation 후 overlap energy, displacement, coordination, force balance, pressure estimator가 plateau 또는 sharp growth를 보이는 상태
- `frustration`: incompatible local constraints 때문에 further compression/rearrangement가 실패하는 상태. graph transition 실패와 motif incompatibility를 중심으로 판정
- `pressure-off`: compression 또는 external pressure-control term을 제거하고 pure exclusion dynamics만 적분
- strictly non-overlapping이고 velocity-free인 상태에서는 pure exclusion pressure-off dynamics가 trivial하다. 따라서 Stage 8은 overlap이 있으면 residual-stress relaxation, overlap이 거의 없으면 small displacement/velocity perturbation response로 해석한다.

필수 metadata:

- `N`, `seed`, `phi`, `box_matrix`, `compression_protocol`
- `phi_target`, `phi_achieved`, `max_overlap`, `mean_overlap`, `overlap_energy_per_vertex`, `fraction_overlapping_pairs`
- `force_model`, `epsilon_contact`, `epsilon_cycle`, `epsilon_motif`
- `epsilon_contact_on`, `epsilon_contact_off`, `contact_rule`
- `phase_enabled`, `lambda_phase`, `phase_seed`, `omega_Q`
- `perturbation_enabled`, `perturbation_amplitude`, `perturbation_seed`

산출물:

- `config_schema.json`
- `metadata_schema.json`
- `protocol_definitions.md`
- HDF5 snapshot read/write smoke test

## 4. Stage 1: 초기 조건 생성

초기 조건 family:

- random gas
- random gas + annealing
- weak HCP/FCC seed
- ideal HCP
- ideal FCC
- dense HCP/FCC benchmark near `phi = 0.74`

목적:

- compression path에서 자연스럽게 생기는 graph와 constructed lattice의 graph를 분리해서 비교한다.

산출물:

- `initial_state_seed_XXXX.h5`
- `hcp_state_phi_XXXX.h5`
- `fcc_state_phi_XXXX.h5`

## 5. Stage 2: Packing-ratio sweep periodic compression

목적:

- 단위공간을 압축하면서 각 `phi`에서 spatial graph를 관측할 snapshot을 생성한다.

압축 target:

```python
phi_target_list = [
    0.20, 0.30, 0.40, 0.50, 0.55, 0.58, 0.60,
    0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.735, 0.740,
]
```

force model baseline:

```text
U_core(r) = 0.5*k*(ell_Q - r)^2, if r < ell_Q
          = 0, otherwise
```

산출물:

- `snapshot_phi_XXXX.h5`
- `energy_curve.csv`
- `pressure_estimator_curve.csv`
- `coordination_curve.csv`
- `overlap_diagnostics_by_phi.csv`
- `phi_target_vs_achieved.csv`

판정:

- `phi_target`은 외부 protocol이 시도한 packing ratio이다.
- `phi_achieved`는 실제 box/snapshot에서 기록된 packing ratio이다.
- `phi > 0.64` 구간은 `max_overlap`, `mean_overlap`, `overlap_energy_per_vertex`가 tolerance 이하일 때만 hard-sphere-like packing으로 인정한다.
- tolerance를 넘으면 achieved hard packing이 아니라 `soft-overcompressed/frustrated` 상태로 라벨링한다.

## 6. Stage 3: Spatial contact graph atlas 및 jamming/frustration 판정

목적:

- 각 packing ratio에서 어떤 공간 그래프가 생기는지 기록한다.
- 단순히 높은 `phi`가 아니라 실제 jammed 또는 frustrated 된 상태를 선별한다.

graph observable:

- `E_phi`
- degree distribution
- connected component size
- clustering coefficient
- graph distance distribution
- contact persistence between neighboring `phi`
- graph edit distance between consecutive `phi`

jamming/frustration metric:

```text
Delta MSD < epsilon_msd
Delta Z < epsilon_Z
Delta E < epsilon_E
P_bulk > P_jam 또는 dP/dphi 급증
graph_edit_rate(phi) -> plateau
Z_mean
Z_rattler_removed
force_balance_residual
bulk_modulus_proxy dP/dphi
failed_rearrangement_count
```

산출물:

- `graph_edges_phi_XXXX.csv`
- `graph_summary_by_phi.csv`
- `graph_transition_by_phi.csv`
- `jammed_state_candidates.csv`
- `frustrated_state_candidates.csv`
- `state_selection_report.md`

## 7. Stage 4: Primitive triangle/square cycle 검출

목적:

- 공간 그래프의 기본 face-like 구조를 파악한다.
- geometric validation 전에 각 candidate vertex를 anchor 기준 local minimum-image patch로 unwrap한다.

Triangle:

- induced `K3`
- non-collinear
- 모든 edge length가 `ell_Q` 근처

Square:

- chordless `i-j-k-l-i`
- diagonal contact 없음
- near-planar
- edge length와 diagonal length가 tolerance 안에 있음

산출물:

- `primitive_triangles_phi_XXXX.csv`
- `primitive_squares_phi_XXXX.csv`
- `cycle_summary_by_phi.csv`

## 8. Stage 5: T/O motif 검출

T-motif 조건:

- four vertices induced `K4`
- six edges near `ell_Q`
- tetrahedron volume `> epsilon_volume`
- four primitive triangle faces
- 내부 chordless square 없음

O-motif 조건:

- six vertices
- internal contact edge count `12`
- exactly three disjoint opposite non-edges
- opposite pair midpoints share one center
- three opposite axes are nearly orthogonal
- three primitive square cycles
- eight primitive triangle cycles
- three primitive square cycles가 12개 O-edge를 정확히 한 번씩 cover
- 각 O-edge가 정확히 두 개 O-triangle에 incident

산출물:

- `T_motifs_phi_XXXX.csv`
- `O_motifs_phi_XXXX.csv`
- `T_density_by_phi.csv`
- `O_density_by_phi.csv`
- `O_cluster_distribution.csv`
- `P_O_by_phi.csv`

## 9. Stage 6: Compression graph/cycle/motif atlas

목적:

- compression 과정의 spatial graph, primitive cycle, T/O motif 내용을 한 atlas로 묶는다.

핵심 질문:

- 각 `phi`에서 graph topology가 어떻게 바뀌는가?
- jamming/frustration 근처에서 triangle/square/T/O 구조가 급증하거나 사라지는가?
- random-compressed graph와 HCP/FCC graph는 어떤 motif signature가 다른가?

산출물:

- `compression_graph_motif_atlas.csv`
- `cycle_motif_by_phi_summary.csv`
- `graph_motif_atlas_report.md`
- summary figures

## 10. Stage 7: Representative state 및 HCP/FCC benchmark 선정

선택 state:

- random typical state at selected `phi`
- jammed state
- frustrated state
- T-rich/O-poor state
- O-positive state
- ideal HCP close-packing benchmark
- ideal FCC close-packing benchmark
- TOTO-labeled HCP/FCC-compatible constructed benchmark

HCP/FCC benchmark 조건:

- ideal close-packing benchmark와 target-phi constructed benchmark를 분리한다.
- 동일 `N`, 동일 periodic cell convention으로 생성하되, `phi_target`과 `phi_achieved/contact scale`을 별도 기록한다.
- periodic lattice를 절단하거나 일부만 채운 box를 축소하지 말고, 호환되지 않는 `N`은 명시적으로 거부한다.
- contact graph, primitive cycle, T/O motif를 같은 pipeline으로 검출한다.

산출물:

- `selected_states_manifest.csv`
- `selected_state_*.h5`
- `hcp_fcc_graph_summary.csv`
- `hcp_fcc_motif_summary.csv`

## 11. Stage 8: Pressure-off exclusion-only dynamics 및 time series

목적:

- 앞서 선택한 state와 HCP/FCC state에서 외부 압축/압력을 제거했을 때 배제력만으로 어떤 움직임이 나오는지 본다.
- 중심 질문은 "얼마나 움직였는가"가 아니라 공간 그래프, cycle, T/O motif가 유지/전환되는가이다.
- 초기 overlap/stress가 있을 때만 `A. Residual-stress relaxation`으로 라벨링한다.
- hard-sphere-like 초기 상태에 명시적 perturbation이 적용된 경우만 `B. Perturbation-response`로 라벨링한다.
- non-overlapping, unperturbed, unforced run은 `none`으로 라벨링한다. exclusion force가 없는 상태에서 scalar phase 할당만으로는 perturbation이 아니다.

baseline dynamics:

- periodic cell fixed
- external compression/pressure term off
- pure exclusion force only
- overdamped 또는 weakly damped dynamics를 protocol label로 분리

핵심 observable:

- contact survival `S_E(t)`
- graph edit distance from `t0`
- motif/cycle birth-death events
- primitive triangle/square count over time
- `N_T(t)`, `N_O(t)`
- `S_T(t)`, `S_O(t)`
- graph transition events

보조 observable:

- displacement MSD
- effective velocity/relaxation energy

산출물:

- `pressure_off_trajectory_*.h5`
- `pressure_off_summary.csv`
- `dynamic_graph_timeseries.csv`
- `dynamic_motif_timeseries.csv`

## 12. Stage 9: Scalar phase proxy coupling 및 lambda/phase-seed sweep

목적:

- 모든 단위공간이 같은 고유 진동수 `omega_Q`를 가지되 초기 위상이 다른 경우 movement, graph transition, motif survival이 어떻게 달라지는지 본다.
- 이 stage는 full SU(2) quaternion dynamics가 아니라 scalar phase proxy이다.
- 목적은 vertex-internal phase mismatch가 exclusion stiffness, graph transition, O-motif rarity/survival을 bias할 수 있는지 보는 controlled pilot이다.

minimal model:

```text
theta_i(t) = omega_Q * t + theta_i(0)
theta_i(0) ~ Uniform(0, 2*pi)
k_ij(t) = k_core * [1 + lambda * cos(theta_i(t) - theta_j(t))]
0 <= lambda < 1
```

실험 조건:

- selected random-compressed states
- selected jammed/frustrated states
- HCP/FCC states
- `lambda = 0.0, 0.1, 0.3, 0.5`
- multiple phase seeds

metric:

- `P(O > 0 | phi, lambda)`
- `E[N_O/N | phi, lambda]`
- `S_E(t | lambda)`
- `S_T(t | lambda)`
- `S_O(t | lambda)`
- phase coherence of all contacts
- phase coherence of O-motif edges

산출물:

- `phase_model_spec.md`
- `phase_sweep_summary.csv`
- `phase_sweep_by_lambda.csv`
- `phase_conditioned_o_rarity.csv`
- `phase_conditioned_survival.csv`

후속 pilot:

```text
Stage 9b: SU(2) relative-frame coupling pilot
h_ij = q_i^{-1} q_j
k_ij = k_core * [1 + lambda * f(h_ij)]
f(h_ij) = 0.5 * Re Tr(h_ij) 같은 단순 similarity function부터 시작
```

## 13. Stage 10: Pressure-off 전 small-vibration perturbation sweep

목적:

- 압력 제거 전에 격자 또는 compressed state 내부에 작은 진동/위치 perturbation을 넣고 같은 dynamics를 반복한다.

perturbation:

```text
x_i <- x_i + A * u_i
u_i: random unit vector or selected lattice vibration mode
A / ell_Q = 1e-4, 1e-3, 1e-2
```

protocol:

- random micro-jitter
- acoustic-like long-wavelength mode
- HCP/FCC phonon-like mode
- phase off / phase on 모두 가능하되 별도 label로 저장

산출물:

- `perturbed_initial_state_*.h5`
- `perturbation_metadata.csv`
- `perturbed_pressure_off_summary.csv`
- `perturbed_dynamic_motif_timeseries.csv`

## 14. Stage 11: Baseline/phase/perturbed/HCP/FCC dynamics 비교

목적:

- dynamics 조건과 초기 구조 조건의 차이를 한 표로 정리한다.

비교군:

- baseline pressure-off
- phase-conditioned pressure-off
- pre-vibrated pressure-off
- phase-conditioned + pre-vibrated pressure-off
- random-compressed
- HCP
- FCC

metric:

- degree distribution distance
- cycle distribution distance
- T/O motif density difference
- graph edit distance to HCP/FCC contact graph
- dynamics survival difference

산출물:

- `dynamics_comparison_summary.csv`
- `random_vs_lattice_graph_comparison.csv`
- `random_vs_lattice_dynamics_comparison.csv`
- `dynamics_comparison_report.md`

## 15. Stage 12: Packing ratio별 dynamics sensitivity

목적:

- 같은 dynamics protocol이 `phi`에 따라 어떻게 달라지는지 본다.

metric:

- `MSD_final(phi)`
- `S_E_final(phi)`
- `S_T_final(phi)`
- `S_O_final(phi)`
- graph transition count vs `phi`

산출물:

- `phi_dynamics_sensitivity.csv`
- `phi_dynamics_sensitivity_report.md`

## 16. Stage 13: Large-N compression graph atlas scaling

목적:

- compression graph atlas가 system size에 대해 안정적인지 본다.

smoke matrix:

```text
N = 64, 128
seeds = 3
phi_targets = [0.30, 0.50, 0.64]
states = [random_compressed, HCP, FCC]
dynamics = [baseline, phase, perturbed]
```

scaling matrix:

```text
Pilot: N = 256, seeds = 10
Production: N = 512, seeds = 20..50
Optimized O-detector 이후: N = 1024, seeds = 10..30
```

O-motif 검출은 brute-force `N^6`를 금지하고, contact graph의 local neighborhood/2-hop candidate pool에서만 후보를 생성한다.

산출물:

- `scaling_graph_atlas_summary.csv`
- `scaling_graph_transition_summary.csv`

## 17. Stage 14: Large-N pressure-off dynamics scaling

목적:

- pressure-off 움직임과 motif survival이 `N`에 따라 유지되는지 본다.

산출물:

- `scaling_pressure_off_summary.csv`
- `scaling_motif_survival_summary.csv`

## 18. Stage 15: O-motif rarity 및 abundance-like observable

목적:

- O-motif가 드문 구조인지, phase/perturbation/HCP/FCC에서 occurrence가 어떻게 달라지는지 정량화한다.

metric:

- `P(O > 0 | phi, N, protocol)`
- `E[N_O/N | phi, N, protocol]`
- `Var(N_O/N | phi, N, protocol)`
- O cluster size distribution
- O survival after pressure-off dynamics

산출물:

- `o_rarity_by_protocol.csv`
- `o_abundance_observable_report.md`

## 19. Stage 16: T/O motif survival 및 transition

목적:

- pressure-off dynamics 중 T/O motif가 보존, 붕괴, 생성, 전환되는 경로를 기록한다.

metric:

- motif birth/death events
- T-to-O adjacency changes
- motif cluster split/merge
- motif lifetime distribution

산출물:

- `motif_lifetime_distribution.csv`
- `motif_transition_events.csv`

## 20. Stage 17: Graph topology transition 분석

목적:

- 단위공간의 움직임이 graph topology 변화로 어떻게 나타나는지 파악한다.

metric:

- edge birth/death rate
- graph edit distance
- component stability
- cycle birth/death rate
- motif transition correlated with graph transition

산출물:

- `graph_topology_transition_summary.csv`
- `graph_transition_events.csv`

## 21. Stage 18: Final tables/figures/report

최종 표:

- `compression_graph_atlas_summary.csv`
- `cycle_motif_by_phi_summary.csv`
- `selected_state_summary.csv`
- `pressure_off_dynamics_summary.csv`
- `phase_perturbation_comparison_summary.csv`
- `random_vs_lattice_summary.csv`
- `scaling_summary.csv`

최종 figure:

- graph topology vs `phi`
- triangle/square cycle count vs `phi`
- T/O density vs `phi`
- pressure-off MSD/energy/contact survival
- phase sweep comparison
- perturbation comparison
- HCP/FCC vs random-compressed comparison

## 22. Stage 19: 최종 해석과 후속 계획 분리

해석 질문:

- 압축된 단위공간은 각 `phi`에서 어떤 spatial graph family를 만드는가?
- jammed/frustrated 상태의 cycle/T/O motif signature는 무엇인가?
- HCP/FCC와 random-compressed graph는 같은 구조 계열인가, 다른 계열인가?
- 압력을 제거했을 때 배제력만으로 graph와 motif가 얼마나 유지되는가?
- quaternion phase와 pressure-off 전 small vibration은 movement와 motif survival을 유의미하게 바꾸는가?

후속 계획:

- boundary, finite-wall, peeling 실험은 이번 결과와 섞지 않고 별도 plan으로 분리한다.
