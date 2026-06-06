# Qaether Packing Simulator - Developer & Agent Guide (`agents.md`)

이 문서는 Qaether 단위공간 압축/그래프/motif/dynamics 시뮬레이터 프로젝트(`Packing_Simulator`)에서 작업을 수행하는 AI 에이전트 및 개발자를 위한 코드 작성, 디버깅, 분석 및 기여의 기본 지침서입니다. 향후 진행할 모든 기능 개발 및 리팩토링은 이 문서의 원칙과 정의를 준수해야 합니다.

---

## 1. 프로젝트 핵심 개념 및 목적

본 시뮬레이터는 Qaether 단위공간의 **패킹 비율(Packing Ratio, $\phi$)**에 따른 **공간 그래프(Spatial Contact Graph)** 형성, **기하학적 모티프(Cycle, T-motif, O-motif)**의 분포 및 **압력 해제 시의 동적 생존력(Pressure-off Survival)**을 탐구합니다.

> **핵심 해석 원칙**: Qaether 시뮬레이션은 "공간 안의 입자 운동"이 아니라, **최소 단위공간의 임베딩 그래프/모티프 응답 실험**입니다. 근본 관측 대상은 좌표($\rho_i$)가 아니라 $G=(V,E)$, primitive cycle, T/O motif 및 그 survival/transition입니다.

### 핵심 물리량 및 Proxy 해석 기준

*   **Qaether 단위공간**: Vertex 또는 단위공간의 중심을 의미합니다.
*   **Effective Exclusion Proxy**: Qaether 단위공간의 배제 체적(Exclusion volume)을 표현하기 위해 구(Sphere) 모델을 차용합니다.
    *   $\ell_Q = 1.0$ (Qaether 기본 단위 길이)
    *   $R_Q = \ell_Q / 2 = 0.5$ (단위공간 유효 반경)
*   **`positions` (좌표)**: 실제 물리적 공간의 입자 위치가 아니라, **접촉 그래프(Contact Graph)를 생성하기 위해 계산적으로 임베딩한 Proxy 좌표**입니다.
*   **`velocities` (속도)**: 실제 미시적 입자의 물리적 속도가 아니라, **수치적 이완(Relaxation) 프로토콜에서 좌표를 업데이트하기 위해 도입한 유효 변수(Effective update variable)**입니다.
*   **패킹 분율(Packing Fraction, $\phi$)**:
    $$\phi = \frac{N \times \frac{4}{3}\pi R_Q^3}{V_{\text{cell}}}$$
    (주기적 경계 조건 periodic unit cell 내의 유효 패킹 분율)

### $\phi_{\text{target}}$ vs $\phi_{\text{achieved}}$ 분리 규칙

압축 실험에서는 반드시 다음 두 값을 분리하여 기록해야 합니다:
*   **$\phi_{\text{target}}$**: 외부 압축 프로토콜이 시도한 패킹 비율.
*   **$\phi_{\text{achieved}}$**: 실제 스냅샷에서 기록된 유효 패킹 비율.

특히 $\phi > 0.64$ 구간에서는 harmonic soft-core 모델 때문에 overlap을 허용한 **soft-sphere compression artifact**가 발생할 수 있으므로, 다음 판정을 적용합니다:

```text
max_overlap <= tolerance AND overlap_energy_per_vertex <= tolerance
  → hard_sphere_like
그 외
  → soft_overcompressed_or_stressed (achieved hard packing이 아님)
```

추가 overlap 진단 지표: `max_overlap`, `mean_overlap`, `overlap_energy_per_vertex`, `fraction_overlapping_pairs`.

### Jamming vs Frustration 정의

이 두 개념은 명확히 분리되어야 합니다:

*   **Jamming** = **기계적 구속 안정성(Mechanical constraint stability)**
    *   판정 지표: overlap energy plateau, displacement(MSD) 감소, coordination number($Z$), force-balance residual, pressure estimator($P$), bulk modulus proxy($dP/d\phi$)
    *   추가 참고: rattler가 있을 수 있으므로 `Z_rattler_removed`를 따로 두는 것이 권장됨
*   **Frustration** = **비호환 로컬 구속에 의한 압축/재배열 실패**
    *   판정 지표: graph transition 실패, motif incompatibility, `failed_rearrangement_count`, `graph_edit_rate` plateau under compression attempts

### Pressure-off Dynamics의 두 가지 해석 모드

현재 force model은 순수 배제력이므로 모든 쌍이 $r_{ij} \ge \ell_Q$이면 힘이 0이 되어 dynamics가 trivial합니다. 따라서 pressure-off 실험은 반드시 다음 두 모드 중 하나로 라벨링해야 합니다:

*   **A. Residual-stress relaxation**: 압축 상태에 미세 overlap 또는 잔류 응력이 남아 있고, pressure-off 후 이것이 이완되는 경우. ($r_{ij} < \ell_Q \Rightarrow F_{ij} \neq 0$)
*   **B. Perturbation-response**: 겹침이 거의 없어 pure exclusion만으로는 움직임이 없으므로, 작은 변위/속도 교란을 인가한 후 graph/motif survival을 관측하는 경우.

> 핵심 관측량의 우선순위는 다음과 같습니다:
> *   **Primary**: $S_E(t)$, graph edit distance, $S_T(t)$, $S_O(t)$, cycle/motif birth-death
> *   **Secondary**: MSD, effective velocity, relaxation energy
>
> 이론 목적상 MSD는 보조 지표입니다. 핵심은 "얼마나 움직였나"가 아니라, **공간 그래프와 motif가 유지되었는가**입니다.

---

## 2. 실험 Phase 구분

전체 Stage는 다음 네 Phase로 그룹핑됩니다. 각 Phase를 순서대로 완료해야 실패 원인 추적이 용이합니다.

| Phase | Stage | 목표 |
|:---|:---|:---|
| **Phase A**: Geometry-only atlas | Stage 0–7 | $\phi \to G_\phi \to C_3, C_4, T, O$ 확정. Quaternion, perturbation, dynamics를 모두 끔 |
| **Phase B**: Pressure-off survival | Stage 8, 10, 11, 12 | $G(t)$, $T(t)$, $O(t)$ survival과 transition 관측 |
| **Phase C**: Phase coupling pilot | Stage 9 | $\lambda$가 O-motif rarity/survival에 영향을 주는지 확인 |
| **Phase D**: Scaling & transition | Stage 13–17 | $N$ 증가 시 결과 안정성 확인 |

---

## 3. 디렉토리 구조 및 주요 모듈 설명

```text
Packing_Simulator/
├── qaether_sim/                # 핵심 시뮬레이션 라이브러리
│   ├── __init__.py             # 패키지 진입점 및 공개 인터페이스 정의
│   ├── config.py               # 실험 파라미터 및 오차 허용치(Tolerances) 정의
│   ├── state.py                # QaetherState 정의 (HDF5 입출력 포함)
│   ├── geometry.py             # PBC(주기적 경계 조건) 기반 거리 및 기하학적 unwrapping
│   ├── forces.py               # 배제력(Harmonic Overlap Force) 및 위상 의존 스티프니스 계산
│   ├── initial_conditions.py   # Random Gas, FCC, HCP 격자 생성기
│   ├── compression.py          # 압축 스윕(Rescaling) 및 이완(Relaxation) 루프
│   ├── contact_graph.py        # NetworkX 기반 접촉 그래프 생성 (일반 및 히스테리시스 규칙)
│   ├── cycles.py               # 원시 사이클(Primitive Triangle/Square) 검출
│   ├── motifs_T.py             # Tetrahedron(T) 모티프 검출
│   ├── motifs_O.py             # Octahedron(O) 모티프 검출
│   ├── jamming.py              # 압축 흔적(MSD, Energy, Pressure) 기반 잼/프러스트레이션 판정
│   ├── lattice_benchmarks.py   # FCC/HCP 레퍼런스 및 최밀 충전 제한 정의
│   ├── bulk_dynamics.py        # 압력 해제(Pressure-off) 상태의 시간 발전 및 생존율 측정
│   ├── phase.py                # 스칼라 위상 변수 할당 프로토콜
│   ├── perturbation.py         # 압력 해제 전 미세 지터(Jitter) 인가
│   └── analysis.py             # 전체 파이프라인 수립 및 실행 제어
├── scripts/                    # 실행 스크립트 및 테스트 유틸리티
│   ├── run_smoke_experiment.py # 스모크 테스트 실행용 스크립트
│   └── run_tests.sh            # 유닛 테스트 및 스모크 테스트 자동화 스크립트
└── tests/                      # 테스트 스위트
    └── test_smoke_pipeline.py  # 모티프 탐지, HDF5 저장, 위상 바이어스 등 유닛 테스트
```

---

## 4. 핵심 모듈별 상세 명세 및 규칙

### 4.1. `qaether_sim/config.py` (`ExperimentConfig`)
*   시뮬레이션 전반의 물리 상수 및 접촉/모티프 기하 판정 입실론(오차 범위)을 보관합니다.
*   **주요 변수**:
    *   `epsilon_contact`: 접촉 그래프 엣지를 형성하는 거리 임계치 ($r_{ij} \le \ell_Q \times (1 + \text{epsilon\_contact})$)
    *   `epsilon_cycle`: primitive cycle 판정을 위한 길이 오차 허용 범위
    *   `epsilon_planar`: 4개 점이 동일 평면에 있는지를 확인하기 위한 SVD 기반 평면 판정 오차 범위
    *   `epsilon_volume`: T-motif의 부피 하한선 (0에 수렴하는 평평한 사면체 방지)
    *   `epsilon_center` & `epsilon_perp`: O-motif 대각축 교차 및 직교성 판정 오차 범위

### 4.2. `qaether_sim/state.py` (`QaetherState`)
*   모든 물리 정보(`positions`, `velocities`, `box`, `theta`)와 메타데이터(`metadata` dict)를 관리합니다.
*   **HDF5 직렬화**: `.save_h5()` 및 `.load_h5()`를 통해 시뮬레이션 상태를 완벽하게 보존 및 로드합니다. 메타데이터는 JSON 문자열로 저장되어 데이터 정합성을 유지합니다.

### 4.3. `qaether_sim/geometry.py`
*   **Minimum Image Convention**: 주기적 경계 조건(PBC) 하에서 가장 가까운 이미지와의 변위 및 거리를 계산합니다.
*   **Coordinate Unwrapping (`centered_points`)**: **가장 중요한 기하 계산 보조 함수**입니다. PBC 경계선에 걸쳐진 여러 입자들을 기준 입자(Anchor)를 중심으로 unwrapping하여 실제 인접 로컬 배치를 복원합니다. **사이클 및 모티프 판정 전 필수적으로 수행**해야 합니다.

### 4.4. `qaether_sim/forces.py`
*   **배제력 포텐셜 (Overlap Force)**:
    $$U_{\text{core}}(r) = \frac{1}{2} k ( \ell_Q - r )^2 \quad (\text{단, } r < \ell_Q)$$
    이 식의 미분값을 이용하여 조화 척력을 인가합니다.
*   **위상 의존 스티프니스 (Scalar Phase Proxy Coupling)**:
    $$k_{ij} = k_{\text{core}} [1 + \lambda \cos(\theta_i - \theta_j)]$$
    스칼라 위상 $\theta$의 조화 결합 정도 $\lambda$에 따라 강성(Stiffness)이 변조됩니다.
    > ⚠️ **이것은 full SU(2) quaternion dynamics가 아닙니다.** vertex-internal phase mismatch가 exclusion stiffness, graph transition, O-motif survival을 bias할 수 있는지를 테스트하기 위한 **controlled scalar proxy**입니다. full SU(2) coupling ($h_{ij} = q_i^{-1} q_j$, $k_{ij} = k_{\text{core}} [1 + \lambda f(h_{ij})]$)은 Stage 9b 후속 pilot으로 분리되어 있습니다.

### 4.5. `qaether_sim/contact_graph.py`
*   **Hysteretic Contact Graph**: 동적 실험 시 급격한 엣지 점멸(Flickering)을 방지하기 위해 생성 임계치(`epsilon_on`)와 제거 임계치(`epsilon_off`, `off > on`)를 다르게 두는 히스테리시스 알고리즘을 지원합니다.

### 4.6. `qaether_sim/cycles.py` & `motifs_T.py` & `motifs_O.py`
*   **Primitive Cycle**: 다른 짧은 지름길이 없는 최소 면 단위(Face-like) 구조를 의미합니다.
*   **T-motif (Tetrahedron)**:
    *   4개 정점의 유도 부분그래프(Induced Subgraph)가 $K_4$ (완전 그래프)를 이룸.
    *   6개 변의 길이가 $\ell_Q$ 부근에 위치함.
    *   사면체 부피가 임계치 이상임 (`tetra_volume > epsilon_volume`).
*   **O-motif (Octahedron)**:
    *   6개 정점의 유도 부분그래프 엣지 수가 12개임.
    *   직교하는 3개의 대각 정점 쌍(Disjoint opposite pairs)이 존재하며, 쌍의 중점들이 한 점에 밀집함 (`epsilon_center`).
    *   대각 정점 쌍이 형성하는 3개 축이 거의 직교함 (`epsilon_perp`).
    *   8개의 원시 삼각형이 존재하며, **각 O-edge가 정확히 2개의 O-triangle에 incident**함.
    *   3개의 원시 사각형이 존재하며, **3개의 square cycle이 12개 O-edge를 정확히 한 번씩(exactly once) cover**함.
*   ⚠️ **성능 주의**: 대규모 시뮬레이션($N \ge 256$)에서 O-motif를 brute-force ($O(N^6)$)로 찾는 것은 리소스 낭비가 매우 큽니다. 반드시 접촉 그래프 상의 **로컬 인접 그룹(2-hop neighborhood)을 사전 필터링**한 후 판정해야 합니다.

### 4.7. `qaether_sim/bulk_dynamics.py`
*   **Pressure-off Survival**: 압축 프로토콜이나 등방 압력을 해제하고 유효 배제력 기반 동역학만을 적분할 때, 접촉 그래프의 엣지 생존율 $S_E(t)$, T-motif 생존율 $S_T(t)$, O-motif 생존율 $S_O(t)$의 시계열을 수집합니다.

### 4.8. `qaether_sim/lattice_benchmarks.py`
*   HCP/FCC 격자 benchmark는 다음 두 종류를 분리하여 관리해야 합니다:
    *   **Ideal close-packing benchmark**: 이론적 최밀 충전($\phi \approx 0.7405$) 근처의 참조용 격자.
    *   **TOTO-labeled constructed benchmark**: Qaether의 T/O incidence condition을 만족하도록 구성된 격자. 일반 FCC/HCP contact graph 전체가 자동으로 TOTO 조건을 만족한다고 간주하면 위험합니다.
*   동일 $N$, 동일 periodic cell convention으로 생성하되, `phi_target`과 `phi_achieved/contact scale`을 별도 기록합니다.

---

## 5. 코딩 지침 및 개발 수칙

개발을 진행할 때 다음 원칙을 반드시 지켜야 합니다.

### 1) 기하학적 정합성 및 PBC 언래핑
*   네트워크상의 노드 좌표를 기하학적으로 처리할 때는 항상 `centered_points`를 통해 local patch로 unwrap한 다음 부피, 평면 여부, 표면적 등을 계산하십시오. 원본 `positions` 값을 PBC 보정 없이 그대로 연산에 투입하면 오탐률이 극도로 증가합니다.

### 2) 성능 최적화
*   거리 정렬 및 접촉 엣지 검출 단계에서 데이터가 커질 경우 `scipy.spatial.cKDTree` 또는 `numba` 가속을 적극적으로 활용하여 2중 루프의 오버헤드를 제어해야 합니다.
*   특히 모티프 탐지 모듈(`motifs_O.py`, `motifs_T.py`)의 경우 결합도가 높은 후보 그룹(Clique) 또는 degree 필터링을 앞단에 배치하여 시간 복잡도를 줄여야 합니다.

### 3) 멱등성 및 원본 코드 스타일 보존
*   프로젝트 코드 내 기존의 문서화 주석(docstring)과 변수 해석 방향을 임의로 훼손하지 마십시오.
*   모든 수정 사항은 기존 유닛 테스트(`tests/test_smoke_pipeline.py`) 및 전체 스모크 테스트(`scripts/run_tests.sh`)를 통과해야 합니다.

### 4) 데이터 스키마 및 입출력 형식 고정
*   시뮬레이션에서 생성하는 `.h5` 파일 및 `.csv` 결과 보고서의 컬럼 스펙은 기존 정의(`config.py` 및 `analysis.py` 참고)를 정확히 따라야 합니다.
*   다양한 조건(Random-compressed, HCP, FCC, Perturbed, Phase-conditioned)에서 파생된 데이터 파일은 명확하게 라벨링하여 저장해야 혼선이 없습니다.

### 5) Large-N 스케일링 전략
*   대규모 실험은 단계적으로 진행합니다:
    *   **Smoke**: $N=64, 128$, seeds = 3
    *   **Pilot**: $N=256$, seeds = 10
    *   **Production**: $N=512$, seeds = 20–50
    *   **최적화 후**: $N=1024$, seeds = 10–30
*   O-motif 검출의 brute-force $N^6$은 절대 금지. contact graph의 2-hop neighborhood에서만 후보를 생성해야 합니다.

---

## 6. 실행 및 검증 방법

새로운 코드를 배포하거나 모듈을 업데이트한 경우 반드시 아래의 절차를 통해 확인합니다.

```bash
# 1. 유닛 테스트 및 스모크 파이프라인 전체 실행
bash scripts/run_tests.sh
```

*   `run_tests.sh`는 `tests/` 내부의 모든 유닛 테스트를 구동하고, `scripts/run_smoke_experiment.py`를 활용해 결과 산출 디렉토리(`results_smoke/`)에 다음과 같은 핵심 분석 레포트들이 잘 나오는지 검증합니다.
    *   `energy_curve.csv` (압축 경로)
    *   `compression_graph_motif_atlas.csv` (그래프 및 모티프 아틀라스)
    *   `state_selection_report.csv` (잼/프러스트레이션 판정)
    *   `pressure_off_summary.csv` (압력 해제 동역학)
    *   `phase_sweep_summary.csv` (위상 모델 바이어스 결과)
    *   `perturbed_pressure_off_summary.csv` (지터 교란 결과)
    *   `hcp_fcc_motif_summary.csv` (결정 대조군 분석)
    *   `smoke_summary.json` (요약 JSON 파일)
