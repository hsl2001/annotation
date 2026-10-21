# `anno` 일반화 유전자 예측 모델 고도화 — 작업 인계 문서 (Handoff Guide)
`"HANDOFF.md를 읽고 Step 1부터 이어서 학습, 하이퍼파라미터/다종 일반화 튜닝, 최종 벤치마크 검증까지 완료해 줘."`
## 1. 프로젝트 목표 및 핵심 설계 원칙

* **최종 목표**: 어떤 유전체(*Oryza sativa* `nip`, *Arabidopsis thaliana* `col`, *Solanum lycopersicum* `slm` 등)를 입력하더라도 **80–90% 이상의 Gene-level F1 score**를 기대할 수 있는 종 불변(species-agnostic) 해석 가능 유전자 예측 모델(`anno.model`) 구축.
* **핵심 설계 제약 (중요)**:
  * **"코딩 엑손이 3bp 주기로 진동하고 위상이 유지된다"는 가설은 사용하지 않음** (`plots-tair12/introns/start.png`, `end.png` 실측 결과 엑손 및 인트론 내부의 CWT 파워는 평탄하며, 유의미한 CWT 신호는 오직 **경계(Boundary)의 국소적 임펄스/골(Trough) 파형**에 집중되어 있음).
  * 모든 특징은 다음 3가지 실측 불변량에 기반함:
    1. **4대 경계(`5ss`, `3ss`, `start_codon`, `stop_codon`)의 국소 베이스라인 정규화 CWT 파형** ($\Delta P_a, D^{(1)}_a, D^{(2)}_a, D^{(3)}_a$)
    2. **경계 서열 모티프 × CWT 파형 교차항 (`Motif × CWT` Cross-Terms)**: 진짜 경계(`Real` `GT`/`AG`/`ATG`/Stop)와 가짜 서열(`Control` decoy `GT`/`AG`/`ATG`/Stop)을 구분하는 `k2 × d1_mean`, `k2 × d2_mean`, `PWM_nuc × d2_mean`
    3. **국소 GC 대비(Local GC Contrast) 및 GC 드리프트 적응 특징**: 종 전체의 평균 GC가 달라도($\text{GC}_{\text{nip}}=44.0\%$, $\text{GC}_{\text{col}}=36.3\%$, $\text{GC}_{\text{slm}}=34.5\%$) 모든 종에서 절대적으로 보존되는 $\text{GC}_{\text{CDS}} > \text{GC}_{\text{local\_bg}} > \text{GC}_{\text{intron}}$ 불변량 활용.

---

## 2. 실행 환경 및 필수 명령어 참고사항

1. **외부 평가 도구 경로 (`gffcompare`, `blastp`)**:
   `accuracy.py`를 실행할 때는 반드시 아래 `PATH`를 먼저 설정해야 합니다.
   ```bash
   export PATH="/Users/hyunsu.lim/micromamba/envs/anno/bin:/Users/hyunsu.lim/micromamba/envs/py/bin:$PATH"
   ```
2. **컴파일 및 단위 테스트**:
   macOS에서 `Makefile`의 `ulimit -s unlimited`가 실패하므로 `make test` 대신 아래 명령어로 빌드 및 테스트를 수행합니다.
   ```bash
   make all test-runner && ./test-runner
   ```
3. **스레드 수 권장사항 (`-p 6`)**:
   현재 머신은 8GB RAM macOS 환경입니다. 학습 시 각 워커 스레드가 약 250MB의 DP 버퍼(`Work`)를 사용하므로, 스왑(compressor thrashing) 없이 최고 속도를 내려면 `-p 6` (또는 `-p 4`)을 사용하는 것이 `-p 10`보다 훨씬 빠릅니다.

---

## 3. 완료된 작업 내역 (Completed Work)

### 3.1. 베이스라인 정밀 측정 및 교차 종 오류 원인 규명 (`baseline-analysis` 완료)

* `truth.nip.chr5.3m`에서 300 iteration 학습한 기존 아키텍처 모델(`/tmp/nip3m_300.model`, `/tmp/nip3m_300.weights.txt`)의 성능:
  * **In-species (`truth.nip.chr5.3m`)**: Gene-level F1 **87.54%** (Precision 92.65%, Sensitivity 82.97%), CDS F1 **77.79%**, Intron F1 **83.20%**
  * **Zero-shot (`truth.col.chr5` 전체, 30.1Mbp, 6,253 genes)**: Gene-level F1 **52.27%** (Precision 70.93%, Sensitivity 41.39%), CDS F1 **61.64%**, Intron F1 **67.32%** (참고: 동일 데이터에서 `augustus-rice`의 Gene-level F1은 **29.72%**)
* **규명된 3대 병목 원인**:
  1. **미관측 후보 특징 동결(`freq == 0` freezing) 버그**:
     `feature_counts()`가 정답 경로(`gold_y`)에서 관측된 특징만 카운트했기 때문에, 비표준 스플라이스 부위(`5ss` offset 0의 `GT`/`GC` 외 14개 다이뉴클레오타이드, `3ss` offset -2의 `AG` 외 15개 다이뉴클레오타이드)와 문법상 불가능한 전이(`4: start->stop`, `11: stop->5ss`, `14: 5ss->3ss`, `18: 3ss->start`, `20: 3ss->stop`)가 `freq == 0`으로 비활성화되어 가중치가 `0.0`에 고정되어 있었음.
  2. **종간 GC 이동에 따른 대규모 유전자 병합(Gene Merging) 오류**:
     `col.chr5` 전체 정답 유전자 6,253개 중 **34.1%인 2,135개의 유전자가 906개의 거대 예측 유전자로 병합(merge)**됨(단독 정상 예측은 2,906개). 원인은 애기장대(`col`)와 토마토(`slm`)의 intergenic 영역이 벼(`nip`)보다 훨씬 AT-rich하여, 벼에서 학습된 모델이 이웃한 두 유전자 사이의 intergenic 간격(중앙값 1,163 bp)을 거대한 인트론(bridging intron)으로 잘못 연결하기 때문임.
  3. **Soft-F1 목적함수의 Intergenic-Intron FP 비대칭성**:
     `contig_gradient_f1()`에서 Soft-F1이 `is_cds(st)`만 추적하고 `gold_state == 0`(intergenic) 위치에 인트론(`is_intron(st)`)을 예측하는 것에 아무런 FP 페널티를 주지 않아, L-BFGS가 intergenic 구간을 인트론으로 덮는 오류를 억제하지 못했음.
* **사후 보정(Post-hoc Calibration) 실험을 통한 가설 검증 완료**:
  재학습 없이 `/tmp/nip3m_300.model`에서 동결된 비표준 스플라이스/전이 페널티와 GC 드리프트만 보정한 결과, `truth.col.chr5` 전체(30.1Mbp) zero-shot Gene-level F1이 **52.27% → 70.61%** (`/tmp/col.3m`에서 **58.18% → 75.77%**), Intron F1이 **67.32% → 78.60%**로 즉시 상승함을 확인함.

### 3.2. 3개 종(`nip`, `col`, `slm`) CWT 경계 파형 전수 스캔 (`cwt-boundary-scan` 완료)

* `nip`, `col`, `slm` 각 3Mbp 구간에서 실제 경계(`Real`) vs 비경계 동형 모티프(`Control`)의 CWT 파워(국소 베이스라인 차감 후) 및 1·2차 미분(`d1`, `d2`)의 종간 최소 효과 크기($\min_{\text{species}} |d|$)를 측정함:
  * **`3ss` (Acceptor)**: offset `-6..-1`에서 $\min |d| = 0.34$ (power) 및 $0.30$ (`d2`), branchpoint 영역인 offset `-26..-24`에서 $\min |d| = 0.17$로 3개 종 모두에서 강한 음(-)의 파워 골(trough)과 양(+)의 곡률(`d2`)이 보존됨.
  * **`5ss` (Donor)**: offset `-4..+4`에서 엑손→인트론 급락에 따른 음(-)의 1차 미분(`d1`, $\min |d| = 0.23$) 및 변곡 곡률(`d2`, $\min |d| = 0.23$)이 3개 종 모두에서 보존됨.
  * **`start_codon`**: offset `-4..+4`에서 5'UTR/intergenic→CDS 진입 시점의 강한 곡률(`d2`, $\min |d| = 0.33$) 및 상승 기울기(`d1`, $\min |d| = 0.27$)가 3개 종 모두에서 보존됨.

### 3.3. `crf.h` 및 `test.c` 신규 특징·교차항·최적화 구현 완료 (`implement-features` 완료)

[`crf.h`](./crf.h) 및 [`test.c`](./test.c)에 아래 변경사항이 모두 반영되어 있으며, **29개 단위 테스트(`test-runner`)를 모두 통과**한 상태입니다:

1. **DP 속도 최적화 및 CWT 캐싱 (`contig_gradient_f1`, `backward_step`)**:
   * 내부 루프의 반복적인 `exp()` 호출을 사전 지수화(`exp_prev`, `exp_psi`, `inv_L`)로 대체하고 극단값에서 자동 fallback하도록 구현.
   * Forward 패스의 `SiteStream`에서 계산한 `CwtCached`를 `Work` 구조체(`k->cwt`)에 캐싱하여 Backward 패스(`boundary_observe_cached`)에서 `cwt_at_pos` 재계산을 완전히 제거.
   * `feature_counts()`에서 `site_at()` 대신 `SiteStream` 슬라이딩 윈도우를 사용하여 초기 특징 스케일 계산 속도를 대폭 단축.
2. **국소 CWT 베이스라인 차감 (`cwt_at_pos`)**:
   * 유전체 contig(`c->genome_length >= 128`)에서 $\Delta P_a(g) = P_a(g) - \frac{1}{2}(P_a(g-16) + P_a(g+16))$를 적용하여 종간 평균 GC 차이에 따른 CWT 상수 오프셋을 제거.
3. **4대 경계 CWT 파형 확장 (`OFF_BP`, `OFF_BD`)**:
   * 기존 `3ss`(`kind == 1`)에만 적용되던 CWT 파워 및 미분 특징을 4대 경계 전체(`5ss`, `3ss`, `start_codon`, `stop_codon`)로 확장하고, 오프셋을 branchpoint + 코어 경계인 `{-26, -24, -6, -4, -2, -1, 0, 2, 4}` (9개 오프셋)로 확장.
4. **해석 가능한 `Motif × CWT` 교차항 추가 (`OFF_BCROSS_D1`, `OFF_BCROSS_D2`, `OFF_BPWM_CWT`)**:
   * `OFF_BCROSS_D1`: 경계 다이뉴클레오타이드 `boundary_k2` $\times$ 스케일 평균 1차 미분(`d1_mean`)
   * `OFF_BCROSS_D2`: 경계 다이뉴클레오타이드 `boundary_k2` $\times$ 스케일 평균 2차 미분(`d2_mean`)
   * `OFF_BPWM_CWT`: 경계 CWT 오프셋의 염기(`boundary_cwt_nuc`) $\times$ `d2_mean`
5. **국소 GC 적응 특징 추가 (`OFF_GC_DIFF`, `OFF_GC_LVL`, `OFF_E3_GC`, `OFF_TRANS_GC`)**:
   * 반경 `GC_RADIUS = 256` (513 bp 윈도우)의 국소 GC 비율 $g_{\text{local}}(t)$을 슬라이딩 윈도우로 계산하여:
     * `gc_contrast` $= \text{GC}(k_3) - g_{\text{local}}(t)$ (`OFF_GC_DIFF`: 7개 방출 상태별 국소 GC 대비)
     * `gc_dev` $= g_{\text{local}}(t) - 0.40$ (`OFF_GC_LVL`, `OFF_E3_GC`, `OFF_TRANS_GC`: 국소 GC 배경에 따른 상태·코돈·전이 적응)
6. **후보 특징 활성화 (`feature_counts`) 및 Intergenic-Intron Soft-F1 페널티 (`f1_pred_weight`)**:
   * 유전체 학습 윈도우(`c.length >= 64`)에서 모든 비-`N` 방출·경계·전이 후보 특징을 활성화하여 비표준 스플라이스 및 불가능 전이에 대해 L-BFGS가 음(-)의 억제 가중치를 학습하도록 수정.
   * `f1_pred_weight()`에서 `length >= 64 && gold_state == 0 && is_intron(st)`일 때 `0.5`의 FP 가중치를 부여하여 유전자 병합(bridging intron)을 직접 억제.
   * (`test.c`의 합성 단위 테스트들은 모두 `length <= 18`, `genome_length <= 72`이므로 기존 단위 테스트 동작과 100% 하위 호환됨.)

---

## 4. 남아 있는 작업 및 수행 절차 (Next Steps for Agent)

다음 에이전트는 아래 순서대로 작업을 이어서 완수하면 됩니다.

### Step 1: 신규 아키텍처 단일 종(`nip.chr5.3m`) 300-iter 학습 및 교차 종 성능 확인

먼저 현재 컴파일된 `./anno`로 `truth.nip.chr5.3m`에서 300 iteration 학습을 완료하고 `nip`, `col`, `slm` 성능을 확인합니다.
(참고: `/tmp/col.3m.{fasta,gff3}`와 `/tmp/slm.10m.{fasta,gff3}`가 이미 생성되어 있다면 그대로 사용하고, 재부팅 등으로 `/tmp`가 비워졌다면 `truth.col.chr5`의 앞 3Mbp, `truth.slm.chr5`의 앞 10Mbp를 잘라 생성할 수 있습니다.)

```bash
export PATH="/Users/hyunsu.lim/micromamba/envs/anno/bin:/Users/hyunsu.lim/micromamba/envs/py/bin:$PATH"
make all test-runner && ./test-runner

# 300회 학습 (약 6~8분 소요, -p 6 권장)
./anno -p 6 -e 300 -m /tmp/v1.model truth.nip.chr5.3m.fasta truth.nip.chr5.3m.gff3 > /tmp/v1.nip3m.gff3 2> /tmp/v1.train.log

# In-species (nip) 및 Zero-shot (col, slm) 평가
python3 accuracy.py truth.nip.chr5.3m.fasta truth.nip.chr5.3m.gff3 /tmp/v1.nip3m.gff3
./anno -p 6 -m /tmp/v1.model truth.col.chr5.fasta > /tmp/v1.col.gff3
python3 accuracy.py truth.col.chr5.fasta truth.col.chr5.gff3 /tmp/v1.col.gff3
./anno -p 6 -m /tmp/v1.model truth.slm.chr5.fasta > /tmp/v1.slm.gff3
python3 accuracy.py truth.slm.chr5.fasta truth.slm.chr5.gff3 /tmp/v1.slm.gff3
```

### Step 2: 하이퍼파라미터 튜닝 및 다종(Multi-Species) 일반화 학습 (`optimize-hyperparams`)

1. **가중치 및 오류 분석**:
   * `anno.weights.txt`를 확인하여 `boundary_k2_cwt_d1`, `boundary_k2_cwt_d2`, `boundary_pwm_cwt`, `gc_context` 그룹의 학습된 가중치 분포와 부호를 점검합니다.
   * 과적합되는 특징 그룹이 있다면 `feature_counts()`의 후보 특징 초기 scale 값(현재 `0.05`~`0.1`) 또는 `CRF_L2_C2` / 그룹별 스케일을 조정합니다.
   * 유전자 병합(Merge) vs 분할(Split) 비율을 확인하여 `f1_pred_weight()`의 intergenic-intron 페널티 계수(현재 `0.5`)를 최적값(예: `0.25` ~ `0.75`)으로 조정합니다.
2. **다종(Multi-Species) 학습 데이터 구성 및 일반화 모델 학습**:
   * 단일 종(`nip` 44% GC)만으로도 GC 적응 특징이 작동하지만, `gc_dev`($g_{\text{local}} - 0.40$)에 대한 선형 계수(`OFF_GC_LVL`, `OFF_E3_GC`, `OFF_TRANS_GC`)가 고-GC(외떡잎식물 `nip`)와 저-GC(쌍떡잎식물 `col`, `slm`) 스펙트럼 전체를 정확히 가로지르도록 학습하려면, **`nip` + `col` (+ 필요시 `slm` 일부)의 대표 구간을 결합한 다종 학습 FASTA/GFF3**에서 학습하는 것이 가장 강력합니다.
   * 이때 각 종의 홀드아웃 구간(예: `nip.chr5`의 3Mbp 이후 구간, `col.chr5`의 3Mbp~30Mbp 홀드아웃 구간, `slm.chr5` 전체 또는 홀드아웃 구간)에서 엄격하게 일반화 Gene-level F1을 검증합니다.

### Step 3: 최종 모델 저장, 테스트 검증 및 정리 (`validate-and-test`)

1. 최종 학습된 일반화 모델을 루트의 `anno.model` 및 `anno.weights.txt`로 저장합니다.
2. `make all test-runner && ./test-runner`를 실행하여 모든 단위 테스트가 통과하는지 확인합니다.
3. `truth.nip.chr5.3m`, `truth.col.chr5`, `truth.slm.chr5`에 대한 최종 Gene-level / CDS / Intron F1 벤치마크 표를 정리하고, `test-runner` 바이너리 및 `/tmp`의 임시 파일들을 정리합니다.
