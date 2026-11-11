# Results Summary — Compound Humid Heat in Peninsular Malaysia

> **Status:** pipeline validation on **synthetic data** (10 stations, 1980–2023, 44 years).
> Every number below will be replaced by real-data output; the structure, statistics and
> reporting conventions are final.
>
> Generated from `method_manual/results/*.csv`. Reproduce with notebooks `S01`–`S14`.

---

## Summary of key findings

| # | Finding | Value |
|---|---|---|
| 1 | Event catalogue | 2,485 events across 10 stations × 3 definitions |
| 2 | Antecedent effect on onset | RR 1.52–6.73 (median 3.05); 28/30 cells significant after BH-FDR |
| 3 | Duration–intensity dependence | Kendall τ = 0.295 (Kuantan, Tmax) |
| 4 | Antecedent effect on dependence | pooled Δτ = 0.026, permutation p = 0.44 (**null**) |
| 5 | Copula selection | Gaussian AIC-best; Gaussian/t/Gumbel within ΔAIC < 2 |
| 6 | 50-year joint return level (AND) | 8–9 days, 2.93–3.15 °C above threshold |
| 7 | Definition sensitivity of ranking | RDI 0.08 (Tmax vs Twmax) against a null mean of 0.66 |
| 8 | Non-stationarity check | apparent +79% event increase vanishes under detrending |

---

## Table 1. Event catalogue by heat definition

| Definition | Events | Mean duration (d) | Max duration (d) | Mean intensity (°C) | Max intensity (°C) | Events per station |
|---|---:|---:|---:|---:|---:|---:|
| Tmax | 1,050 | 4.42 | 18 | 1.60 | 4.00 | 90–119 |
| Tmin | 523 | 4.11 | 22 | 2.11 | 5.50 | 44–59 |
| Twmax | 912 | 4.34 | 16 | 1.83 | 4.58 | 80–108 |
| **All** | **2,485** | **4.33** | **22** | **1.79** | **5.50** | — |

*Event definition: calendar-day 90th percentile (1985–2014 baseline, ±7-day window, type-8
quantile); minimum duration 3 days; runs separated by ≤2 days merged. Intensity = peak
exceedance above the day-specific threshold; severity = sum of daily exceedances.*

Per-station detail: `table1_event_summary.csv`.

---

## Table 2. Relative risk of heatwave onset under dry antecedent conditions

RR = P(onset | dry, ARP < 30) / P(onset | wet, ARP > 50), with exposure measured in days.
CI from the log-RR delta method. p values are two-sided day-of-year permutation tests
(1,000 draws); *p*<sub>BH</sub> is Benjamini–Hochberg adjusted across all 30 cells.

| Station | Definition | n | RR | 95% CI | p (perm) | p (BH) | Sig. |
|---|---|---:|---:|---|---:|---:|:--:|
| Alor Star | Tmax | 99 | 2.90 | 1.86–4.53 | 0.001 | 0.001 | * |
| Bayan Lepas | Tmax | 112 | 2.95 | 1.93–4.51 | 0.001 | 0.001 | * |
| Cameron Highlands | Tmax | 112 | 3.15 | 2.04–4.86 | 0.001 | 0.001 | * |
| Chuping | Tmax | 97 | 3.36 | 2.12–5.33 | 0.001 | 0.001 | * |
| Ipoh | Tmax | 119 | 2.77 | 1.81–4.23 | 0.001 | 0.001 | * |
| Kluang | Tmax | 105 | 3.12 | 1.99–4.89 | 0.001 | 0.001 | * |
| Kota Bharu | Tmax | 97 | 1.99 | 1.26–3.16 | 0.030 | 0.032 | * |
| Kuala Krai | Tmax | 108 | 3.06 | 1.94–4.83 | 0.001 | 0.001 | * |
| Kuantan | Tmax | 90 | 4.17 | 2.51–6.91 | 0.001 | 0.001 | * |
| Langkawi | Tmax | 111 | 3.61 | 2.30–5.66 | 0.001 | 0.001 | * |
| Alor Star | Tmin | 46 | 1.57 | 0.84–2.91 | 0.059 | 0.061 | |
| Bayan Lepas | Tmin | 57 | 2.31 | 1.30–4.09 | 0.003 | 0.003 | * |
| Cameron Highlands | Tmin | 59 | 4.01 | 2.10–7.66 | 0.001 | 0.001 | * |
| Chuping | Tmin | 45 | 4.87 | 2.38–9.99 | 0.001 | 0.001 | * |
| Ipoh | Tmin | 59 | 2.24 | 1.27–3.94 | 0.003 | 0.003 | * |
| Kluang | Tmin | 54 | 6.14 | 2.94–12.82 | 0.001 | 0.001 | * |
| Kota Bharu | Tmin | 49 | 1.52 | 0.80–2.89 | 0.245 | 0.245 | |
| Kuala Krai | Tmin | 57 | 4.16 | 2.25–7.71 | 0.001 | 0.001 | * |
| Kuantan | Tmin | 44 | 6.73 | 3.10–14.59 | 0.001 | 0.001 | * |
| Langkawi | Tmin | 53 | 3.03 | 1.59–5.77 | 0.002 | 0.003 | * |
| Alor Star | Twmax | 80 | 2.14 | 1.32–3.45 | 0.003 | 0.003 | * |
| Bayan Lepas | Twmax | 103 | 2.92 | 1.89–4.54 | 0.001 | 0.001 | * |
| Cameron Highlands | Twmax | 90 | 3.12 | 1.95–4.97 | 0.001 | 0.001 | * |
| Chuping | Twmax | 89 | 3.84 | 2.30–6.40 | 0.001 | 0.001 | * |
| Ipoh | Twmax | 108 | 2.10 | 1.36–3.25 | 0.008 | 0.009 | * |
| Kluang | Twmax | 87 | 2.42 | 1.51–3.89 | 0.001 | 0.001 | * |
| Kota Bharu | Twmax | 80 | 2.68 | 1.63–4.39 | 0.002 | 0.003 | * |
| Kuala Krai | Twmax | 95 | 2.48 | 1.57–3.92 | 0.001 | 0.001 | * |
| Kuantan | Twmax | 82 | 3.69 | 2.20–6.17 | 0.001 | 0.001 | * |
| Langkawi | Twmax | 98 | 3.69 | 2.27–6.00 | 0.001 | 0.001 | * |

**Median RR by definition:** Tmax 3.09 (1.99–4.17) · Tmin 3.52 (1.52–6.73) · Twmax 2.80 (2.10–3.84).

Event-attribute comparisons (duration and intensity, dry vs wet) are in
`table2b_ro2_attributes.csv`; **none** survives BH correction (0/30 for both attributes),
consistent with an occurrence-only antecedent effect.

---

## Table 3. Copula selection and tail dependence

Kuantan, Tmax, n = 90 events. Two-step IFM with mid-interval PIT (τ = 0.2951).
Parameters by maximum pseudo-likelihood.

| Family | k | log L | AIC | ΔAIC | τ (fitted) | λ<sub>L</sub> | λ<sub>U</sub> |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Gaussian** | 1 | 8.64 | **−15.29** | 0.00 | 0.302 | 0.000 | 0.000 |
| t (df = 20.07) | 2 | 8.74 | −13.49 | 1.80 | 0.300 | 0.011 | 0.011 |
| Gumbel | 1 | 8.11 | −14.21 | 1.08 | 0.256 | 0.000 | 0.325 |
| Frank | 1 | 7.42 | −12.84 | 2.45 | 0.282 | 0.000 | 0.000 |
| Joe | 1 | 7.07 | −12.14 | 3.15 | 0.200 | 0.000 | 0.383 |
| Clayton | 1 | 4.86 | −7.72 | 7.57 | 0.295 | 0.437 | 0.000 |

**Model uncertainty:** the ΔAIC < 2 set {Gaussian, t, Gumbel} spans λ<sub>U</sub> = 0 to 0.325.
All results downstream use the ensemble, not the single best model.

*Parametric-bootstrap GoF (CvM, 150 reps):* Gumbel p = 0.74 · t p = 0.47 · Gaussian p = 0.43 ·
Joe p = 0.32 · Frank p = 0.25 · Clayton p = 0.043 (the only rejection). With n = 90 the test has
little power, so a non-rejection should be read as "not contradicted", not as "correct".

---

## Table 4. Does antecedent state change duration–intensity dependence?

Per-cell estimates require ≥15 events in each state; 24 of 30 cells qualify.

| Quantity | Value |
|---|---|
| Pooled Δτ (sample-size weighted) | **0.0259** |
| Permutation p (two-sided, 1,000 draws) | **0.442** |
| Null SD of pooled Δτ | 0.0346 |
| Median per-cell SE(Δτ) | 0.182 |
| Detectable \|Δτ\| at 80% power (2 × median SE) | ≈ 0.36 |
| Cells with Δτ > 0 | 13 / 24 |
| Cells whose 95% CI excludes 0 | 0 / 24 |

**Interpretation:** no evidence that antecedent moisture alters event structure, with power
only for large effects. This contrasts sharply with RO2 (RR 3–7 in the same stratification),
supporting an onset-only mechanism.

---

## Table 5. Joint return periods — Kuantan, Tmax

Event rate λ = 90 events / 44 years = **2.045 events yr⁻¹**; T(years) = T(events)/λ.
Quantiles solved on u = v (equal marginal rarity); margins inverted through the
truncated-geometric (duration) and gamma (intensity) distributions.

| T (yr) | AND: D (d) / I (°C) | OR: D / I | Kendall: D / I | q (AND) | q (OR) | q (Kendall) |
|---:|---|---|---|---:|---:|---:|
| 5 | 5 / 2.08 | 8 / 2.84 | 6 / 2.24 | 0.775 | 0.945 | 0.827 |
| 10 | 6 / 2.36 | 9 / 3.18 | 7 / 2.51 | 0.860 | 0.973 | 0.894 |
| 25 | 7 / 2.69 | 11 / 3.61 | 8 / 2.85 | 0.925 | 0.990 | 0.947 |
| 50 | 8 / 2.93 | 13 / 3.91 | 9 / 3.11 | 0.954 | **0.995** | 0.968 |

**Ensemble (ΔAIC < 2), AND definition:**

| T (yr) | Gaussian | t | Gumbel | Range |
|---:|---|---|---|---|
| 10 | 6 d / 2.36 °C | 6 d / 2.36 °C | 6 d / 2.43 °C | 6 d / 2.36–2.43 °C |
| 50 | 8 d / 2.93 °C | 8 d / 2.96 °C | 9 d / 3.15 °C | 8–9 d / 2.93–3.15 °C |

> ⚠️ **Extrapolation flag.** The empirical PIT ceiling is 1 − 1/(n+1) = 0.989, so
> q<sub>OR</sub> = 0.995 at T = 50 yr lies beyond the data. Return levels at T ≥ 25 years are
> model-based extrapolations and must be reported as such.

---

## Table 6. Definition sensitivity of the station risk ranking

Primary metric: annual mean severity (°C·days yr⁻¹). RDI = Σ|r<sub>A</sub> − r<sub>B</sub>| / (n²/2);
0 = identical, 1 = reversed. **Null mean for n = 10 independent rankings: 0.66.**

| Metric | Tmax vs Tmin | Tmax vs Twmax | Tmin vs Twmax |
|---|---|---|---|
| Annual severity | 0.32 (ρ = 0.79) | **0.08 (ρ = 0.98)** | 0.32 (ρ = 0.79) |
| Event rate | 0.20 (ρ = 0.92) | 0.28 (ρ = 0.83) | 0.40 (ρ = 0.70) |
| Mean intensity | 0.76 (ρ = −0.12) | 0.72 (ρ = 0.09) | 0.64 (ρ = 0.02) |
| Mean duration | 0.56 (ρ = 0.37) | 0.64 (ρ = 0.15) | 0.48 (ρ = 0.39) |

**Year-block bootstrap (2,000 resamples) for the severity ranking:**

| Pair | Point | Bootstrap median | 95% interval |
|---|---:|---:|---:|
| Tmax vs Tmin | 0.32 | 0.32 | 0.12–0.56 |
| Tmax vs Twmax | 0.08 | 0.20 | 0.08–0.44 |
| Tmin vs Twmax | 0.32 | 0.36 | 0.12–0.64 |

**Top-3 stability:** Ipoh and Bayan Lepas appear in the top three under all three definitions
(2/3 shared). Kuantan ranks 10th under all three.

---

## Table 7. Sensitivity to methodological choices

| Specification | n events | Median duration (d) | Mean intensity (°C) | τ | RR (Kuantan Tmax) | RDI (severity) |
|---|---:|---:|---:|---:|---:|---:|
| **Baseline** (p90, ±7 d, min 3 d, merge ≤2 d, 1985–2014) | 2,485 | 4 | 1.79 | 0.295 | 4.17 | 0.08 |
| Half-window ±14 d | 2,510 | 4 | 1.81 | 0.344 | 4.81 | 0.16 |
| 95th percentile | 874 | 3 | 1.55 | 0.239 | 7.29 | 0.28 |
| Min 5 d, no merging | 705 | 5 | 2.13 | 0.243 | 3.61 | 0.36 |
| Base period 1991–2020 | 2,110 | 4 | 1.77 | 0.259 | 4.21 | 0.24 |
| Detrended threshold | 2,411 | 4 | 1.75 | 0.218 | 4.96 | 0.20 |

Relative spread (max − min)/mean: RDI 1.27 · n events 0.98 · RR 0.76 · median duration 0.50 ·
τ 0.47 · mean intensity 0.32.

**Reading:** signs and orders of magnitude are stable (RR > 1, τ > 0, RDI ≪ 0.66 in every
variant); magnitudes are not. Report directions and ranges, not point estimates.

### Non-stationarity

| Period | Events per decade, fixed baseline | Events per decade, detrended threshold |
|---|---:|---:|
| 1980–2001 | 0.72 | 0.99 |
| 2002–2023 | 1.28 | 1.01 |

The apparent +79% increase is entirely attributable to baseline drift (regional trend
+0.018 °C yr⁻¹, +0.77 °C over the record).

---

## Figures

All figures are in `method_manual/figs/` (PNG, 200 dpi).

| File | Content |
|---|---|
| `S01_acf.png` | Autocorrelation of the daily series |
| `S01_missing_map.png` | Missing-data map by station and year |
| `S01_seasonal_cycle.png` | Seasonal cycle of Tmax / Tmin / Twmax |
| `S01_trend.png` | Annual series and trend |
| `S02_threshold_curve.png` | Calendar-day threshold curve |
| `S02_seasonal_bias.png` | Why a calendar-day window is required |
| `S03_event_timeline.png` | Runs, declustering and event attributes |
| `S04_definition_attributes.png` | Event attributes under three definitions |
| `S05_event_counts.png` | Event counts by station and definition |
| `S05_event_calendar.png` | Event timeline (Gantt) |
| `S05_events_by_decade.png` | Decadal event counts |
| `S06_arp_distribution.png` | ARP distribution, events vs base days |
| `S06_arp_vs_cdd.png` | Joint ARP–CDD structure |
| `S06_state_composition.png` | Dry/neutral/wet composition of events |
| `S07_forest_rr.png` | Relative risk forest plot (30 cells) |
| `S07_permutation_null.png` | Day-of-year permutation null |
| `S08_raw_scatter.png` | Duration–intensity on the original scale |
| `S08_copula_scale.png` | Same events on the copula scale |
| `S08_family_shapes.png` | Six copula families at equal τ |
| `S08_observed_tail.png` | Observed upper-corner enrichment |
| `S09_pit_duration.png` | Average-rank vs randomised PIT |
| `S09_pit_intensity.png` | PIT of the intensity margin |
| `S10_fit_check.png` | Observed vs simulated copula |
| `S11_delta_tau_forest.png` | Per-cell Δτ with bootstrap CI |
| `S11_perm_delta.png` | Permutation null for pooled Δτ |
| `S12_return_period_contours.png` | AND joint return-period contours |
| `S13_rank_bump.png` | Station rank across three definitions |
| `S13_rdi_null.png` | RDI null distribution for n = 10 |
| `S14_ns_baseline.png` | Fixed vs detrended baseline |
| `S14_enso_proxy.png` | Warm-year proxy analysis |
| `S14_sensitivity.png` | Six-specification sensitivity panel |

## Reporting notes for the manuscript

1. **Effect sizes before p values.** Every inferential statement pairs an estimate with a CI.
2. **The RO3 null is a result, not a failure.** Report Δτ with its detectable-effect threshold.
3. **λ<sub>U</sub> is a model property.** Never report it without naming the copula ensemble.
4. **T ≥ 25 yr is extrapolation.** Flag it in the table caption, not only the text.
5. **Ranking conclusions are conditional.** RDI moves from 0.08 to 0.36 across defensible
   specifications (Table 7); phrase the conclusion accordingly.
