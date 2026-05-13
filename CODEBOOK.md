# Codebook

JACSIS variable mapping for the analysis script (`analysis.R`).

This codebook documents how the JACSIS 2022 / 2023 / 2024 / 2025 survey items are mapped onto the analytic variables used in the manuscript. It is intended for researchers granted access to the JACSIS data who wish to reproduce or extend the analysis.

Three panel CSVs are required by the script (each is a wide-format merge on `Monitor_ID`, with each item suffixed by the wave year):

| File                                  | Waves merged       | Used by                                       |
|---|---|---|
| `data/jacsis_2022_2025.csv`           | 2022 + 2025        | Sens 3 (ANCOVA), Sup Table 3                  |
| `data/jacsis_2023_2024_2025.csv`      | 2023 + 2024 + 2025 | Main (DiD), Sens 2 (ANCOVA), Table 1, Sup Tables 1A/1B/4 |
| `data/jacsis_2022_2023_2025.csv`      | 2022 + 2023 + 2025 | Sens 1 (DiD), Sup Table 2                     |

Cross-wave Q-code shifts are common in JACSIS: the **same construct uses different Q-codes in different waves**. Watch in particular for **employment** (`Q5.1` in 2022/2024 but `Q6.1` in 2023) and **household income** (`Q87.1` in 2022, `Q90.1` in 2023, `Q80.1` in 2024, `Q77.1` in 2025). The script's `derive_confounders()` helper uses a per-wave switch.

---

## 1. Exposure (2025 wave only)

The generative-AI use module was first fielded in the 2025 wave. There are no 2022 / 2023 / 2024 Q-codes for the items in this section.

### `Q37S1_2025` — when did you first start using generative AI? (6 categories, single-answer)

| Code | Meaning | Used by the analysis as |
|---|---|---|
| 1 | Never used | Control in Main / Sens 1 / Sens 2; D = 0 in Sens 3 |
| 2 | Used in the past, but quit | **Excluded** from all tiers (heterogeneous, small) |
| 3 | Started Nov 2022 – Dec 2023 | Excluded from Main / Sens 1 / Sens 2; D = 1 in Sens 3 |
| 4 | Started Jan – Dec 2024       | Excluded from Main / Sens 2; **treatment** in Sens 1; D = 1 in Sens 3 |
| 5 | Started Jan 2025+           | **Treatment** in Main / Sens 1 / Sens 2; D = 1 in Sens 3 |
| 6 | Started Jul 2025+           | **Treatment** in Main / Sens 1 / Sens 2; D = 1 in Sens 3 |

### `Q37S3.1`–`Q37S3.9_2025` — purpose-specific frequency of AI use (5-pt Likert, asked of users)

| Sub-code | Purpose label (ja) | English | Composite |
|---|---|---|---|
| Q37S3.1 | 文章作成・SNS投稿等 | writing / SNS composition | creative |
| Q37S3.2 | 翻訳               | translation              | creative |
| Q37S3.3 | 情報の検索・参照   | lookup / reference       | daily    |
| Q37S3.4 | 画像・動画生成     | image / video generation | creative |
| Q37S3.5 | 学習・勉強         | learning / study         | creative |
| Q37S3.6 | 日常生活の計画     | daily-life planning      | daily    |
| Q37S3.7 | 健康情報           | health information       | daily    |
| Q37S3.8 | 会話・雑談         | conversation             | social   |
| Q37S3.9 | 感情面のサポート   | emotional support        | social   |

Response scale (all 9 sub-items, identical): 1 = まったく利用していない (never); 2 = 月に1回未満 (less than once a month); 3 = 月に数回程度 (a few times a month); 4 = 週に1～数回程度; 5 = ほぼ毎日 (almost every day).

(Q37S3.10「その他」is excluded as a residual catch-all.)

### Derived exposures

| Variable | Definition | Scale |
|---|---|---|
| `use`       | binary `D` indicator = 1 if `Q37S1_2025 ∈ {3, 4, 5, 6}` else 0 (tier-specific cohort filter applied upstream — see Q37S1 table above) | 0 / 1 |
| `intensity` | mean of `Q37S3.1`–`Q37S3.9_2025`. Non-users (Q37S1 ∈ {1, 2}) are anchored at 1 ("never"). | 1–5 |
| `creative`  | item mean of `Q37S3.{1, 2, 4, 5}_2025`. Non-users anchored at 1. | 1–5 |
| `daily`     | item mean of `Q37S3.{3, 6, 7}_2025`. Non-users anchored at 1.    | 1–5 |
| `social`    | item mean of `Q37S3.{8, 9}_2025`. Non-users anchored at 1.       | 1–5 |

`intensity` and the three purpose composites are never combined in one model because `intensity = (4·creative + 3·daily + 2·social) / 9` is a perfect linear combination of the three.

The three-factor purpose structure is validated by exploratory factor analysis (maximum-likelihood extraction with oblimin rotation, parallel-analysis n-factor = 3). See `output/manuscript/diagnostics/efa_{loadings,fit}.csv`.

---

## 2. Outcomes (13 continuous outcomes across 7 wellbeing domains)

The 13 outcomes are drawn from the VanderWeele Flourishing Index (FI) and standalone instruments. Each outcome is measured at every wave used by the relevant tier; Q-code prefixes shift across waves.

| # | id (variable in script) | Domain | Outcome | 2022 | 2023 | 2024 | 2025 | Scale |
|--:|---|---|---|---|---|---|---|---|
| 1 | `d1_pm_avg`     | 1. Mental & Physical Health | Physical & mental health (FI #3 + #4 mean)   | mean(Q78.3, Q78.4)   | mean(Q86.3, Q86.4)   | mean(Q76.3, Q76.4)   | mean(Q74.3, Q74.4)   | 0–10 (higher = better) |
| 2 | `d1_k6`         | 1. Mental & Physical Health | K6 distress (sum of 6 items, recoded 0–4)    | Q68.1–6              | Q75.1–6              | Q65.1–6              | Q65.1–6              | 0–24 (higher = worse) |
| 3 | `d2_hls_avg`    | 2. Happiness & Life Satisfaction | Happiness & life satisfaction (FI #1 + #2 mean) | mean(Q78.1, Q78.2) | mean(Q86.1, Q86.2) | mean(Q76.1, Q76.2) | mean(Q74.1, Q74.2) | 0–10 |
| 4 | `d3_mp_avg`     | 3. Meaning & Purpose        | Meaning & purpose (FI #5 + #14 mean)         | mean(Q78.5, Q78.14)  | mean(Q86.5, Q86.14)  | mean(Q76.5, Q76.14)  | mean(Q74.5, Q74.14)  | 0–10 |
| 5 | `d3_ikigai`     | 3. Meaning & Purpose        | Ikigai (FI #6, standalone)                   | Q78.6                | Q86.6                | Q76.6                | Q74.6                | 0–10 |
| 6 | `d4_cv_avg`     | 4. Character & Virtue       | Character & virtue (FI #8 + #15 mean)        | mean(Q78.8, Q78.15)  | mean(Q86.8, Q86.15)  | mean(Q76.8, Q76.15)  | mean(Q74.8, Q74.15)  | 0–10 |
| 7 | `d5_rel_avg`    | 5. Social Wellbeing         | Relations (FI #9 + #10 mean)                 | mean(Q78.9, Q78.10)  | mean(Q86.9, Q86.10)  | mean(Q76.9, Q76.10)  | mean(Q74.9, Q74.10)  | 0–10 |
| 8 | `d5_lsns_total` | 5. Social Wellbeing         | LSNS-6 total (sum of 6 items minus 6)        | Q18.1–6              | Q18.1–6              | Q17.1–6              | Q20.1–6              | 0–30 (higher = better) |
| 9 | `d5_ucla`       | 5. Social Wellbeing         | UCLA-3 loneliness (sum of 3 items, recoded 0–3) | Q68S1.1–3         | Q76.1–3              | Q66.1–3              | Q66.1–3              | 0–9 (higher = worse) |
| 10 | `d6_worry_avg` | 6. Financial & Material Stability | Financial & safety worry (FI #11 + #12 mean) | mean(Q78.11, Q78.12) | mean(Q86.11, Q86.12) | mean(Q76.11, Q76.12) | mean(Q74.11, Q74.12) | 0–10 (higher = worse) |
| 11 | `d7_sleep`     | 7. Health Behavior          | Sleep hours/day        | Q28.9 | Q31.9 | Q28.9 | Q32.9 | hours/day (band → midpoint; see §3) |
| 12 | `d7_sitting`   | 7. Health Behavior          | Sitting hours/day      | Q28.5 | Q31.5 | Q28.5 | Q32.5 | hours/day (band → midpoint) |
| 13 | `d7_walking`   | 7. Health Behavior          | Walking hours/day      | Q28.6 | Q31.6 | Q28.6 | Q32.6 | hours/day (band → midpoint) |

### 2.1 FI items — stem and scale

**Stem (all four waves, identical wording):** 「以下の質問について0から10段階で当てはまるものを選んでください。」 ("For each of the following statements, please choose the response that best matches you on a 0–10 scale.")

**Choices (all waves):** integer 0, 1, 2, …, 10 (11 levels). Higher = more of the construct.

**Sub-item index across waves (the same `.k` means the same construct in every wave):**

| FI .k | Construct |
|---|---|
| 1  | overall life satisfaction |
| 2  | happiness |
| 3  | physical health |
| 4  | mental health |
| 5  | meaningful daily activities |
| 6  | ikigai (生きがい) |
| 8  | delayed gratification / patience |
| 9  | relationships with friends |
| 10 | relationships with family / intimate |
| 11 | financial worry |
| 12 | safety / crime worry |
| 14 | sense of life purpose |
| 15 | service to society / contribution |

### 2.2 K6 distress (sum of 6 items)

Stem (all waves, identical): 「直近30日間に、どれくらいの頻度で次のことがありましたか。」

Choices (5-pt, all waves): 1 = いつも, 2 = たいてい, 3 = ときどき, 4 = 少しだけ, 5 = まったくない.

Recode: raw → `max(0, min(4, 5 − raw))` (so 0 = まったくない, 4 = いつも). Score = sum of 6 recoded items, range 0–24.

### 2.3 LSNS-6 total (sum of 6 items)

Stem (all waves): 「家族や親戚、近くに住んでいる人を含むあなたの友人全体について考えます。下記の質問に最もあてはまる回答を選んでください。」

Choices (6-pt, all waves): 1 = いない (0人), 2 = 1人, 3 = 2人, 4 = 3–4人, 5 = 5–8人, 6 = 9人以上.

Recode: raw → raw − 1 (so item range 0–5). Score = sum of 6 recoded items, range 0–30.

### 2.4 UCLA-3 loneliness (sum of 3 items)

Stem (all waves): 「直近30日間に、どれくらいの頻度で次のことがありましたか。」

Choices (4-pt, all waves): 1 = 常にある, 2 = 時々ある, 3 = ほとんどない, 4 = 決してない.

Recode: raw → `max(0, min(3, 4 − raw))` (so 0 = 決してない, 3 = 常にある). Score = sum of 3 recoded items, range 0–9.

---

## 3. Time-use items (sleep / sitting / walking / smartphone / PC-tablet)

The time-use grid question changed scale in 2025: earlier waves use a 12-level frequency band; 2025 splits the upper bands and uses 16 levels. The script first re-bins 2025 onto the 12-level scale, then converts the band to midpoint hours/day.

### 3.1 Q-code shifts

| Construct | 2022 | 2023 | 2024 | 2025 |
|---|---|---|---|---|
| Sleep             | Q28.9  | Q31.9  | Q28.9  | Q32.9  |
| Sitting           | Q28.5  | Q31.5  | Q28.5  | Q32.5  |
| Walking           | Q28.6  | Q31.6  | Q28.6  | Q32.6  |
| Smartphone time   | Q28.13 | Q31.13 | Q28.13 | Q32.11 |
| PC / tablet time  | Q28.14 | Q31.14 | Q28.14 | Q32.12 |

(Smartphone and PC/tablet time are used as **confounders** at the tier's baseline wave only, not as outcomes.)

### 3.2 Raw band scales

**2022 / 2023 / 2024 (12 levels):**
1. なし (0時間)
2. 30分未満
3. 30分程度
4. 1時間
5. 2時間
6. 3時間
7. 4–5時間
8. 6–7時間
9. 8–9時間
10. 10–11時間
11. 12時間以上
12. わからない

**2025 (16 levels, finer-grained for 4–11 h):**
codes 1–6 align with the earlier scale; 7 = 4h, 8 = 5h, 9 = 6h, 10 = 7h, 11 = 8h, 12 = 9h, 13 = 10h, 14 = 11h, 15 = 12時間以上, 16 = わからない.

### 3.3 Harmonisation pipeline

`recode_timeuse_2025()` collapses 2025 levels back onto the 12-level scale (`7+8 → 7`, `9+10 → 8`, `11+12 → 9`, `13+14 → 10`, `15 → 11`, `16 → 12`).

`band_to_hours()` then maps the 12-level band to numeric hours/day midpoints:

| band | hours/day |
|---:|---:|
| 1 (なし)        | 0    |
| 2 (<30 min)     | 0.25 |
| 3 (30 min)      | 0.5  |
| 4 (1 h)         | 1    |
| 5 (2 h)         | 2    |
| 6 (3 h)         | 3    |
| 7 (4–5 h)       | 4.5  |
| 8 (6–7 h)       | 6.5  |
| 9 (8–9 h)       | 8.5  |
| 10 (10–11 h)    | 10.5 |
| 11 (12 h+)      | 12   |
| 12 (わからない) | NA   |

### 3.4 わからない handling

- As an **outcome**: rows where Y_pre or Y_post is わからない drop from the time-use regression (`is.finite(y_pre) & is.finite(y_post)` filter).
- As a **covariate** (other-outcome adjustment for non-time-use focal outcomes): NA → 0 via `fill_zero_pre()`.
- In **descriptive tables**: わからない is excluded from the displayed mean (the supplement notes how many were excluded).

---

## 4. Confounders (9 variables, drawn from the tier's baseline wave)

The same construct is taken from the appropriate baseline wave for each tier:

| Tier         | Baseline wave |
|---|---|
| Main, Sens 2 | 2024 |
| Sens 1       | 2023 |
| Sens 3       | 2022 |

### 4.1 Sex (binary)

| Wave | Column      | Choices                       | Coding                              |
|---|---|---|---|
| 2022 | `SEX_2022`  | 1 = 男性, 2 = 女性             | `is_female = 1` if `SEX == 2` else 0 |
| 2023 | `SEX_2023`  | (same)                        | (same)                              |
| 2024 | `SEX_2024`  | (same)                        | (same)                              |

### 4.2 Age (9 dummies + reference ≥ 65)

`AGE_<wave>` is numeric (years). Nine dummy variables are created for the 5-year bands 18–24, 25–29, 30–34, 35–39, 40–44, 45–49, 50–54, 55–59, 60–64; reference = ≥ 65.

### 4.3 Education (2 dummies + reference)

Stem (all waves): 「Q21.学歴についてお答えください。それぞれについて最後に卒業された、または在学中の学校を教えてください。」

| Wave | Q-code         | Choices |
|---|---|---|
| 2022 | `Q21.1_2022`   | 1=中学校, 2=私立高校, 3=国立・公立高校, 4=専門学校, 5=短大・高専, 6=私立大学, 7=国立大学, 8=公立大学, 9=大学院, 10=その他, 11=分からない |
| 2023 | `Q21.1_2023`   | (same 11 levels) |
| 2024 | `Q21.1_2024`   | (same) |

Coding: `edu_univ = 1` if `Q21.1 ∈ {6, 7, 8}`; `edu_grad = 1` if `== 9`; reference ("less than university") = codes 1–5 + Unknown (10, 11) + NA.

### 4.4 Employment (5 dummies + reference) — ⚠️ Q-code differs in 2023

Stem (all waves): 「仕事（休業中の仕事も含む）の状況にあてはまるものを1つ選んでください。」

| Wave | Q-code             | Choices |
|---|---|---|
| 2022 | `Q5.1_2022`        | 1=会社の役員, 2=自営業主, 3=フリーランス, 4=自家営業の手伝い, 5=正社員（管理職）, 6=正社員（管理職以外）, 7=派遣社員, 8=契約社員・嘱託, 9=アルバイト・パート, 10=プラットフォーム単発, 11=内職, 12=アルバイト学生, 13=仕事をしていない学生, 14=リタイア, 15=専業主婦・主夫, 16=無職 |
| 2023 | **`Q6.1_2023`** ⚠️ | (same 16 levels). NOTE: `Q5.1_2023` is a vaccine-attitudes question, not employment. |
| 2024 | `Q5.1_2024`        | (same 16 levels) |

Coding: `emp_exec = 1` if `== 1`; `emp_self = 1` if `∈ {2, 3, 4}`; `emp_nonreg = 1` if `∈ {7, 8, 9, 10, 11}`; `emp_student = 1` if `∈ {12, 13}`; `emp_notwork = 1` if `∈ {14, 15, 16}` or NA; reference = Regular employee (codes 5, 6).

### 4.5 Household income (4 dummies + Unknown + reference) — ⚠️ Q-code differs in 2023 and 2024

Stem (all waves): 「あなたの下記の金額はおおよそどのくらいですか。」

| Wave | Q-code               | Choices |
|---|---|---|
| 2022 | `Q87.1_2022`         | 20-band ordinal: 1 = 0円, 2 = 50万円未満, …, 18 = 2,000万円以上, 19 = 答えたくない, 20 = 分からない |
| 2023 | **`Q90.1_2023`** ⚠️ | (same 20-level scale). NOTE: `Q87.1_2023` is **adverse childhood experiences**, not income. |
| 2024 | **`Q80.1_2024`** ⚠️ | (same 20-level scale). |

Coding (fixed yen-band cuts):
- `inc_lt2m`   (reference) = codes 1–4 (< 2,000,000 yen)
- `inc_2_6m`   = codes 5–8 (2M to < 6M)
- `inc_6_10m`  = codes 9–12 (6M to < 10M)
- `inc_10m_plus` = codes 13–18 (≥ 10M)
- `inc_unknown` = codes 19, 20, or NA

### 4.6 Marital (binary)

Stem (all waves): 「現在、配偶者（夫または妻）は、いますか。配偶者には、事実上夫婦として生活しているが、婚姻届を提出していない場合や同性のパートナーも含みます。」

| Wave | Q-code     | Coding |
|---|---|---|
| 2022 | `Q2_2022`  | `married = 1` if `Q2 ∈ {1, 2, 3}` else 0 |
| 2023 | `Q2_2023`  | (same) |
| 2024 | `Q2_2024`  | (same) |

### 4.7 Living alone (binary)

Stem (all waves): 「ふだん一緒にお住まいで、かつ、生計を共にしている方（世帯員）は、あなたを含めて何人ですか。」

| Wave | Q-code        | Choices         | Coding |
|---|---|---|---|
| 2022 | `Q1.1_2022`   | numeric (人数)  | `living_alone = 1` if `Q1.1 == 1` else 0 |
| 2023 | `Q1.1_2023`   | (same)          | (same) |
| 2024 | `Q1.1_2024`   | (same)          | (same) |

### 4.8 Smartphone time (continuous, hours/day)

| Wave | Q-code        | Raw choices |
|---|---|---|
| 2022 | `Q28.13_2022` | 12-level band (see §3.2) |
| 2023 | `Q31.13_2023` | (same) |
| 2024 | `Q28.13_2024` | (same) |

Coding: `band_to_hours()` (§3.3). わからない (band 12) → NA, then `fill_zero_pre()` replaces NA with 0 for analytic use.

### 4.9 PC / tablet time (continuous, hours/day)

| Wave | Q-code        | Raw choices |
|---|---|---|
| 2022 | `Q28.14_2022` | (same 12-level band) |
| 2023 | `Q31.14_2023` | (same) |
| 2024 | `Q28.14_2024` | (same) |

Coding: identical to §4.8.

### 4.10 Pre-outcome adjustment (12 other-outcome levels at the baseline wave)

For the focal outcome `k`, the covariate vector for tier `t` includes:
- the 12 other-outcome levels (i.e., the 13 outcomes minus `k`), materialised at tier `t`'s baseline wave;
- for **Sens 2** and **Sens 3** (ANCOVA), the focal outcome's own baseline `Y_k_<pre>` enters as a 13th covariate (with its slope estimated freely);
- for **Main** and **Sens 1** (first-difference DiD), the focal `Y_k_<pre>` is differenced out into ΔY and is not entered as a covariate.

Pre-outcome NA → 0 (consistent with the JACSIS mandatory-item convention: all outcome items are mandatory at every wave, so genuine NAs are rare).

---

## 5. Statistical specification (summary)

For full detail, see the manuscript Methods. The script's regressions are:

- **Main / Sens 1 (DiD):**  ΔY_k = α + τ · exposure + Y_pre,(−k)′η + X_pre′β + ε.
- **Sens 2 / Sens 3 (ANCOVA):**  Y_post,k = α + τ · exposure + ρ · Y_pre,k + Y_pre,(−k)′η + X_pre′β + ε.

with:
- 13 outcomes × 5 exposures fit per tier (Pattern A0 fits `use`; Pattern A fits `intensity`; Pattern B jointly fits `creative + daily + social`);
- OLS estimation with HC1 robust standard errors (`sandwich::vcovHC(., type = "HC1")`);
- per-exposure Benjamini–Hochberg FDR at q < 0.10 across the 13-outcome family;
- outcome-standardised coefficients reported alongside raw: `estimate_oz = estimate / sd_pre`, where `sd_pre` is the baseline-wave SD of the focal outcome on the tier's analytic cohort;
- E-values (VanderWeele & Ding 2017) computed for every results row via the Chinn (2000) `RR ≈ exp(0.91 × |d|)` approximation.

### Users-only sensitivity (Pattern A + Pattern B restricted to AI users)

In addition to the main tier-level fits, each tier also fits Pattern A (`intensity`) and Pattern B (`creative + daily + social` jointly) restricted to the tier-specific user subset. This captures within-user dose–response without the user-vs-non-user contrast bundled into the full-sample fits. The outcome-standardisation `sd_pre` for these models is computed on the user subset only.

---

## 6. Output column glossary (per-tier `results_<tier>_<phase>.csv`)

Each tier writes a combined results CSV with 65 rows (13 outcomes × 5 exposures) — except the users-only sensitivities, which have 52 rows (13 outcomes × 4 exposures; Pattern A0 is degenerate among users).

| Column | Meaning |
|---|---|
| `outcome`         | Outcome variable id (`d1_pm_avg`, `d1_k6`, …, `d7_walking`). |
| `exposure`        | One of `use`, `intensity`, `creative`, `daily`, `social`. |
| `estimate`        | Regression coefficient on the raw outcome scale. |
| `se`              | HC1 robust standard error. |
| `z`, `p`          | Wald z-statistic and two-sided p-value. |
| `ci_lo`, `ci_hi`  | 95% CI on the raw scale (`estimate ± 1.96 × se`). |
| `n`               | Analytic sample size for the row. |
| `kind`            | `"continuous"` for all 13 outcomes in this release. |
| `family`          | `A0_<phase>` / `A_<phase>` / `B_<phase>` (with `<phase>` = `did`, `pretrends`, `ancova`). |
| `p_bh`            | Benjamini–Hochberg q-value within the per-exposure 13-test family. |
| `bh_sig`          | `TRUE` if `p_bh < 0.10`. |
| `sd_pre`          | Baseline-wave SD of the focal outcome on the tier's analytic cohort. |
| `estimate_oz`     | Outcome-standardised coefficient = `estimate / sd_pre`. |
| `se_oz`, `ci_lo_oz`, `ci_hi_oz` | Same rescaling applied to SE and CI bounds. |
| `e_value`         | E-value (VanderWeele & Ding 2017) at the point estimate. |
| `e_value_ci`      | E-value at the CI bound closer to the null; = 1.00 when the 95% CI crosses zero. |
