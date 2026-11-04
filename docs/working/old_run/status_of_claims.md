# Status of each claim (earlier validation)

From [`0_overview.md`](https://github.com/dm-k1/objective-bayes/blob/main/docs/0_overview.md) as it was before the redesign. The check
numbers (P1, A4, C7, ...) are those of the earlier studies, which remain on
the `main` branch; [`C_study_design.md`](../../appendix/C_study_design.md) maps each to its new number.
This table will be rewritten from the redesigned validation.


Two labels are used. **Checked** means a numerical check passes on the current
code, in the settings the study covers; no claim is proven analytically here.
**Run** means a simulation study was run on the current code, and its result
is given. Details are in [`validation_report.md`](validation_report.md). All results are
from the full run of 4 November 2026 (commit 21d42dc), in which 137 of 138
checks pass. The one failure is a finding about flat R-squared (validation 10,
P4; see [`comparison.md`](comparison.md), section 5).

| Claim | Flat partial r | Flat R | Flat R-squared | Where |
|---|---|---|---|---|
| Valid correlation matrices give the three versions | checked (F1) | checked (F2) | checked (F3) | validation 00 |
| The Stan model samples the prior correctly | checked (01-P5–P7, KS ≤ .012) | checked (01-P1–P4, KS ≤ .009) | checked (01-P8–P10, KS ≤ .009) | validation 01, P1–P10 |
| The Stan posterior is calibrated ($k = 8$) | checked (SBC *p* ≥ .02) | checked (SBC *p* ≥ .38) | checked (SBC *p* ≥ .22) | validation 01, S |
| The whole-model Bayes factor is exact | checked (A1–A4) | checked (A1–A4) | checked (A1–A4) | validation 02 |
| Units and recoding don't change the answer | checked (A2) | checked (A2) | checked (A2) | validation 02 |
| Predictive matching (BF = 1 at $n = k + 1$) | checked (B1) | checked (B1) | checked (B1) | validation 02 |
| Information consistency, for $n \ge k + 3$ | checked (B2, B2b) | checked (B2, B2b) | checked (B2, B2b) | validation 02 |
| Intrinsic consistency | checked (B5) | checked (B5) | checked (B5) | validation 02 |
| Model-selection consistency, fixed $k$ | checked: $\tfrac q2 \log n$ in every case | checked: ½ log n against the intercept-only model, $\tfrac q2 \log n$ otherwise | checked: $\log n$ against the intercept-only model, $\tfrac q2 \log n$ otherwise (slopes 1.01, 0.96, 0.99 and 0.55, 1.55, 2.95) | validation 02, B3, B3b, B4 |
| Model-selection consistency, $k$ growing with $n$ | checked (G1, G2): evidence grows fast (log BF₀₁ 408 at $n = 2000$, $k = n/4$) | checked, but slow under the null (log BF₀₁ 1.9) | checked, but slow under the null (log BF₀₁ 3.7) | validation 08 |
| What singles the prior out | checked (C7) | checked (V5) | – | validation 05, 09 |
| The conditional test is computed correctly | checked (C1, C3–C5) | checked (C1, C2, C5) | checked (C1, C2, C5) | validation 05, C1–C5 |
| The block conditional test is computed correctly | checked (C9–C11) | checked (C9–C11) | checked (C9–C11) | validation 05, C9–C11 |
| Separating nulls from small effects, against JZS | run: equal or better (AUC +.001 to +.045, mean over $n$ and pattern per $k$ and correlation) | run: the same as flat partial r (conditional) | run: the same as flat partial r (conditional) | validation 04, 06 |
| Confirming a true null, one coefficient | run: 86–88%, for any $k$ | run: 55% at $k = 2$, 6% at $k = 20$ | run: 86% at $k = 2$, 26% at $k = 20$ | validation 06 |
| Calling a small effect (0.1) null | run: 52–60% | run: 19% at $k = 2$, 0.5% at $k = 20$ | run: 52% at $k = 2$, 4% at $k = 20$ | validation 06 |
| Calibrated on its own prior | checked (P1) | checked (P1) | checked (P1) | validation 10 |
| Variable selection | run: fewest false inclusions of the three, misses more (mean 0.36 of 8 nulls included, 2.09 of 4 real missed) | run: many more inclusions, fewer misses (1.14, 1.55) | run: nearly the same as flat R (1.13, 1.56) | validation 09 |
| Estimation: prediction error, interval coverage | run: better than OLS everywhere, coverage ≥ .916 | run: lowest error at small $n$, coverage down to .76 | run: in between, coverage down to .82 | validation 07 |

---

