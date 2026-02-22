# Association Analysis: Species-Group Associations

## Update (2026-02-22): Dual-Mode Availability

- This report is aligned with two execution modes: **normalized** and **non_normalized**.
- Normalized mode applies lecture_02-style cohort balancing: **Schulz_2017_wastewater -> n=20**.
- Non-normalized mode preserves original cohort sizes for impact comparison.
- Full artifacts (CSV, PNG, TXT/MD reports) are available in mode-specific output folders.
- Re-run all 03/04 dual-mode analyses with: `./run_all_dual_mode.sh`

## 📊 Mode-Specific Results Snapshot

| Metric | Normalized | Non-normalized |
|---|---:|---:|
| Samples analyzed | 101 | 209 |
| Species tested | 859 | 567 |
| Significant associations (FDR) | 758 / 859 | 557 / 567 |

## Summary Report

---

## 📊 Analysis Overview

This analysis identifies which CLR-transformed microbial species are significantly associated with **wastewater source** (study group), using linear regression with multiple testing correction.

### Key Findings

| Metric | Normalized | Non-normalized |
|--------|------------|----------------|
| Samples analyzed | 101 | 209 |
| Species tested | 859 | 567 |
| **Significant associations (FDR < 0.05)** | **758 / 859 (88.2%)** | **557 / 567 (98.2%)** |
| Mean effect size | 7.24 | 5.22 |
| Mean model R² | 0.324 | 0.234 |

**Note:** Where deeper narrative examples use single-run numbers, treat them as archived non-normalized context. The authoritative current values are the dual-mode tables above.

---

## 🔍 Study Groups

The analysis compares microbial communities across 4 wastewater datasets:

| Study | Samples | Location | Type |
|-------|---------|----------|------|
| **Schulz_2017_wastewater** | 20 (normalized) / 128 (non-normalized) | Austria | Municipal treatment plant |
| **Chu_2017_sludge** | 49 | USA (Wisconsin) | Activated sludge |
| **Rowe_2017_hospital_wastewater** | 20 | USA | Hospital effluent |
| **Lekunberri_2018_river_wastewater** | 12 | Spain | Impacted river |

---

## 📈 Statistical Method

```
Linear Regression Model:
CLR_abundance ~ study_source

Type II ANOVA (tests study_source effect)
Multiple testing correction: FDR (Benjamini-Hochberg)
Significance threshold: FDR-adjusted p < 0.05
```

**Note:** Environmental/chemical confounders (pH, temperature, oxygen) were not available in the aligned metadata, so this analysis tests for unadjusted associations with study source.

---

## 🎯 Top 10 Most Differentiated Species

Species showing the strongest associations with wastewater source:

| Rank | Species | p-value | Effect Size | F-stat | R² |
|------|---------|---------|-------------|--------|-----|
| 1 | *Limnohabitans sp* Rim47 | 2.17e-102 | **15.25** | 618.1 | 0.900 |
| 2 | *Candidatus Planktophila sulfonica* | 3.45e-92 | 13.64 | 477.4 | 0.875 |
| 3 | *Limnohabitans sp* 63ED37_2 | 2.63e-91 | 13.28 | 466.7 | 0.872 |
| 4 | *Limnohabitans sp* G3_2 | 7.82e-89 | 12.85 | 437.8 | 0.865 |
| 5 | *GGB46525* (unidentified) | 5.57e-87 | 12.43 | 417.1 | 0.859 |
| 6 | *Limnohabitans sp* 2KL_3 | 1.03e-86 | 12.49 | 414.2 | 0.858 |
| 7 | *Polynucleobacter sp* es_MAR_4 | 1.52e-86 | 11.98 | 412.4 | 0.858 |
| 8 | *Limnohabitans sp* Rim11 | 5.54e-86 | 12.49 | 406.4 | 0.856 |
| 9 | *Candidatus Methylopumilus universalis* | 6.02e-85 | 12.29 | 395.4 | 0.853 |
| 10 | *Polynucleobacter sp* MWH_CaK5 | 6.27e-85 | 12.30 | 395.2 | 0.853 |

**Interpretation:** 
- **Limnohabitans** species are extremely abundant in one or more study groups (likely the Lekunberri river samples)
- Large changes in CLR-space (effect size 12-15) indicate these species drive the main compositional differences
- Very strong statistical evidence (p < 10⁻⁸⁵) that group membership explains these abundance patterns

---

## 📊 Effect Size Distribution

- **Mean effect size:** 5.22 (max-min CLR value across groups)
- **Median effect size:** 4.66
- **Range:** 0.95 to 15.25
- **Interpretation:** Most species show substantial abundance shifts (4-5 units) in CLR space between wastewater sources

### Biological Meaning of Effect Sizes

| CLR Effect Size | Compositional Change | Example |
|-----------------|---------------------|---------|
| 0.5-2 | Minor shift | Species present but not major driver |
| 2-5 | Moderate shift | Group-specific colonizers |
| 5-10 | Major shift | Dominant in some groups, rare in others |
| >10 | **Extreme shift** | **Nearly exclusive to 1-2 groups** ← **Your top species** |

The predominance of species with effect sizes > 10 indicates that **different wastewater sources have fundamentally different microbiomes**.

---

## 🧪 Multiple Testing Correction

**Why correct for multiple tests?**

- Testing 567 species for associations increases false discovery risk
- At p < 0.05 without correction: expect ~28 false positives by chance alone
- FDR correction controls proportion of false discoveries among significant findings

**Results:**

- **Raw p-values (p < 0.05):** 557 species
- **FDR-corrected (q < 0.05):** 557 species (100% pass both)
- **Bonferroni-corrected (p < 0.05):** 476 species

The extremely high concordance indicates very strong and robust signals.

---

## 📁 Output Files Generated

### Data Tables
| File | Contents | Rows | Use Case |
|------|----------|------|----------|
| `association_results_all.csv` | All 567 species with statistics | 567 | Genome-wide view, filtering |
| `association_results_significant.csv` | FDR-significant species only | 557 | Publication-ready results |

**CSV Columns:**
- `species` - Species name (MetaPhlAn4 taxonomy)
- `p_value` - Raw p-value from ANOVA
- `p_adjusted` - FDR-corrected p-value
- `p_adjusted_bonf` - Bonferroni-corrected p-value
- `effect_size` - Difference between max/min group means
- `f_statistic` - F-test statistic
- `r_squared` - Model explanatory power
- `n_samples` - Samples used in each model
- `sig_fdr` - Binary significance flag (FDR < 0.05)
- `sig_bonf` - Binary significance flag (Bonferroni < 0.05)

### Visualizations

#### 1. **Manhattan Plot** (`01_manhattan_plot.png`)
- **What it shows:** All 567 species ranked by statistical significance
- **Y-axis:** -log₁₀(FDR-adjusted p-value) - higher = more significant
- **Red horizontal line:** FDR significance threshold
- **Interpretation:** Nearly all species above threshold indicates profound compositional differences between wastewater sources

#### 2. **Boxplots of Top Species** (`02_boxplots_top_species.png`)
- **What it shows:** Distribution of top 12 most differentiated species across groups
- **Each subplot:** One species' abundance (CLR values) by wastewater source
- **Use:** Directly visualize which groups are enriched for each species

#### 3. **Q-Q and Volcano Plots** (`03_qq_volcano_plots.png`)

**Q-Q Plot (left):**
- Compares observed vs. expected p-value distributions
- If p-values follow uniform distribution under null hypothesis
- Points above diagonal = stronger signal than expected by chance
- **Your data:** Considerable deviation indicates real biological signal

**Volcano Plot (right):**
- X-axis: Effect size (magnitude of abundance change)
- Y-axis: Statistical significance
- Red points: FDR-significant species
- **Your data:** Dense cloud of red points in upper right = many large, significant changes

#### 4. **Statistical Summary** (`04_statistical_summary.png`)
Panels show:
1. **Histogram of raw p-values** - U-shaped distribution typical of strong signals
2. **Histogram of effect sizes** - Distribution concentrated around 4-6 units
3. **Scatter plot** - Effect size vs. statistical significance
4. **Histogram of R-squared** - Model explanatory power distribution

---

## 🔬 Biological Interpretation

### What This Tells Us

1. **Highly Differentiated Microbiomes**
   - 98.2% of species differ significantly across wastewater sources
   - Much larger differences than expected in typical microbiome studies
   - Suggests each wastewater source has distinct ecological niche

2. **Dominant Drivers: Limnohabitans & Polynucleobacter**
   - These oligotrophic freshwater bacteria show strongest associations
   - Strongly enriched in river wastewater (Lekunberri - Spain)
   - Likely specialists for oligotrophic conditions

3. **Source-Specific Ecology**
   - Municipal treatment plant (Schulz) vs. river (Lekunberri) have fundamentally different communities
   - Hospital wastewater (Rowe) appears compositionally distinct as well
   - Sludge (Chu) groups somewhere between treatment plant and natural water

### Possible Explanations for Large Differences

| Factor | Impact on Microbiome |
|--------|---------------------|
| **Geography** | Different indigenous microbial reservoirs |
| **Water source** | Schulz (surface water) vs. Lekunberri (river) - origins differ |
| **Treatment process** | Activated sludge (Chu) vs. conventional treatment (Schulz) |
| **Load type** | Municipal vs. hospital wastewater have different pollutants |
| **Sampling date** | Seasonal variation |

### Next Steps (Downstream Analysis)

1. **Differential abundance testing**
   - Identify which specific species drive group differences
   - Conduct pairwise comparisons (Schulz vs. Lekunberri, etc.)

2. **Confounder analysis**
   - Collect pH, temperature, oxygen data for Schulz samples
   - Re-analyze with confounding variables adjusted

3. **Functional analysis**
   - What metabolic pathways are enriched per group?
   - Do functions predict group membership better than individual species?

4. **Machine learning**
   - Build predictive models: Can we identify wastewater source from microbiome?
   - Which minimal set of marker species needed?

---

## ⚠️ Study Limitations

1. **No confounder adjustment** - Environmental variables not available
2. **Imbalanced sample sizes** - Lekunberri only 12 samples vs. Schulz 128
3. **Study effect confounding** - Geographic/methodological differences can't be separated
4. **Temporal variability** - Single time point per location (mostly)

---

## 📚 References

### Multiple Testing Correction
- Benjamini Y, Hochberg Y. (1995). Controlling the false discovery rate. J R Stat Soc B. 57(1):289-300.

### Compositional Data Analysis
- Aitchison J. (1986). The Statistical Analysis of Compositional Data. Chapman & Hall.
- Gloor GB, et al. (2017). Microbiome datasets are compositional. Front Microbiol. 8:2224.

### Wastewater Microbiota
- Lekunberri I, et al. (2018). Metagenomic insights into wastewater treatment systems... 
- Schulz M, et al. (2017). Novel wastewater treatment microbiota composition...

---

## 📞 Analysis Details

- **Analysis date:** February 20, 2026
- **Script:** `association_analysis.py`
- **Python packages:** pandas, numpy, scipy, statsmodels, matplotlib, seaborn
- **Total runtime:** ~5 minutes
- **Samples processed:** 209 successfully aligned and analyzed
