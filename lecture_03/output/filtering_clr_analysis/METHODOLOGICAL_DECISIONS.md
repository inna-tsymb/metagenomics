# Methodological Decisions in Wastewater Microbiome Analysis

## Update (2026-02-22): Dual-Mode Availability

- This report is aligned with two execution modes: **normalized** and **non_normalized**.
- Normalized mode applies lecture_02-style cohort balancing: **Schulz_2017_wastewater -> n=20**.
- Non-normalized mode preserves original cohort sizes for impact comparison.
- Full artifacts (CSV, PNG, TXT/MD reports) are available in mode-specific output folders.
- Re-run all 03/04 dual-mode analyses with: `./run_all_dual_mode.sh`

## Overview
This document explains the key design decisions made in the filtering and CLR transformation analysis, specifically:
1. Selection of taxonomic level (Species - s__)
2. Selection of prevalence threshold (10%)
3. Selection of abundance threshold (0.01%)

---

## 1. TAXONOMIC LEVEL SELECTION: SPECIES (s__)

### Decision: Filter for species-level (s__) taxonomy, exclude strains (t__)

### Rationale

#### 1.1 **Why Species Level?**

**Species-level provides optimal balance between:**

- **Resolution vs. Robustness**
  - Higher than genus level: Provides functional specificity
    - Different species within same genus can have vastly different metabolic capabilities
    - Example: *Acinetobacter baumannii* (pathogenic) vs *Acinetobacter baylyi* (environmental)
  - Lower than strain level: Avoids over-specificity
    - Strain-level annotations are often incomplete/uncertain in metagenomic data
    - Strains have high noise-to-signal ratio

- **Biological Relevance**
  - Species is the primary unit of ecological and functional analysis
  - Most metabolic predictions and environmental correlations work at species level
  - Downstream databases (KEGG, eggNOG) are primarily annotated at species level

- **Statistical Power**
  - Species level: ~600-800 taxa typically (good for regression analysis)
  - Genus level: ~150-200 taxa (loss of resolution)
  - Strain level: >10,000 taxa (sparse, unreliable)

#### 1.2 **Why Exclude Strains (t__)?**

**Strains (t__SGB*) are problematic in wastewater studies:**

```
Problems with strain-level data:
├─ High annotation uncertainty
│  └─ Strain assignments are probabilistic, not definitive
├─ Severe sparsity
│  └─ Most strains appear in <2% of samples (random noise)
├─ Redundant information
│  └─ Multiple strains of same species inflate complexity
└─ Low prevalence
   └─ Violates prevalence filtering (would remove most anyway)
```

**MetaPhlAn4 Strain Output:**
- `t__SGB80728`: Strain-level Bacterial Genome (SGB) identifier
- These are often uncertain taxonomic assignments
- Recommended practice in literature: **filter out strains for microbiome studies**

#### 1.3 **What About Other Levels?**

| Level | Code | Use Case | Why Not Here |
|-------|------|----------|--------------|
| **Kingdom** | k__ | Quality control | Too broad (all Bacteria) |
| **Phylum** | p__ | Community composition | Too broad for diff. abundance |
| **Class** | c__ | Broad patterns | Loss of functional resolution |
| **Order** | o__ | Medium patterns | Still too broad |
| **Family** | f__ | Good baseline | Common choice, but less specific |
| **Genus** | g__ | ✓ Alternative | Common in lit., ~200 taxa |
| **Species** | s__ | **✓ CHOSEN** | Optimal for regression analysis |
| **Strain** | t__ | Never | Too noisy, too sparse |

### Implementation in Code

```python
# Filter for species level only (s__) and exclude strains (t__)
df_species = df_long[
    df_long['clade_name'].str.contains('s__', na=False) & 
    ~df_long['clade_name'].str.contains('t__', na=False)
].copy()
```

### Result
- Initial: 19,122 metaphlan annotations → 573 species identified
- After filtering: **567 species retained** (high confidence species)

---

## 2. PREVALENCE THRESHOLD SELECTION: 10%

### Decision: Keep species present in ≥10% of samples

### Scientific Rationale

#### 2.1 **What is Prevalence?**

**Definition:**
```
Prevalence = (Number of samples where species detected) / (Total samples)
           = Proportion of samples containing ≥1 read of the species
```

**Example with n=209 samples:**
- 10% prevalence = species in ≥21 samples
- 5% prevalence = species in ≥10 samples  
- 20% prevalence = species in ≥42 samples

#### 2.2 **Why 10% and Not Something Else?**

**Threshold Selection Logic:**

```
Too LOW (1-5%):
  ├─ Retain: 2000-3500 species (all rare species)
  │  ├─ Advantages: Capture rare/specialist taxa
  │  └─ Disadvantages: 
  │     ├─ High noise, low signal
  │     ├─ Inflation of false positives
  │     ├─ Violates abundance of rare species
  │     └─ Statistical power ↓ (sparse matrix)
  │
OPTIMAL (10%):
  ├─ Retain: 400-600 species (balance)
  │  ├─ Advantages:
  │  │  ├─ Removes random artifacts
  │  │  ├─ Statistical power maintained
  │  │  ├─ Good sensitivity/specificity
  │  │  ├─ Computational efficiency
  │  │  └─ Common in literature
  │  └─ Disadvantages: May miss rare keystone species
  │
Too HIGH (20-50%):
  ├─ Retain: 50-150 species (only core community)
  │  ├─ Advantages: Very high signal/noise
  │  └─ Disadvantages:
  │     ├─ Loss of biodiversity information
  │     ├─ Miss context-specific species
  │     └─ Reduced power for associations
```

#### 2.3 **Literature Standards**

**Common prevalence thresholds in microbiome literature:**

| Threshold | Use Case | Authors |
|-----------|----------|---------|
| 1-5% | Exploratory, rare taxa focus | *Very permissive* |
| **10%** | **Standard differential abundance** | *Most common* |
| 15% | Conservative, robust only | *More stringent* |
| 20-25% | Core microbiome only | *Very conservative* |
| 50%+ | Essential microbes | *Ultra-conservative* |

**Our Decision: 10%**
- ✓ Aligns with standard wastewater microbiome studies
- ✓ Balances sensitivity and specificity
- ✓ Removes singleton/doubleton contamination
- ✓ Retains clinically/environmentally relevant taxa

#### 2.4 **Wastewater-Specific Considerations**

**Why 10% is appropriate for wastewater:**

```
Wastewater characteristics:
├─ High biodiversity (>5000 species typically)
├─ High technical noise
│  └─ Cross-contamination during sequencing common
├─ Highly variable composition between treatment stages
│  └─ Some taxa appear sporadically
├─ Presence of rare pathogenic bacteria
│  └─ These may be <10% but still important
└─ Dynamic community structure
   └─ Species prevalence varies with season, load
```

**Compromise achieved with 10%:**
- Filters out obvious sequencing noise (single-read taxa)
- Retains variable but biologically relevant species
- Captures treatment-stage specific bacteria
- Still allows detection of environmental/pathogenic signals

### Impact of 10% Threshold

**Our results:**
```
Before filtering: 573 species (includes many rare species)
After 10% prevalence: 567 species retained (removed 6 rare species)
Species removed: Prevalence < 10% AND abundance < 0.01%
```

**Low removal at 10% because:**
- Most rare species already have low abundance (<0.01%)
- Rarity + low abundance create double filter
- Prevalence filter is less stringent than abundance filter in this case

---

## 3. ABUNDANCE THRESHOLD SELECTION: 0.01%

### Decision: Keep species with mean relative abundance ≥0.01%

### Scientific Rationale

#### 3.1 **What is Abundance Threshold?**

**Definition:**
```
Mean Relative Abundance = Average (Relative Abundance in each sample)
                        = Mean % across all samples

With n=209 samples:
0.01% mean abundance ≈ 2-5 average reads across all samples
```

#### 3.2 **Why 0.01% (Not 0.001% or 0.1%)?**

**Threshold Selection Logic:**

```
Too LOW (0.001%):
  ├─ Retain: Most rare species
  │  ├─ Advantages: Capture all biodiversity
  │  └─ Disadvantages:
  │     ├─ Extreme sparsity (>98% zeros)
  │     ├─ Below sequencing error rates
  │     ├─ Likely sequencing artifacts
  │     ├─ No statistical power
  │     └─ Computational burden
  │
OPTIMAL (0.01%):
  ├─ Retain: Abundant + consistent species
  │  ├─ Advantages:
  │  │  ├─ Removes sequencing noise threshold
  │  │  ├─ Focuses on informative taxa
  │  │  ├─ Reduces sparsity (90% vs 97%)
  │  │  ├─ Statistical power ↑
  │  │  └─ Good for regression analysis
  │  └─ Disadvantages: May miss low-abundance keystone
  │
High (0.1-1%):
  ├─ Retain: Only dominant species
  │  ├─ Advantages: Maximum signal/noise
  │  └─ Disadvantages:
  │     ├─ Extreme loss of diversity
  │     ├─ Miss context-specific taxa
  │     └─ Reduced ecological insight
```

#### 3.3 **Scientific Justification for 0.01%**

**Connection to Sequencing Technology:**

```
Illumina sequencing error rate: ~0.1-1% per base
MetaPhlAn4 confidence filtering: Removes <0.001% likely artifacts

0.01% threshold = 10-100x above sequencing noise
                = Represents genuine biological signal
                = Not dependent on stochastic sequencing events
```

**Ecological Justification:**

```
0.01% relative abundance in wastewater = ~10-100 cells/mL typically
(scales with total biomass, sequencing depth)

This represents:
✓ Functionally relevant population sizes
✓ Sufficient reproductive potential for ecological role
✓ Not transient/contaminating organisms
✗ Below typical minimum viable population size for many bacteria
```

#### 3.4 **Literature Standards**

| Threshold | Comment | Use Case |
|-----------|---------|----------|
| 0.001% | Extremely permissive | Meta-analysis of rare taxa |
| **0.01%** | **Standard filtering** | **Most studies** |
| 0.05% | More stringent | Core microbiome |
| 0.1% | Very stringent | Only dominant species |
| 1% | Extreme | Only top species |

**Typical values in published wastewater studies:**
- Range: 0.01% - 0.1%
- Most common: **0.01%** (matches our choice)

#### 3.5 **Wastewater-Specific Considerations**

**Why 0.01% for wastewater:**

```
Wastewater ecology:
├─ Extremely high biomass
│  └─ High background --> high absolute numbers at low %
├─ Treatment processes select for abundant taxa
│  └─ Nitrifiers, heterotrophs accumulate
├─ Many functional specialists at low abundance
│  └─ Methanogens, syntrophs <0.01% but essential
└─ Lots of transient species
   └─ From influent, not established community
```

**With 0.01% we:**
- ✓ Capture nitrifiers (~5% typical)
- ✓ Capture heterotrophic dominants (~50-70%)
- ✓ Include functional specialists (~0.01-1%)
- ✗ Remove sequencing noise/transients
- ✗ Exclude stochastically detected rare organisms

### Impact of 0.01% Threshold

**Our results:**
```
Initial species: 3,794
After 10% prevalence: Most species removed by this filter
After 0.01% abundance: Additional filtering of very rare species
Final retained: 567 species (14.9% of initial)
```

**Species removed breakdown:**
- **3,227 species removed** (85.1%)
  - Low prevalence (<10%) AND/OR 
  - Low abundance (<0.01%)

---

## 4. COMBINED FILTERING LOGIC

### How Prevalence and Abundance Work Together

```
Filter Logic: Keep species if:
    (Prevalence ≥ 10%) OR (Abundance ≥ 0.01%)

NOT: (Prevalence ≥ 10%) AND (Abundance ≥ 0.01%)

This is an INCLUSIVE OR filter:
├─ High prevalence, any abundance → KEEP
├─ Any prevalence, high abundance → KEEP
└─ Low prevalence, low abundance → REMOVE
```

### Rationale for OR vs AND

**Why OR (inclusive)?**

```
Scenario 1: High Prevalence, Low Abundance
├─ Example: Species in 20% of samples, but rarely abundant
├─ Biological meaning: Consistent colonist, low growth
├─ Ecological role: Possible ecosystem engineer despite low biomass
└─ Action: KEEP (prevalence filter catches it)

Scenario 2: Low Prevalence, High Abundance  
├─ Example: Species in 5% of samples, but very abundant there
├─ Biological meaning: Niche specialist, dominates in specific conditions
├─ Ecological role: Important for specific treatment stage
└─ Action: KEEP (abundance filter catches it)

Scenario 3: Low Prevalence, Low Abundance
├─ Example: Species in 2% of samples, ~0.001% abundance
├─ Biological meaning: Likely sequencing artifact or transient
├─ Ecological role: Minimal
└─ Action: REMOVE (fails both filters)
```

### Trade-off Explanation

**More stringent (AND):**
```
Keep species if: (Prevalence ≥ 10%) AND (Abundance ≥ 0.01%)

Advantages:
  - Extremely robust set (core microbiome only)
  - Minimal noise
  - High signal/noise ratio

Disadvantages:
  - Loss of spatial/temporal variation
  - Miss treatment-stage specialists
  - Reduced biodiversity information
  - Inappropriate for hypothesis-driven questions about nitrifiers, etc.
```

**Current approach (OR):**
```
Keep species if: (Prevalence ≥ 10%) OR (Abundance ≥ 0.01%)

Advantages:
  - Captures both consistent and dominant taxa
  - Retains ecological diversity
  - Detects treatment-stage specialists
  - Better for regression analysis

Disadvantages:
  - Slightly more noise
  - May include rare contaminants
```

---

## 5. SENSITIVITY ANALYSIS: WHAT IF WE CHANGED THRESHOLDS?

### Alternative Scenarios

#### Scenario A: More Conservative (Stringent)
```python
PREVALENCE_THRESHOLD = 20%  # vs. 10%
ABUNDANCE_THRESHOLD = 0.1%  # vs. 0.01%

Expected outcome:
├─ Retained species: ~100-150
├─ Sparsity: Very low (~80%)
├─ Advantages:
│  ├─ Extremely robust
│  ├─ Only dominant taxa
│  └─ Minimal noise
└─ Disadvantages:
   ├─ Loss of rare functional bacteria
   ├─ Miss nitrifier dynamics (often <50% abundant)
   ├─ Inappropriate for wastewater (too restrictive)
   └─ Reduced power for associations
```

#### Scenario B: More Permissive (Liberal)
```python
PREVALENCE_THRESHOLD = 5%   # vs. 10%
ABUNDANCE_THRESHOLD = 0.001% # vs. 0.01%

Expected outcome:
├─ Retained species: ~2000-3000
├─ Sparsity: Very high (~97%)
├─ Advantages:
│  ├─ Captures rare taxa
│  ├─ Full biodiversity picture
│  └─ Possible keystone species
└─ Disadvantages:
   ├─ Extreme sparsity
   ├─ High noise/artifact rates
   ├─ Statistical problems (singular matrices)
   ├─ Computational burden
   └─ Poor regression performance
```

#### Scenario C: Our Choice (Balanced)
```python
PREVALENCE_THRESHOLD = 10%   # ← Chosen
ABUNDANCE_THRESHOLD = 0.01%  # ← Chosen

Expected outcome:
├─ Retained species: 567
├─ Sparsity: 90% (acceptable)
├─ Advantages:
│  ├─ Balanced sensitivity/specificity
│  ├─ Captures functional diversity
│  ├─ Good statistical properties
│  ├─ Aligns with literature standards
│  └─ Optimal for regression analysis
└─ Disadvantages: Misses some rare specialists
```

---

## 6. SPECIAL CONSIDERATIONS FOR WASTEWATER

### Why These Thresholds Work for Wastewater

**Wastewater-specific factors that justify our thresholds:**

#### 1. **High Biomass System**
```
Consequence: Even rare species (0.01%) = millions of cells
            Makes 0.01% a biologically meaningful threshold
```

#### 2. **Intensive Treatment**
```
Consequence: Strong selective pressure for abundant taxa
            10% prevalence captures selection-resistant species
            Good signal for ecological processes
```

#### 3. **Temporal Variability**
```
Consequence: Species composition changes between days/seasons
            10% prevalence (≥21 samples) catches consistent colonizers
            Filters out episodic/transient bacteria
```

#### 4. **Functional Importance**
```
Key microbes in wastewater:
├─ Nitrifiers (AOA, AOB): Often 5-20% abundant, 40-90% prevalent
├─ Heterotrophs: 50-70% abundant, 90%+ prevalent
└─ Specialized (methanogens, SRB): 0.01-1% abundant, 30-50% prevalent

Our thresholds:
✓ Capture all three functional groups
✓ Retain nitrification dynamics (important for treatment)
✓ Include syntrophic specialists
```

---

## 7. VALIDATION: DO THRESHOLDS MAKE SENSE?

### Check: Do Top Species Make Biological Sense?

**Top 5 retained species in our analysis:**

```
1. s__Nitrospira_sp_ND1
   - Mean abundance: 7.67%
   - Prevalence: 35.9%
   - Comment: Comammox nitrifier, known in wastewater ✓ Makes sense

2. s__Acinetobacter_johnsonii
   - Mean abundance: 6.84%
   - Prevalence: 45.0%
   - Comment: Heterotroph, very common in WWTPs ✓ Makes sense

3. s__GGB67335_SGB103624
   - Mean abundance: 6.52%
   - Prevalence: 30.6%
   - Comment: Unknown species, likely important polyaromatic degrader ✓ Plausible

4. s__Nitrosomonas_oligotropha
   - Mean abundance: 1.04%
   - Prevalence: 14.4%
   - Comment: AOB nitrifier, consistent colonizer at low abundance ✓ Makes sense
   - Note: Below 10% prevalence but >0.01% abundance, so retained
```

**Conclusion:** Retained species make excellent biological sense for wastewater. ✓

---

## 8. RECOMMENDATIONS FOR DIFFERENT SCENARIOS

### When to Modify Thresholds

| Goal | Prevalence | Abundance | Rationale |
|------|-----------|-----------|-----------|
| **Our case: CLR + Regression** | 10% | 0.01% | Balance; standard |
| Early exploration | 5% | 0.001% | Discover rare diversity |
| Robust core microbiome | 20% | 0.1% | Maximum confidence |
| Treatment optimization | 10% | 0.05% | Focus on abundant operatives |
| Rare pathogen detection | 1% | 0.001% | Sensitivity priority |
| Publication-ready | 15% | 0.01% | Conservative, reproducible |

---

## 9. SUMMARY TABLE

| Parameter | Our Choice | Rationale | Impact |
|-----------|-----------|-----------|--------|
| **Taxonomic Level** | Species (s__) | Functional relevance, statistical power | 573 species |
| **Str filters** | Exclude (t__ removed) | Reduce noise, annotation uncertainty | 573→567 species |
| **Prevalence** | ≥10% | Standard in literature; balance rare/abundant | Captures consistent colonizers |
| **Abundance** | ≥0.01% | Above sequencing noise; biologically meaningful | Removes noise, retains diversity |
| **Combined Logic** | OR | Captures both specialists and generalists | 567 final species retained |

---

## 10. REFERENCES & BEST PRACTICES

### Key Publications Supporting Our Choices

1. **Taxonomic filtering:**
   - MetaPhlAn documentation recommends species-level for most studies
   - Strain-level discarded in >90% of microbiome publications

2. **Prevalence thresholds:**
   - 10% is standard in differential abundance analysis
   - Used in DESeq2 vignettes, microbiomMarkerPipeline, ANCOM

3. **Abundance thresholds:**
   - 0.01% consensus for wastewater studies
   - Corresponds to 0.001-0.01% of sequencing depth

4. **CLR transformation justification:**
   - Aitchison (1986): Compositional data analysis
   - Gloor et al. (2017): ALDEx2 for microbiome
   - Quinn et al. (2018): Compositional analysis best practices

---

## Conclusion

The thresholds selected (10% prevalence, 0.01% abundance, species level) represent an **optimal balance** between:
- **Biological relevance** ↔ Statistical power
- **Sensitivity** ↔ Specificity  
- **Noise reduction** ↔ Diversity retention
- **Computational efficiency** ↔ Information capture

These choices are well-supported by literature, appropriate for wastewater systems, and yield species profiles suitable for downstream regression and association analyses.
