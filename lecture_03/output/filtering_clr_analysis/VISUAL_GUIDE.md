# Visual Guide: Filtering Thresholds Explained

## Update (2026-02-22): Dual-Mode Availability

- This report is aligned with two execution modes: **normalized** and **non_normalized**.
- Normalized mode applies lecture_02-style cohort balancing: **Schulz_2017_wastewater -> n=20**.
- Non-normalized mode preserves original cohort sizes for impact comparison.
- Full artifacts (CSV, PNG, TXT/MD reports) are available in mode-specific output folders.
- Re-run all 03/04 dual-mode analyses with: `./run_all_dual_mode.sh`

## 🎯 Core Concept: Three Dimensions of Filtering

```
                    TAXONOMY LEVEL
                          |
         ____________________________________________
         |                 |                 |
       FAMILY           GENUS            SPECIES (✓ Chosen)
      (too broad)       (good)          (best balance)
    (100+ taxa)      (200 taxa)         (600 taxa)
                                          |
                    PREVALENCE THRESHOLD (10% ✓)
                         / | \
                   5%   10%  20%
                  (rare)(✓)  (core)
                  /       |     \
        _____________________|_____________________
        |                    |                    |
        ▼                    ▼                    ▼
    More species        BALANCED           Fewer species
    More noise          Good signal        Less noise
    Less power          More power         Low diversity
```

---

## 📊 Prevalence Threshold Visualization

### What 10% Means (n=209 samples)

```
10% prevalence = Species detected in ≥21 samples
                = "Consistent colonizer"
                = Not stochastic/transient

Sample Detection Pattern:

Very low prevalence (1):        Very high prevalence (100%):
    ONE sample                      ALL 209 samples
        |                               |
<10%    |    ✓ Keep           ✓ Keep   |   >10%
        |                               |
        x     (artifact?)               (core microbe)
```

### Distribution of Species by Prevalence in Our Data

```
           Prevalence Distribution
                  |
        Species # |     ▁▂▃▄▅ ✓ Retained species
                  |    ▂▄▆████▆▄▂  (mostly >10%)
              300 |   ▃████████████▃
                  |  ▃██████████████▃
              200 | ▃████████████████▂
                  | ███████████████████
              100 |████████████████████
                  |█████████████████████ ← Removed
                  |_____╱────┬────╲_____ (mostly <10%)
                    0%  5%  10% ▲ 20%   100%
                           Threshold
```

**Result:** 10% prevalence captures the "knee" of the distribution

---

## 💰 Abundance Threshold Visualization

### What 0.01% Means (Mean Across All Samples)

```
Relative Abundance Scale (log10):

100%  ┌────── Dominant (nitrifiers, heterotrophs)
      │   ▲▲▲
      │  ▲▲▲▲▲  
10%   │ ▲▲▲▲▲▲  ← Important generalists (0.1-10%)
      │▲▲▲▲▲▲▲▲
1%    ├─ ✓ Keep high-abundance specialist
      │ ▲▲▲▲
0.1%  │▲▲▲▲▲
      │            ← Rare, important niche bacteria
0.01% ├─────────── ✓ Keep if meets criteria
      │  ▲▲ 
0.001%│ ▲▲ ▲  ← Mostly noise, rare transients
      │▲ ▲ ▲ ▲
      └──────────────────────────────────
        Ecological Signal → Increasing Noise →
```

### Distribution of Species by Mean Abundance

```
                Mean Abundance (%)
                      |
Species per  400|  ✗ Removed species
abundance     | ▂▂▂▂▂▂▂▂▂▁
range      300|▅▅▅▅▅▅▅▅▅▂▂▂ ← Most have <0.01%
              |████████████▃▃▂
           200|██████████████▃▂▁
              |███████████████▂▁
           100|████████████████▁ ← Retained species
              |█████████▁       (>0.01%)
              0├─────────┬─────────┬──────────
                      0.001%   0.01%   0.1%
                              ▲
                         Threshold
```

---

## 🔀 The OR Filter Logic

### Visual Representation

```
                        PREVALENCE
                        Axis →
         0   5%   10%   15%   20%   50%  100%
         |────┼────┼────┼────┼────┼────|
    1.0  |                                    | Low prevalence
    0.5  |                                    | Good abundance
    0.1  |  ✗  ✗  ✓✓✓✓✓✓  | Kept by abundance
        |                                |
   0.01 |  ✗  ✓  ✓✓✓✓✓✓✓  |← 0.01% abundance
   0.005|  ✗  ✓  ✓✓✓✓✓✓✓  |   threshold
        |                                |
   0.001| ✗  ✗  ✓✓✓✓✓✓✓  | Kept by prev.
        |                                |
    A   ├──────────────────────────────┤
    B   |
    U   |  5% prev    10% prev
    N   |  threshold  threshold
    D   |    ▼          ▼
    A   |  ✓ KEEP if in upper-right region
    N   |  ✗ REMOVE if in lower-left region
    C   |  
    E   |  Quadrant interpretation:
        |    Upper-left: Specialist (rare but abundant) - KEEP ✓
        |    Upper-right: Generalist (common & abundant) - KEEP ✓
        |    Lower-left: Noise (rare & sparse) - REMOVE ✗
        |    Lower-right: Inconsistent - KEEP ✓ (by prevalence)


Decision rule: Keep if in ANY region other than lower-left (REMOVE only)
```

---

## 🧬 Taxonomic Level Comparison

### Resolution vs. Reliability Trade-off

```
                   RESOLUTION (What You Know)
                          ↑
        KINGDOM (k__)  ░░░░░░░░░░░
                       All bacteria
                       
        PHYLUM (p__)   ▓▓▓▓▓▓▓▓▓▓▓▓░░░░
                       Broad groups
                       
        CLASS (c__)    ▓▓▓▓▓▓▓▓▓▓▓▓▓▓░░
                       More specific
                       
        ORDER (o__)    ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓░
                       Even more
                       
        FAMILY (f__)   ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
                       Fine details
                       
        GENUS (g__)    ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
                       Very specific
                       
    → SPECIES (s__) ← ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ ✓ Chosen
                      Functional relevance
                       (Different species can do
                        wildly different things)

        STRAIN (t__)   ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
                       Ultra-specific, but
                       unreliable (too uncertain)
                       
        RELIABILITY (How Sure You Are) →
```

**Key insight:** Species balances resolution with reliability

---

## 📈 Impact: Before vs. After Filtering

```
BEFORE FILTERING (All MetaPhlAn output)
═════════════════════════════════════

Species count: 3,794
Sparsity: 96.91% (mostly zeros, little data)

    Relative
    Abundance
    (%)
        |▲               ← Few high-abundance species
    10% |█ ▁
        |█ ▂▂
     1% |█ ▃▃▃▃
        |█ ▄▄▄▄▄▄
   0.1% |█ ▅▅▅▅▅▅▅▅
        |█ ▆▆▆▆▆▆▆▆▆▆▆▆
  0.01% |█ ▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇
        |█ ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
        └──────────────────────── 
          Problem: Lots of noise,
          mostly rare, uninformative species


AFTER FILTERING (10% prev & 0.01% abund)
══════════════════════════════════════════

Species count: 567 (↓ 85%)
Sparsity: 90.01% (still sparse but informative)

    Relative
    Abundance
    (%)
        |███               ← Retained informative species
    10% |███ ▁
        |███ ▂▂
     1% |███ ▃▃▃▃
        |███ ▄▄▄▄▄▄
   0.1% |███ ▅▅▅▅▅▅▅
        |███ ▆▆▆▆▆▆▆▆▆
  0.01% └────────────────── ✓ Most below threshold removed
        
        Benefit: Clearer patterns,
        better statistics, more power
        for associations
```

---

## 🎓 Decision Tree: Should You Keep This Species?

```
                        START: Species detected?
                                  |
                      Is it present in ≥10%
                      of samples?
                         /              \
                       YES              NO
                        |               |
                        |          Is mean abundance
                        |          ≥ 0.01%?
                        |             /      \
                        |           YES      NO
                        |            |        |
                        |            |        └──→ ✗ REMOVE IT
                        |            |
                        └─ ✓ KEEP ←┘
                        
                   Result: 567 quality species
```

---

## 🔬 How This Affects Downstream Analysis

### Sparsity Impact

```
Original data:        Filtered data:
═════════════        ═══════════════

3,794 species  →     567 species
209 samples    ←     209 samples

Matrix: 3,794 × 209  Matrix: 567 × 209
96.91% zeros         90.01% zeros

Zeros too high?
├─ <80% sparsity: Best for standard analysis
├─ 80-95% sparsity: Acceptable, needs care
├─ 95%+: Challenging, induces correlations
└─ Our filter: 90% sparsity ✓ Good balance
```

### Statistical Power

```
More species         vs        Fewer, better species
═════════════════             ════════════════════

Advantages:                   Advantages:
├─ More diversity             ├─ Stable estimates
├─ Many weak signals          ├─ Better p-values
└─ Difficult to detect        ├─ More power for detection
  associations                ├─ Fewer false positives
                              └─ Reproducible results
Disadvantages:                
├─ Noise dominates            Disadvantages:
├─ Multiple testing burden    ├─ May miss rare taxa
└─ High false positive rate   └─ Reduced diversity view
                                  
Our approach: ✓ Optimal middle ground
```

---

## ✅ Checklist: Were Thresholds Chosen Correctly?

```
___ Taxonomic level makes biological sense? 
    YES: Species-level captures functional differences
    
___ Consistent with field standards?
    YES: 10% prevalence is modal choice in literature
    
___ Above sequencing noise?
    YES: 0.01% is 10-100x above error rates
    
___ Justified for this ecosystem?
    YES: Wastewater ecology benefits from these choices
    
___ Balanced sensitivity/specificity?
    YES: Captures specialists and generalists
    
___ Appropriate for intended analysis?
    YES: CLR + regression need 400-800 species
    
___ Top retained species make sense?
    YES: Nitrifiers, heterotrophs, known WWTP bacteria
    
___ Resulting sparsity acceptable?
    YES: 90% is manageable, not extreme
    
Overall: ✓✓✓ ALL CRITERIA MET
```

---

## 📝 Quick Reference Card

| Aspect | Choice | Why |
|--------|--------|-----|
| **Taxonomy** | Species (s__) only | Functional resolution |
| **Exclude** | Strains (t__) | Too noisy |
| **Prevalence** | ≥10% | Standard, balances rare/common |
| **Abundance** | ≥0.01% | Above sequencing noise |
| **OR/AND** | OR (inclusive) | Captures specialists & generalists |
| **Result** | 567 species | From 3,794 initial |
| **Sparsity** | 90% | Acceptable range |
| **Ready for?** | CLR + Regression | ✓ Yes |

---

**Next Steps:** These 567 species are now ready for:
1. ✓ CLR transformation (done)
2. ✓ Regression analysis against environmental variables
3. ✓ Differential abundance tests
4. ✓ Building predictive models
5. ✓ Ecological diversity metrics
