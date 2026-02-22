#!/usr/bin/env python3
"""
SENSITIVITY ANALYSIS: META-ANALYSIS WITH BALANCED COHORTS
==========================================================
Re-run meta-analysis with downsampled cohorts to test robustness.

Strategy:
- Use balanced data from lecture_03 (n=12 per cohort, 48 total)
- Re-run identical meta-analysis pipeline
- Compare results with original meta-analysis (209 samples)
- Assess concordance and effect size correlation
"""

import pandas as pd
import numpy as np
import os
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Configuration
RANDOM_SEED = 42
FDR_ALPHA = 0.05
OUTPUT_DIR = Path('output/sensitivity_meta_balanced')
OUTPUT_DIR = Path(os.getenv('OUTPUT_DIR', str(OUTPUT_DIR)))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Input files - use balanced data from lecture_03 sensitivity analysis
CLR_FILE = os.getenv('CLR_FILE', '../lecture_03/output/sensitivity_balanced/abundance_clr_balanced.csv')
METADATA_FILE = os.getenv('METADATA_FILE', '../lecture_03/output/sensitivity_balanced/metadata_balanced.csv')
ORIGINAL_META_RESULTS = os.getenv('ORIGINAL_META_RESULTS', 'output/meta_analysis/meta_analysis_results_all.csv')

print("="*70)
print("SENSITIVITY ANALYSIS: META-ANALYSIS WITH BALANCED COHORTS")
print("="*70)
print(f"Strategy: Use balanced cohorts (n=12 each) from lecture_03")
print(f"Random seed: {RANDOM_SEED}")
print()

# ============================================================================
# STEP 1: LOAD BALANCED DATA
# ============================================================================
print("STEP 1: LOAD BALANCED DATA FROM LECTURE_03")
print("-" * 70)

# Load balanced CLR-transformed data
clr_df = pd.read_csv(CLR_FILE, index_col=0)
if clr_df.index.duplicated().any():
    n_dup = clr_df.index.duplicated().sum()
    print(f"  Removing {n_dup} duplicate samples from CLR table...")
    clr_df = clr_df[~clr_df.index.duplicated(keep='first')]
print(f"  ✓ Loaded balanced CLR data: {clr_df.shape[0]} samples × {clr_df.shape[1]} species")

# Load balanced metadata
metadata = pd.read_csv(METADATA_FILE, index_col=0)
if metadata.index.duplicated().any():
    n_dup = metadata.index.duplicated().sum()
    print(f"  Removing {n_dup} duplicate samples from metadata...")
    metadata = metadata[~metadata.index.duplicated(keep='first')]
print(f"  ✓ Loaded balanced metadata: {len(metadata)} samples")

# Verify alignment
assert all(clr_df.index == metadata.index), "Sample mismatch!"
print(f"  ✓ Data aligned: {len(clr_df)} samples")

# Show cohort distribution
print("\n  Balanced cohort composition:")
cohort_counts = metadata['study_code'].value_counts().sort_index()
for cohort, count in cohort_counts.items():
    print(f"    • {cohort}: {count} samples")

n_cohorts = len(cohort_counts)
print(f"  Total cohorts: {n_cohorts}")

# ============================================================================
# STEP 2: BATCH EFFECT ASSESSMENT (BALANCED)
# ============================================================================
print("\nSTEP 2: BATCH EFFECT ASSESSMENT (BALANCED COHORTS)")
print("-" * 70)

# PCA - adjust n_components for small sample size
n_samples = len(clr_df)
n_components_pca = min(20, n_samples - 1)  # Can't exceed n_samples-1
scaler = StandardScaler()
clr_scaled = scaler.fit_transform(clr_df)
pca = PCA(n_components=n_components_pca)
pca_coords = pca.fit_transform(clr_scaled)

# Variance explained
var_exp = pca.explained_variance_ratio_ * 100
print(f"  PC1 explains {var_exp[0]:.1f}% variance")
print(f"  PC2 explains {var_exp[1]:.1f}% variance")
print(f"  PC1+PC2 explain {var_exp[0] + var_exp[1]:.1f}% variance")

# Test batch effects on first 5 PCs
print("\n  Batch effect tests (cohort vs PC):")
batch_results = []
n_test = min(5, pca_coords.shape[1])
for i in range(n_test):
    pc_data = pd.DataFrame({
        'PC': pca_coords[:, i],
        'cohort': metadata['study_code']
    })
    
    # ANOVA
    groups = [pc_data[pc_data['cohort'] == c]['PC'].values 
              for c in cohort_counts.index]
    f_stat, p_val = stats.f_oneway(*groups)
    
    sig = "YES" if p_val < 0.05 else "NO"
    print(f"    PC{i+1} (var={var_exp[i]:.1f}%): F={f_stat:.2f}, p={p_val:.2e} [{sig}]")
    
    batch_results.append({
        'PC': f'PC{i+1}',
        'variance_explained': var_exp[i],
        'F_statistic': f_stat,
        'p_value': p_val,
        'significant': p_val < 0.05
    })

batch_df = pd.DataFrame(batch_results)
batch_df.to_csv(OUTPUT_DIR / 'batch_effect_tests_balanced.csv', index=False)

# Save PCA coordinates
n_save = min(10, pca_coords.shape[1])  # Save up to 10 PCs
pca_df = pd.DataFrame(
    pca_coords[:, :n_save],
    index=clr_df.index,
    columns=[f'PC{i+1}' for i in range(n_save)]
)
pca_df['cohort'] = metadata['study_code']
pca_df.to_csv(OUTPUT_DIR / 'pca_coordinates_balanced.csv')

# ============================================================================
# STEP 3: META-ANALYSIS (BALANCED COHORTS)
# ============================================================================
print("\nSTEP 3: META-ANALYSIS (BALANCED COHORTS)")
print("-" * 70)

meta_results = []

for i, species in enumerate(clr_df.columns, 1):
    if i % 100 == 0:
        print(f"  Progress: {i}/{len(clr_df.columns)} species...")
    
    # Get abundance for this species
    abundance = clr_df[species]
    
    # Calculate overall mean
    overall_mean = abundance.mean()
    
    # Calculate cohort-specific effects
    cohort_effects = []
    cohort_ses = []
    cohort_ns = []
    
    for cohort in cohort_counts.index:
        cohort_mask = metadata['study_code'] == cohort
        cohort_abundance = abundance[cohort_mask]
        
        n = len(cohort_abundance)
        cohort_mean = cohort_abundance.mean()
        cohort_std = cohort_abundance.std()
        
        # Effect size: deviation from overall mean
        effect = cohort_mean - overall_mean
        
        # Standard error
        se = cohort_std / np.sqrt(n)
        
        cohort_effects.append(effect)
        cohort_ses.append(se)
        cohort_ns.append(n)
    
    cohort_effects = np.array(cohort_effects)
    cohort_ses = np.array(cohort_ses)
    cohort_ns = np.array(cohort_ns)
    
    # Inverse-variance weights
    weights = 1 / (cohort_ses ** 2)
    weights[np.isinf(weights)] = 0  # Handle zero SE
    
    # Pooled effect (fixed-effects)
    if weights.sum() > 0:
        pooled_effect = np.sum(weights * cohort_effects) / np.sum(weights)
        pooled_se = np.sqrt(1 / np.sum(weights))
    else:
        pooled_effect = 0
        pooled_se = np.inf
    
    # Z-test
    if pooled_se > 0 and not np.isinf(pooled_se):
        z_score = pooled_effect / pooled_se
        p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
    else:
        z_score = 0
        p_value = 1.0
    
    # Heterogeneity: Q statistic and I²
    Q = np.sum(weights * (cohort_effects - pooled_effect) ** 2)
    df = len(cohort_effects) - 1
    Q_pvalue = 1 - stats.chi2.cdf(Q, df) if df > 0 else 1.0
    I2 = max(0, ((Q - df) / Q) * 100) if Q > 0 else 0
    
    # Effect range
    effect_range = cohort_effects.max() - cohort_effects.min()
    
    meta_results.append({
        'species': species,
        'n_cohorts': n_cohorts,
        'pooled_effect': pooled_effect,
        'pooled_se': pooled_se,
        'z_score': z_score,
        'p_value': p_value,
        'Q_statistic': Q,
        'Q_pvalue': Q_pvalue,
        'I2_heterogeneity': I2,
        'min_effect': cohort_effects.min(),
        'max_effect': cohort_effects.max(),
        'effect_range': effect_range
    })

print(f"  Completed: {len(meta_results)} species tested")

# Create results dataframe
results_balanced = pd.DataFrame(meta_results)

# FDR correction
results_balanced['p_adjusted'] = stats.false_discovery_control(
    results_balanced['p_value'],
    method='bh'
)
results_balanced['significant'] = results_balanced['p_adjusted'] < FDR_ALPHA

# Sort by p-value
results_balanced = results_balanced.sort_values('p_value')

# Save balanced results
results_balanced.to_csv(OUTPUT_DIR / 'meta_analysis_balanced_all.csv', index=False)
results_balanced[results_balanced['significant']].to_csv(
    OUTPUT_DIR / 'meta_analysis_balanced_significant.csv',
    index=False
)

n_sig_balanced = results_balanced['significant'].sum()
pct_sig_balanced = (n_sig_balanced / len(results_balanced)) * 100

print(f"\n  Results (balanced n=48):")
print(f"    • Species tested: {len(results_balanced)}")
print(f"    • Significant (FDR<0.05): {n_sig_balanced}/{len(results_balanced)} ({pct_sig_balanced:.1f}%)")
print(f"    • Median I² heterogeneity: {results_balanced['I2_heterogeneity'].median():.1f}%")

# ============================================================================
# STEP 4: COMPARE WITH ORIGINAL META-ANALYSIS
# ============================================================================
print("\nSTEP 4: COMPARE WITH ORIGINAL META-ANALYSIS (n=209)")
print("-" * 70)

# Load original meta-analysis results
results_original = pd.read_csv(ORIGINAL_META_RESULTS)
if results_original.duplicated(subset=['species']).any():
    n_dup = results_original.duplicated(subset=['species']).sum()
    print(f"  Removing {n_dup} duplicate species rows from original meta-analysis results...")
    results_original = results_original.sort_values('p_value').drop_duplicates(subset=['species'], keep='first')
print(f"  ✓ Loaded original meta-analysis: {len(results_original)} species")

n_sig_original = results_original['significant'].sum()
pct_sig_original = (n_sig_original / len(results_original)) * 100

print(f"  Results (original n=209):")
print(f"    • Species tested: {len(results_original)}")
print(f"    • Significant (FDR<0.05): {n_sig_original}/{len(results_original)} ({pct_sig_original:.1f}%)")
print(f"    • Median I² heterogeneity: {results_original['I2_heterogeneity'].median():.1f}%")

# Merge results
comparison = results_original.merge(
    results_balanced,
    on='species',
    suffixes=('_original', '_balanced')
)

# Concordance analysis
both_sig = (comparison['significant_original'] & comparison['significant_balanced']).sum()
original_only = (comparison['significant_original'] & ~comparison['significant_balanced']).sum()
balanced_only = (~comparison['significant_original'] & comparison['significant_balanced']).sum()
neither_sig = (~comparison['significant_original'] & ~comparison['significant_balanced']).sum()

print(f"\n  Concordance:")
print(f"    • Both significant: {both_sig}/{len(comparison)} ({both_sig/len(comparison)*100:.1f}%)")
print(f"    • Original only: {original_only} ({original_only/len(comparison)*100:.1f}%)")
print(f"    • Balanced only: {balanced_only} ({balanced_only/len(comparison)*100:.1f}%)")
print(f"    • Neither: {neither_sig} ({neither_sig/len(comparison)*100:.1f}%)")

# Effect size correlation
eff_corr = comparison['pooled_effect_original'].corr(comparison['pooled_effect_balanced'])
print(f"\n  Pooled effect correlation: r = {eff_corr:.3f}")

# P-value correlation (log scale)
comparison['log_p_original'] = -np.log10(comparison['p_value_original'] + 1e-320)
comparison['log_p_balanced'] = -np.log10(comparison['p_value_balanced'] + 1e-320)
p_corr = comparison['log_p_original'].corr(comparison['log_p_balanced'])
print(f"  P-value correlation: r = {p_corr:.3f}")

# I² correlation
i2_corr = comparison['I2_heterogeneity_original'].corr(comparison['I2_heterogeneity_balanced'])
print(f"  I² heterogeneity correlation: r = {i2_corr:.3f}")

# Save comparison
comparison.to_csv(OUTPUT_DIR / 'comparison_original_vs_balanced_meta.csv', index=False)
print(f"\n  ✓ Saved comparison to comparison_original_vs_balanced_meta.csv")

# ============================================================================
# STEP 5: VISUALIZATIONS
# ============================================================================
print("\nSTEP 5: CREATE COMPARISON VISUALIZATIONS")
print("-" * 70)

# Set style
sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 300

# ==================
# FIGURE 1: Meta-analysis comparison overview
# ==================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: Pooled effect correlation
ax = axes[0, 0]
colors = ['red' if (s1 and s2) else 'gray'
          for s1, s2 in zip(comparison['significant_original'],
                           comparison['significant_balanced'])]
ax.scatter(comparison['pooled_effect_original'],
          comparison['pooled_effect_balanced'],
          alpha=0.5, s=20, c=colors, edgecolors='none')
max_val = max(abs(comparison['pooled_effect_original']).max(),
              abs(comparison['pooled_effect_balanced']).max())
ax.plot([-max_val, max_val], [-max_val, max_val], 'k--', lw=1, alpha=0.5, label='y=x')
ax.set_xlabel('Pooled Effect - Original (n=209)', fontsize=10)
ax.set_ylabel('Pooled Effect - Balanced (n=48)', fontsize=10)
ax.set_title(f'A) Pooled Effect Correlation (r={eff_corr:.3f})', fontsize=11, fontweight='bold')
ax.legend(['y=x', 'Both sig', 'Not both'], fontsize=8)
ax.grid(True, alpha=0.3)

# Panel B: P-value correlation
ax = axes[0, 1]
ax.scatter(comparison['log_p_original'],
          comparison['log_p_balanced'],
          alpha=0.5, s=20, c=colors, edgecolors='none')
max_val = max(comparison['log_p_original'].max(), comparison['log_p_balanced'].max())
ax.plot([0, max_val], [0, max_val], 'k--', lw=1, alpha=0.5)
ax.set_xlabel('-log10(p) - Original (n=209)', fontsize=10)
ax.set_ylabel('-log10(p) - Balanced (n=48)', fontsize=10)
ax.set_title(f'B) P-value Correlation (r={p_corr:.3f})', fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3)

# Panel C: Concordance bar chart
ax = axes[1, 0]
categories = ['Both\nsignificant', 'Original\nonly', 'Balanced\nonly', 'Neither']
counts = [both_sig, original_only, balanced_only, neither_sig]
colors_bar = ['#2ecc71', '#e74c3c', '#3498db', '#95a5a6']
bars = ax.bar(categories, counts, color=colors_bar, alpha=0.7, edgecolor='black', linewidth=1)
ax.set_ylabel('Number of Species', fontsize=10)
ax.set_title('C) Significance Concordance', fontsize=11, fontweight='bold')
ax.grid(axis='y', alpha=0.3)
# Add counts on bars
for bar, count in zip(bars, counts):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{count}\n({count/len(comparison)*100:.1f}%)',
            ha='center', va='bottom', fontsize=9)

# Panel D: I² comparison
ax = axes[1, 1]
ax.scatter(comparison['I2_heterogeneity_original'],
          comparison['I2_heterogeneity_balanced'],
          alpha=0.4, s=20, c=colors, edgecolors='none')
ax.plot([0, 100], [0, 100], 'k--', lw=1, alpha=0.5)
ax.set_xlabel('I² Heterogeneity - Original', fontsize=10)
ax.set_ylabel('I² Heterogeneity - Balanced', fontsize=10)
ax.set_title(f'D) Heterogeneity Correlation (r={i2_corr:.3f})', fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.set_xlim([0, 100])
ax.set_ylim([0, 100])

plt.tight_layout()
plt.savefig(OUTPUT_DIR / '01_meta_sensitivity_comparison.png', dpi=300, bbox_inches='tight')
print("  ✓ Saved 01_meta_sensitivity_comparison.png")
plt.close()

# ==================
# FIGURE 2: Batch effects comparison
# ==================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: PCA (balanced)
ax = axes[0]
cohort_colors = {'Chu_2017_sludge': '#e74c3c',
                'Lekunberri_2018_river_wastewater': '#3498db',
                'Rowe_2017_hospital_wastewater': '#2ecc71',
                'Schulz_2017_wastewater': '#9b59b6'}

for cohort in cohort_counts.index:
    mask = pca_df['cohort'] == cohort
    ax.scatter(pca_df.loc[mask, 'PC1'],
              pca_df.loc[mask, 'PC2'],
              label=cohort.replace('_', ' '),
              alpha=0.7, s=60,
              color=cohort_colors.get(cohort, 'gray'))

ax.set_xlabel(f'PC1 ({var_exp[0]:.1f}% variance)', fontsize=10)
ax.set_ylabel(f'PC2 ({var_exp[1]:.1f}% variance)', fontsize=10)
ax.set_title('A) PCA - Balanced Cohorts (n=12 each)', fontsize=11, fontweight='bold')
ax.legend(fontsize=7, loc='best')
ax.grid(True, alpha=0.3)

# Panel B: Significance rates comparison
ax = axes[1]
datasets = ['Original\n(n=209)', 'Balanced\n(n=48)']
sig_rates = [pct_sig_original, pct_sig_balanced]
bars = ax.bar(datasets, sig_rates, color=['#9b59b6', '#e67e22'], alpha=0.7, edgecolor='black', linewidth=1)
ax.set_ylabel('% Significant Species (FDR<0.05)', fontsize=10)
ax.set_title('B) Meta-Analysis Significance Rate', fontsize=11, fontweight='bold')
ax.set_ylim([0, 100])
ax.grid(axis='y', alpha=0.3)
# Add percentages
for bar, rate in zip(bars, sig_rates):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{rate:.1f}%',
            ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / '02_batch_and_significance.png', dpi=300, bbox_inches='tight')
print("  ✓ Saved 02_batch_and_significance.png")
plt.close()

# ============================================================================
# STEP 6: SUMMARY REPORT
# ============================================================================
print("\nSTEP 6: GENERATE SUMMARY REPORT")
print("-" * 70)

report = f"""SENSITIVITY ANALYSIS REPORT: META-ANALYSIS WITH BALANCED COHORTS
{'='*70}

OBJECTIVE
---------
Test robustness of meta-analysis results to sample size imbalance by
using balanced cohorts (n=12 each) from lecture_03 sensitivity analysis.

SAMPLE SIZE COMPARISON
---------------------
Original meta-analysis (lecture_04):
  • Total samples: 209 (from lecture_03 unbalanced)
  • Schulz: 128 samples (61%)
  • Chu: 49 samples (23%)
  • Rowe: 20 samples (10%)
  • Lekunberri: 12 samples (6%)
  • Imbalance ratio: 10.7x

Balanced meta-analysis:
  • Total samples: 48
  • All cohorts: 12 samples each (25%)
  • Imbalance ratio: 1.0x (perfectly balanced)

BATCH EFFECTS (BALANCED COHORTS)
--------------------------------
PCA variance explained:
  • PC1: {var_exp[0]:.1f}%
  • PC2: {var_exp[1]:.1f}%
  • PC1+PC2: {var_exp[0] + var_exp[1]:.1f}%

ANOVA tests (cohort effect on PCs):
"""

for i, row in batch_df.iterrows():
    sig = "YES" if row['significant'] else "NO"
    report += f"  • {row['PC']}: F={row['F_statistic']:.2f}, p={row['p_value']:.2e} [{sig}]\n"

report += f"""

META-ANALYSIS RESULTS
---------------------
Original (n=209):
  • Species tested: {len(results_original)}
  • Significant (FDR<0.05): {n_sig_original}/{len(results_original)} ({pct_sig_original:.1f}%)
  • Median I² heterogeneity: {results_original['I2_heterogeneity'].median():.1f}%

Balanced (n=48):
  • Species tested: {len(results_balanced)}
  • Significant (FDR<0.05): {n_sig_balanced}/{len(results_balanced)} ({pct_sig_balanced:.1f}%)
  • Median I² heterogeneity: {results_balanced['I2_heterogeneity'].median():.1f}%

CONCORDANCE ANALYSIS
--------------------
Significance concordance:
  • Both significant: {both_sig}/{len(comparison)} ({both_sig/len(comparison)*100:.1f}%)
  • Original only: {original_only} ({original_only/len(comparison)*100:.1f}%)
  • Balanced only: {balanced_only} ({balanced_only/len(comparison)*100:.1f}%)
  • Neither: {neither_sig} ({neither_sig/len(comparison)*100:.1f}%)

Correlation metrics:
  • Pooled effect correlation: r = {eff_corr:.3f}
  • P-value correlation: r = {p_corr:.3f}
  • I² heterogeneity correlation: r = {i2_corr:.3f}

TOP 10 SPECIES - ORIGINAL META-ANALYSIS RANKING
-----------------------------------------------
"""

for idx, row in results_original.head(10).iterrows():
    species_name = row['species']
    orig_p = row['p_value']
    orig_eff = row['pooled_effect']
    orig_i2 = row['I2_heterogeneity']
    
    # Find in balanced results
    bal_row = results_balanced[results_balanced['species'] == species_name]
    if len(bal_row) > 0:
        bal_p = bal_row['p_value'].values[0]
        bal_eff = bal_row['pooled_effect'].values[0]
        bal_i2 = bal_row['I2_heterogeneity'].values[0]
        bal_sig = "YES" if bal_row['significant'].values[0] else "NO"
        bal_rank = results_balanced[results_balanced['species'] == species_name].index[0] + 1
    else:
        bal_p = np.nan
        bal_eff = np.nan
        bal_i2 = np.nan
        bal_sig = "N/A"
        bal_rank = "N/A"
    
    report += f"\n{idx+1}. {species_name}\n"
    report += f"   Original: effect={orig_eff:.2f}, p={orig_p:.2e}, I²={orig_i2:.1f}%\n"
    report += f"   Balanced: effect={bal_eff:.2f}, p={bal_p:.2e}, I²={bal_i2:.1f}%, rank={bal_rank}, sig={bal_sig}\n"

report += f"""

KEY FINDINGS
------------
1. Significance rate: {pct_sig_original:.1f}% → {pct_sig_balanced:.1f}% with balanced samples
   • {'Decreased' if pct_sig_balanced < pct_sig_original else 'Increased'} due to {'reduced' if pct_sig_balanced < pct_sig_original else 'changed'} statistical power

2. Concordance: {both_sig/len(comparison)*100:.1f}% of species significant in BOTH analyses
   • {'High' if both_sig/len(comparison) > 0.7 else 'Moderate' if both_sig/len(comparison) > 0.5 else 'Low'} concordance indicates {'robust' if both_sig/len(comparison) > 0.7 else 'moderately robust' if both_sig/len(comparison) > 0.5 else 'sample size dependent'} findings

3. Pooled effect correlation: r = {eff_corr:.3f}
   • {'Strong' if abs(eff_corr) > 0.8 else 'Moderate' if abs(eff_corr) > 0.5 else 'Weak'} correlation
   • Effect sizes {'highly consistent' if abs(eff_corr) > 0.8 else 'moderately consistent' if abs(eff_corr) > 0.5 else 'vary'} between analyses

4. I² heterogeneity correlation: r = {i2_corr:.3f}
   • Heterogeneity patterns {'preserved' if i2_corr > 0.5 else 'change'} with balanced samples

INTERPRETATION
--------------
"""

if both_sig/len(comparison) > 0.7 and abs(eff_corr) > 0.7:
    report += "✓ Meta-analysis results are ROBUST to sample size imbalance.\n"
    report += "✓ Pooled effects are consistent regardless of balancing.\n"
    report += "✓ Cross-cohort patterns genuine, not driven by Schulz dominance.\n"
    report += "✓ Original meta-analysis (n=209) is valid and preferred (higher power).\n"
elif both_sig/len(comparison) > 0.5:
    report += "⚠ Meta-analysis shows MODERATE robustness to sample imbalance.\n"
    report += "⚠ Focus on species significant in BOTH analyses for highest confidence.\n"
    report += "⚠ Some original findings may be influenced by Schulz cohort size.\n"
else:
    report += "⚠ Meta-analysis results are SENSITIVE to sample size imbalance.\n"
    report += "⚠ Balanced results may better reflect true cross-cohort patterns.\n"
    report += "⚠ Original results may overweight Schulz cohort contributions.\n"

report += f"""

RECOMMENDATION
--------------
"""

if both_sig/len(comparison) > 0.7 and abs(eff_corr) > 0.7:
    report += "✓ Use original meta-analysis (n=209) as primary results.\n"
    report += f"✓ {both_sig} species significant in both = high-confidence meta-biomarkers.\n"
    report += "✓ Report sensitivity analysis as validation of robustness.\n"
elif both_sig/len(comparison) > 0.5:
    report += "⚠ Report both analyses with interpretation caveats.\n"
    report += f"⚠ Prioritize {both_sig} species significant in both analyses.\n"
    report += "⚠ Note sample size effects in limitations.\n"
else:
    report += "⚠ Consider balanced results for more generalizable conclusions.\n"
    report += "⚠ Original results may reflect Schulz-specific patterns.\n"

report += f"""

FILES GENERATED
---------------
• meta_analysis_balanced_all.csv - All species (balanced, n=48)
• meta_analysis_balanced_significant.csv - Significant species only
• comparison_original_vs_balanced_meta.csv - Side-by-side comparison
• batch_effect_tests_balanced.csv - Batch effect ANOVA results
• pca_coordinates_balanced.csv - PCA coordinates
• 01_meta_sensitivity_comparison.png - 4-panel comparison
• 02_batch_and_significance.png - PCA + significance rates
• sensitivity_meta_report.txt - This report

Analysis completed: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
"""

# Save report
report_path = OUTPUT_DIR / 'sensitivity_meta_report.txt'
with open(report_path, 'w') as f:
    f.write(report)

print(f"  ✓ Saved sensitivity_meta_report.txt")
print()
print("="*70)
print("SENSITIVITY ANALYSIS COMPLETE")
print("="*70)
print(f"\nKey Result: {both_sig/len(comparison)*100:.1f}% concordance between meta-analyses")
print(f"Effect correlation: r = {eff_corr:.3f}")
print(f"Interpretation: {'Meta-analysis is ROBUST' if both_sig/len(comparison) > 0.7 and abs(eff_corr) > 0.7 else 'Meta-analysis is MODERATELY ROBUST' if both_sig/len(comparison) > 0.5 else 'Meta-analysis is SENSITIVE to sample size'}")
print(f"\nAll outputs saved to: {OUTPUT_DIR}/")
print()
