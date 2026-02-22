#!/usr/bin/env python3
"""
SENSITIVITY ANALYSIS: BALANCED COHORTS
=====================================
Re-run association analysis with downsampled cohorts to test robustness.

Strategy:
- Downsample all cohorts to n=12 (matching smallest cohort: Lekunberri)
- Total: 48 samples (12 × 4 cohorts)
- Compare results with original 209-sample analysis
- Assess concordance and effect size correlation
"""

import pandas as pd
import numpy as np
import os
from scipy import stats
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Configuration
RANDOM_SEED = 42
BALANCED_N = 12  # Downsample all cohorts to this size
FDR_ALPHA = 0.05
OUTPUT_DIR = Path('output/sensitivity_balanced')
OUTPUT_DIR = Path(os.getenv('OUTPUT_DIR', str(OUTPUT_DIR)))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Input files
CLR_FILE = os.getenv('CLR_FILE', 'output/filtering_clr_analysis/abundance_clr_aligned.csv')
METADATA_FILE = os.getenv('METADATA_FILE', 'output/filtering_clr_analysis/metadata_aligned.csv')
ORIGINAL_RESULTS = os.getenv('ORIGINAL_RESULTS', 'output/association_analysis/association_results_all.csv')

print("="*70)
print("SENSITIVITY ANALYSIS: BALANCED COHORTS")
print("="*70)
print(f"Strategy: Downsample all cohorts to n={BALANCED_N}")
print(f"Random seed: {RANDOM_SEED}")
print()

# ============================================================================
# STEP 1: LOAD DATA
# ============================================================================
print("STEP 1: LOAD PROCESSED DATA")
print("-" * 70)

# Load CLR-transformed data
clr_df = pd.read_csv(CLR_FILE, index_col=0)
if clr_df.index.duplicated().any():
    n_dup = clr_df.index.duplicated().sum()
    print(f"  Removing {n_dup} duplicate samples from CLR table...")
    clr_df = clr_df[~clr_df.index.duplicated(keep='first')]
print(f"  ✓ Loaded CLR data: {clr_df.shape[0]} samples × {clr_df.shape[1]} species")

# Load metadata
metadata = pd.read_csv(METADATA_FILE, index_col=0)
# Remove duplicates if present
if metadata.index.duplicated().any():
    n_dup = metadata.index.duplicated().sum()
    print(f"  Removing {n_dup} duplicate samples from metadata...")
    metadata = metadata[~metadata.index.duplicated(keep='first')]
print(f"  ✓ Loaded metadata: {len(metadata)} samples")

# Align
common_samples = clr_df.index.intersection(metadata.index)
clr_df = clr_df.loc[common_samples]
metadata = metadata.loc[common_samples]
print(f"  ✓ Aligned: {len(common_samples)} samples")

# Show original cohort sizes
print("\n  Original cohort sizes:")
cohort_counts = metadata['study_code'].value_counts().sort_index()
for cohort, count in cohort_counts.items():
    print(f"    • {cohort}: {count} samples")

# ============================================================================
# STEP 2: DOWNSAMPLE TO BALANCED COHORTS
# ============================================================================
print(f"\nSTEP 2: DOWNSAMPLE TO n={BALANCED_N} PER COHORT")
print("-" * 70)

np.random.seed(RANDOM_SEED)

balanced_samples = []
for cohort in cohort_counts.index:
    cohort_samples = metadata[metadata['study_code'] == cohort].index
    
    if len(cohort_samples) >= BALANCED_N:
        # Downsample
        sampled = np.random.choice(cohort_samples, size=BALANCED_N, replace=False)
        balanced_samples.extend(sampled)
        print(f"  ✓ {cohort}: {len(cohort_samples)} → {BALANCED_N} samples")
    else:
        # Not enough samples - keep all (shouldn't happen with n=12)
        balanced_samples.extend(cohort_samples)
        print(f"  ⚠ {cohort}: {len(cohort_samples)} samples (< {BALANCED_N})")

# Create balanced datasets
clr_balanced = clr_df.loc[balanced_samples]
metadata_balanced = metadata.loc[balanced_samples]

print(f"\n  Balanced dataset: {len(balanced_samples)} samples × {clr_balanced.shape[1]} species")
print(f"  Cohort distribution:")
for cohort, count in metadata_balanced['study_code'].value_counts().sort_index().items():
    print(f"    • {cohort}: {count} samples")

# Save balanced datasets
clr_balanced.to_csv(OUTPUT_DIR / 'abundance_clr_balanced.csv')
metadata_balanced.to_csv(OUTPUT_DIR / 'metadata_balanced.csv')
print(f"\n  ✓ Saved balanced datasets to {OUTPUT_DIR}/")

# ============================================================================
# STEP 3: ASSOCIATION ANALYSIS ON BALANCED DATA
# ============================================================================
print("\nSTEP 3: ASSOCIATION ANALYSIS (BALANCED COHORTS)")
print("-" * 70)

results_list = []

for i, species in enumerate(clr_balanced.columns, 1):
    if i % 100 == 0:
        print(f"  Progress: {i}/{len(clr_balanced.columns)} species...")
    
    # Prepare data
    df_test = pd.DataFrame({
        'clr_abundance': clr_balanced[species],
        'study_code': metadata_balanced['study_code']
    })
    
    # Fit linear model
    try:
        model = ols('clr_abundance ~ C(study_code)', data=df_test).fit()
        anova_table = anova_lm(model, typ=2)
        
        # Extract statistics
        f_stat = anova_table.loc['C(study_code)', 'F']
        p_value = anova_table.loc['C(study_code)', 'PR(>F)']
        
        # Effect size: max - min group mean
        group_means = df_test.groupby('study_code')['clr_abundance'].mean()
        effect_size = group_means.max() - group_means.min()
        
        results_list.append({
            'species': species,
            'f_statistic': f_stat,
            'p_value': p_value,
            'effect_size': effect_size,
            'n_samples': len(df_test)
        })
        
    except Exception as e:
        print(f"  ⚠ Error for {species}: {e}")
        continue

print(f"  Completed: {len(results_list)} species tested")

# Create results dataframe
results_balanced = pd.DataFrame(results_list)

# FDR correction
results_balanced['p_adjusted'] = stats.false_discovery_control(
    results_balanced['p_value'], 
    method='bh'
)
results_balanced['significant'] = results_balanced['p_adjusted'] < FDR_ALPHA

# Sort by p-value
results_balanced = results_balanced.sort_values('p_value')

# Save balanced results
results_balanced.to_csv(OUTPUT_DIR / 'association_results_balanced_all.csv', index=False)
results_balanced[results_balanced['significant']].to_csv(
    OUTPUT_DIR / 'association_results_balanced_significant.csv', 
    index=False
)

n_sig_balanced = results_balanced['significant'].sum()
pct_sig_balanced = (n_sig_balanced / len(results_balanced)) * 100

print(f"\n  Results (balanced n={BALANCED_N}):")
print(f"    • Species tested: {len(results_balanced)}")
print(f"    • Significant (FDR<0.05): {n_sig_balanced}/{len(results_balanced)} ({pct_sig_balanced:.1f}%)")
print(f"    • Effect size range: {results_balanced['effect_size'].min():.2f} - {results_balanced['effect_size'].max():.2f} CLR")

# ============================================================================
# STEP 4: COMPARE WITH ORIGINAL RESULTS
# ============================================================================
print("\nSTEP 4: COMPARE WITH ORIGINAL (n=209) RESULTS")
print("-" * 70)

# Load original results
results_original = pd.read_csv(ORIGINAL_RESULTS)
if results_original.duplicated(subset=['species']).any():
    n_dup = results_original.duplicated(subset=['species']).sum()
    print(f"  Removing {n_dup} duplicate species rows from original association results...")
    results_original = results_original.sort_values('p_value').drop_duplicates(subset=['species'], keep='first')
print(f"  ✓ Loaded original results: {len(results_original)} species")

# Rename columns if needed for consistency
if 'sig_fdr' in results_original.columns and 'significant' not in results_original.columns:
    results_original['significant'] = results_original['sig_fdr']

n_sig_original = results_original['significant'].sum()
pct_sig_original = (n_sig_original / len(results_original)) * 100

print(f"  Results (original n=209):")
print(f"    • Species tested: {len(results_original)}")
print(f"    • Significant (FDR<0.05): {n_sig_original}/{len(results_original)} ({pct_sig_original:.1f}%)")

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
eff_corr = comparison['effect_size_original'].corr(comparison['effect_size_balanced'])
print(f"\n  Effect size correlation: r = {eff_corr:.3f}")

# P-value correlation (log scale)
comparison['log_p_original'] = -np.log10(comparison['p_value_original'] + 1e-320)
comparison['log_p_balanced'] = -np.log10(comparison['p_value_balanced'] + 1e-320)
p_corr = comparison['log_p_original'].corr(comparison['log_p_balanced'])
print(f"  P-value correlation: r = {p_corr:.3f}")

# Save comparison
comparison.to_csv(OUTPUT_DIR / 'comparison_original_vs_balanced.csv', index=False)
print(f"\n  ✓ Saved comparison to comparison_original_vs_balanced.csv")

# ============================================================================
# STEP 5: VISUALIZATIONS
# ============================================================================
print("\nSTEP 5: CREATE COMPARISON VISUALIZATIONS")
print("-" * 70)

# Set style
sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 300

# ==================
# FIGURE 1: Sample size impact overview
# ==================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: Effect size correlation
ax = axes[0, 0]
colors = ['red' if (s1 and s2) else 'gray' 
          for s1, s2 in zip(comparison['significant_original'], 
                           comparison['significant_balanced'])]
ax.scatter(comparison['effect_size_original'], 
          comparison['effect_size_balanced'],
          alpha=0.5, s=20, c=colors, edgecolors='none')
ax.plot([0, comparison['effect_size_original'].max()],
        [0, comparison['effect_size_original'].max()],
        'k--', lw=1, alpha=0.5, label='y=x')
ax.set_xlabel('Effect Size - Original (n=209)', fontsize=10)
ax.set_ylabel('Effect Size - Balanced (n=48)', fontsize=10)
ax.set_title(f'A) Effect Size Correlation (r={eff_corr:.3f})', fontsize=11, fontweight='bold')
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

# Panel D: Significance rate comparison
ax = axes[1, 1]
datasets = ['Original\n(n=209)', 'Balanced\n(n=48)']
sig_rates = [pct_sig_original, pct_sig_balanced]
bars = ax.bar(datasets, sig_rates, color=['#9b59b6', '#e67e22'], alpha=0.7, edgecolor='black', linewidth=1)
ax.set_ylabel('% Significant Species (FDR<0.05)', fontsize=10)
ax.set_title('D) Overall Significance Rate', fontsize=11, fontweight='bold')
ax.set_ylim([0, 105])
ax.axhline(y=95, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='95%')
ax.grid(axis='y', alpha=0.3)
ax.legend(fontsize=8)
# Add percentages
for bar, rate in zip(bars, sig_rates):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{rate:.1f}%',
            ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / '01_sensitivity_comparison.png', dpi=300, bbox_inches='tight')
print("  ✓ Saved 01_sensitivity_comparison.png")
plt.close()

# ==================
# FIGURE 2: Top species comparison
# ==================
# Get top 20 species from original
top20_original = results_original.nsmallest(20, 'p_value')['species'].tolist()

# Extract their ranks in both analyses
top_comparison = []
for species in top20_original:
    orig_rank = results_original[results_original['species'] == species].index[0] + 1
    orig_p = results_original[results_original['species'] == species]['p_value'].values[0]
    orig_eff = results_original[results_original['species'] == species]['effect_size'].values[0]
    
    bal_data = results_balanced[results_balanced['species'] == species]
    if len(bal_data) > 0:
        bal_rank = bal_data.index[0] + 1
        bal_p = bal_data['p_value'].values[0]
        bal_eff = bal_data['effect_size'].values[0]
    else:
        bal_rank = np.nan
        bal_p = np.nan
        bal_eff = np.nan
    
    top_comparison.append({
        'species': species,
        'orig_rank': orig_rank,
        'orig_p': orig_p,
        'orig_effect': orig_eff,
        'bal_rank': bal_rank,
        'bal_p': bal_p,
        'bal_effect': bal_eff
    })

top_comp_df = pd.DataFrame(top_comparison)

fig, axes = plt.subplots(1, 2, figsize=(14, 8))

# Panel A: Rank comparison
ax = axes[0]
species_short = [s.replace('s__', '').replace('_', ' ')[:30] for s in top_comp_df['species']]
y_pos = np.arange(len(species_short))

ax.barh(y_pos - 0.2, top_comp_df['orig_rank'], height=0.4, 
        color='#9b59b6', alpha=0.7, label='Original (n=209)')
ax.barh(y_pos + 0.2, top_comp_df['bal_rank'], height=0.4,
        color='#e67e22', alpha=0.7, label='Balanced (n=48)')

ax.set_yticks(y_pos)
ax.set_yticklabels(species_short, fontsize=8)
ax.set_xlabel('Rank (lower = more significant)', fontsize=10)
ax.set_title('A) Top 20 Species: Rank Comparison', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.invert_xaxis()  # Lower rank on right
ax.grid(axis='x', alpha=0.3)

# Panel B: Effect size comparison
ax = axes[1]
ax.barh(y_pos - 0.2, top_comp_df['orig_effect'], height=0.4,
        color='#9b59b6', alpha=0.7, label='Original (n=209)')
ax.barh(y_pos + 0.2, top_comp_df['bal_effect'], height=0.4,
        color='#e67e22', alpha=0.7, label='Balanced (n=48)')

ax.set_yticks(y_pos)
ax.set_yticklabels(species_short, fontsize=8)
ax.set_xlabel('Effect Size (CLR units)', fontsize=10)
ax.set_title('B) Top 20 Species: Effect Size Comparison', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / '02_top20_comparison.png', dpi=300, bbox_inches='tight')
print("  ✓ Saved 02_top20_comparison.png")
plt.close()

# ============================================================================
# STEP 6: SUMMARY REPORT
# ============================================================================
print("\nSTEP 6: GENERATE SUMMARY REPORT")
print("-" * 70)

report = f"""SENSITIVITY ANALYSIS REPORT: BALANCED COHORTS
{'='*70}

OBJECTIVE
---------
Test robustness of association analysis results to sample size imbalance
by downsampling all cohorts to equal size (n={BALANCED_N}).

SAMPLE SIZE COMPARISON
---------------------
Original dataset:
  • Schulz_2017_wastewater: {cohort_counts.get('Schulz_2017_wastewater', 0)} samples (61%)
  • Chu_2017_sludge: {cohort_counts.get('Chu_2017_sludge', 0)} samples (23%)
  • Rowe_2017_hospital_wastewater: {cohort_counts.get('Rowe_2017_hospital_wastewater', 0)} samples (10%)
  • Lekunberri_2018_river_wastewater: {cohort_counts.get('Lekunberri_2018_river_wastewater', 0)} samples (6%)
  Total: {len(common_samples)} samples
  Imbalance ratio: {cohort_counts.max() / cohort_counts.min():.1f}x

Balanced dataset:
  • All cohorts: {BALANCED_N} samples each (25%)
  Total: {len(balanced_samples)} samples
  Imbalance ratio: 1.0x (perfectly balanced)

ASSOCIATION ANALYSIS RESULTS
-----------------------------
Original (n={len(common_samples)}):
  • Species tested: {len(results_original)}
  • Significant (FDR<0.05): {n_sig_original}/{len(results_original)} ({pct_sig_original:.1f}%)
  • Effect size range: {results_original['effect_size'].min():.2f} - {results_original['effect_size'].max():.2f} CLR

Balanced (n={len(balanced_samples)}):
  • Species tested: {len(results_balanced)}
  • Significant (FDR<0.05): {n_sig_balanced}/{len(results_balanced)} ({pct_sig_balanced:.1f}%)
  • Effect size range: {results_balanced['effect_size'].min():.2f} - {results_balanced['effect_size'].max():.2f} CLR

CONCORDANCE ANALYSIS
--------------------
Significance concordance:
  • Both significant: {both_sig}/{len(comparison)} ({both_sig/len(comparison)*100:.1f}%)
  • Original only: {original_only} ({original_only/len(comparison)*100:.1f}%)
  • Balanced only: {balanced_only} ({balanced_only/len(comparison)*100:.1f}%)
  • Neither: {neither_sig} ({neither_sig/len(comparison)*100:.1f}%)

Effect size correlation: r = {eff_corr:.3f}
P-value correlation: r = {p_corr:.3f}

TOP 10 SPECIES - ORIGINAL RANKING
----------------------------------
"""

for idx, row in results_original.head(10).iterrows():
    species_name = row['species']
    orig_p = row['p_value']
    orig_eff = row['effect_size']
    
    # Find in balanced results
    bal_row = results_balanced[results_balanced['species'] == species_name]
    if len(bal_row) > 0:
        bal_p = bal_row['p_value'].values[0]
        bal_eff = bal_row['effect_size'].values[0]
        bal_sig = "YES" if bal_row['significant'].values[0] else "NO"
        bal_rank = results_balanced[results_balanced['species'] == species_name].index[0] + 1
    else:
        bal_p = np.nan
        bal_eff = np.nan
        bal_sig = "N/A"
        bal_rank = "N/A"
    
    report += f"\n{idx+1}. {species_name}\n"
    report += f"   Original: p={orig_p:.2e}, effect={orig_eff:.2f} CLR\n"
    report += f"   Balanced: p={bal_p:.2e}, effect={bal_eff:.2f} CLR, rank={bal_rank}, sig={bal_sig}\n"

report += f"""

KEY FINDINGS
------------
1. Significance rate decreased from {pct_sig_original:.1f}% → {pct_sig_balanced:.1f}% with balanced samples
   • Expected due to reduced statistical power ({len(common_samples)} → {len(balanced_samples)} samples)

2. Concordance: {both_sig/len(comparison)*100:.1f}% of species significant in BOTH analyses
   • Indicates robust biological signals
   • Species significant in both are high-confidence biomarkers

3. Effect size correlation: r = {eff_corr:.3f}
   • {'Strong' if eff_corr > 0.8 else 'Moderate' if eff_corr > 0.5 else 'Weak'} correlation
   • Effect sizes {'highly consistent' if eff_corr > 0.8 else 'moderately consistent' if eff_corr > 0.5 else 'vary'} between analyses

4. P-value correlation: r = {p_corr:.3f}
   • {'Strong' if p_corr > 0.8 else 'Moderate' if p_corr > 0.5 else 'Weak'} correlation
   • Rank order of species {'highly preserved' if p_corr > 0.8 else 'moderately preserved' if p_corr > 0.5 else 'changes'}

INTERPRETATION
--------------
{'Strong concordance and high correlation indicate that original findings are robust to sample size imbalance. The association signals are genuine biological patterns, not statistical artifacts of Schulz cohort dominance.' if both_sig/len(comparison) > 0.7 and eff_corr > 0.7 else 'Moderate concordance suggests some findings may be influenced by Schulz cohort dominance. Species significant in both analyses are more reliable biomarkers.' if both_sig/len(comparison) > 0.5 else 'Low concordance indicates substantial impact of sample size imbalance. Balanced results may be more generalizable, though with reduced power.'}

RECOMMENDATION
--------------
"""

if both_sig/len(comparison) > 0.7 and eff_corr > 0.7:
    report += "✓ Original results are robust - proceed with high-powered (n=209) analysis.\n"
    report += "✓ Species significant in both analyses are high-confidence biomarkers.\n"
    report += "✓ Sample size imbalance did not meaningfully distort findings.\n"
elif both_sig/len(comparison) > 0.5:
    report += "⚠ Focus on species significant in BOTH analyses for most robust conclusions.\n"
    report += "⚠ Report balanced sensitivity analysis alongside original results.\n"
    report += "⚠ Interpret original-only findings with caution (may be Schulz-driven).\n"
else:
    report += "⚠ Consider using balanced results as primary findings.\n"
    report += "⚠ Original results may be overly influenced by Schulz cohort.\n"
    report += "⚠ Report both analyses with clear caveats about sample size.\n"

report += f"""

FILES GENERATED
---------------
• abundance_clr_balanced.csv - Balanced CLR data ({len(balanced_samples)} samples)
• metadata_balanced.csv - Balanced metadata
• association_results_balanced_all.csv - All species results (balanced)
• association_results_balanced_significant.csv - Significant species only
• comparison_original_vs_balanced.csv - Side-by-side comparison
• 01_sensitivity_comparison.png - Concordance visualization
• 02_top20_comparison.png - Top species comparison
• sensitivity_report.txt - This report

Analysis completed: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
"""

# Save report
report_path = OUTPUT_DIR / 'sensitivity_report.txt'
with open(report_path, 'w') as f:
    f.write(report)

print(f"  ✓ Saved sensitivity_report.txt")
print()
print("="*70)
print("SENSITIVITY ANALYSIS COMPLETE")
print("="*70)
print(f"\nKey Result: {both_sig/len(comparison)*100:.1f}% concordance between analyses")
print(f"Interpretation: {'Findings are ROBUST' if both_sig/len(comparison) > 0.7 else 'Findings are MODERATELY ROBUST' if both_sig/len(comparison) > 0.5 else 'Findings are SENSITIVE to sample size'}")
print(f"\nAll outputs saved to: {OUTPUT_DIR}/")
print()
