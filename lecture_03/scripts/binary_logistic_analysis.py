#!/usr/bin/env python3
"""
Binary Presence/Absence Analysis: Logistic Regression for Sparse Taxa
======================================================================
Analyzes species presence/absence (binary data) to identify associations with
wastewater groups, accounting for sparse taxa that may have zero inflation.

Method:
- Convert CLR-transformed abundance to binary (present/absent)
- Logistic regression: presence ~ study_code
- Multiple testing correction: FDR (Benjamini-Hochberg)
- Comparison with linear regression results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import shutil
from scipy import stats
from scipy.special import expit, logit
from sklearn.linear_model import LogisticRegression
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod import families
from statsmodels.formula.api import glm
from statsmodels.stats.multitest import multipletests
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================
def _env_bool(name, default):
    value = os.getenv(name)
    if value is None:
        return default
    return value.strip().lower() in {'1', 'true', 'yes', 'on'}


ABUNDANCE_FILE = os.getenv('ABUNDANCE_FILE', 'output/filtering_clr_analysis/abundance_clr_aligned.csv')
METADATA_FILE = os.getenv('METADATA_FILE', 'output/filtering_clr_analysis/metadata_aligned.csv')
LINEAR_RESULTS_FILE = os.getenv('LINEAR_RESULTS_FILE', 'output/association_analysis/association_results_all.csv')
OUTPUT_DIR = os.getenv('OUTPUT_DIR', 'output/binary_logistic_analysis')

GROUP_VAR = 'study_code'

# Sample-size normalization (lecture_02 style)
APPLY_SAMPLE_SIZE_NORMALIZATION = _env_bool('APPLY_SAMPLE_SIZE_NORMALIZATION', True)
NORMALIZE_STUDY = os.getenv('NORMALIZE_STUDY', 'Schulz_2017_wastewater')
NORMALIZED_N = int(os.getenv('NORMALIZED_N', '20'))
RANDOM_SEED = int(os.getenv('RANDOM_SEED', '42'))

# Multiple testing correction
FDR_THRESHOLD = 0.05

# Visualization parameters
TOP_N_SPECIES = 12

RUN_STAMP = datetime.now().strftime('%Y%m%d_%H%M%S')
HISTORY_DIR = f'{OUTPUT_DIR}/report_history'


def snapshot_report(file_path):
    if os.path.exists(file_path):
        base_name = os.path.basename(file_path)
        stem, ext = os.path.splitext(base_name)
        archived = f'{HISTORY_DIR}/{stem}_{RUN_STAMP}{ext}'
        shutil.copy2(file_path, archived)
        print(f"✓ Archived report: {archived}")

# ============================================================================
# LOAD DATA
# ============================================================================
print("Loading data...")
abundance = pd.read_csv(ABUNDANCE_FILE, index_col=0)
metadata = pd.read_csv(METADATA_FILE, index_col=0)
linear_results = pd.read_csv(LINEAR_RESULTS_FILE)

# Clean metadata (remove duplicates)
metadata_clean = metadata[~metadata.index.duplicated(keep='first')].copy()

# Align samples
common_samples = abundance.index.intersection(metadata_clean.index)
abundance = abundance.loc[common_samples]
metadata = metadata_clean.loc[common_samples]

print(f"Aligned samples: {len(common_samples)}")

# Sample-size normalization to reduce cohort imbalance
normalization_note = "No sample-size normalization applied"
if APPLY_SAMPLE_SIZE_NORMALIZATION and GROUP_VAR == 'study_code':
    target_mask = metadata[GROUP_VAR] == NORMALIZE_STUDY
    target_n = int(target_mask.sum())
    if target_n > NORMALIZED_N:
        normalized_subset = metadata[target_mask].sample(n=NORMALIZED_N, random_state=RANDOM_SEED)
        metadata = pd.concat([metadata[~target_mask], normalized_subset]).sort_index()
        abundance = abundance.loc[metadata.index]
        normalization_note = f"Downsampled {NORMALIZE_STUDY}: {target_n} -> {NORMALIZED_N}"
    else:
        normalization_note = (
            f"Normalization skipped for {NORMALIZE_STUDY} "
            f"(available: {target_n}, threshold: {NORMALIZED_N})"
        )

print(f"Sample-size normalization: {normalization_note}")
print("Updated group distribution:")
print(metadata[GROUP_VAR].value_counts())

# ============================================================================
# CONVERT TO BINARY PRESENCE/ABSENCE
# ============================================================================
print("\n" + "="*70)
print("CONVERTING TO BINARY PRESENCE/ABSENCE DATA")
print("="*70)

# Use detection threshold: species present if CLR > -inf (i.e., any non-negative value in original space)
# In CLR space, we use a more conservative threshold: abundance > -1 (very weak presence)
binary_abundance = (abundance > -1).astype(int)

print(f"\nBinary abundance table shape: {binary_abundance.shape}")

# Calculate prevalence (% of samples with presence)
prevalence = binary_abundance.mean() * 100

print(f"\nPrevalence statistics:")
print(f"  Mean prevalence: {prevalence.mean():.1f}%")
print(f"  Median prevalence: {prevalence.median():.1f}%")
print(f"  Range: {prevalence.min():.1f}% - {prevalence.max():.1f}%")

# Identify sparse species
sparse_threshold = 10  # Present in < 10% of samples
sparse_species = prevalence[prevalence < sparse_threshold].index.tolist()
print(f"\nSparse species (< {sparse_threshold}% prevalence): {len(sparse_species)} out of {len(binary_abundance.columns)}")

# Report group-specific prevalence
print(f"\nPrevalence by group:")
for group in metadata[GROUP_VAR].unique():
    mask = metadata[GROUP_VAR] == group
    group_binary = binary_abundance[mask]
    group_prev = group_binary.mean() * 100
    print(f"  {group}: {group_prev.mean():.1f}% mean prevalence")

# ============================================================================
# LOGISTIC REGRESSION ANALYSIS
# ============================================================================
print("\n" + "="*70)
print("LOGISTIC REGRESSION ANALYSIS FOR EACH SPECIES")
print("="*70)

species = binary_abundance.columns.tolist()
n_species = len(species)

# Store results
results = []

print(f"\nFitting logistic models for {n_species} species...")

for i, sp in enumerate(species):
    if (i + 1) % 50 == 0:
        print(f"  {i+1}/{n_species} species completed")
    
    try:
        # Prepare data for this species
        model_data = pd.DataFrame({
            GROUP_VAR: metadata[GROUP_VAR],
            'presence': binary_abundance[sp]
        }).dropna()
        
        if len(model_data) < 10:  # Need minimum samples
            continue
        
        # Check prevalence - skip if no variation (all 0 or all 1)
        if model_data['presence'].std() == 0:
            continue
        
        # Fit logistic regression model
        # Using GLM with binomial family for more detailed statistics
        formula = f'presence ~ C({GROUP_VAR})'
        model = glm(formula, data=model_data, family=families.Binomial()).fit()
        
        # Extract summary statistics
        p_values = model.pvalues
        
        # Get p-value for group effect (average of all group comparisons)
        group_params = [p for i, p in enumerate(p_values) if f'C({GROUP_VAR})' in model.model.exog_names[i]]
        
        if len(group_params) == 0:
            continue
        
        # Use minimum p-value from group comparisons (most significant)
        p_value = min(group_params) if len(group_params) > 0 else 1.0
        
        # Calculate odds ratios for each group
        coefs = model.params
        odds_ratios = np.exp(coefs)
        
        # Calculate effect size: difference in prevalence between most/least prevalent groups
        prev_by_group = model_data.groupby(GROUP_VAR)['presence'].agg(['mean', 'sum', 'count'])
        effect_size = prev_by_group['mean'].max() - prev_by_group['mean'].min()
        
        # Calculate prevalence in each group
        max_prev_group = prev_by_group['mean'].idxmax()
        min_prev_group = prev_by_group['mean'].idxmin()
        
        results.append({
            'species': sp,
            'p_value': p_value,
            'aic': model.aic,
            'bic': model.bic,
            'effect_size': effect_size,
            'max_prevalence': prev_by_group['mean'].max(),
            'max_prev_group': max_prev_group,
            'min_prevalence': prev_by_group['mean'].min(),
            'min_prev_group': min_prev_group,
            'n_samples': len(model_data),
            'overall_prevalence': model_data['presence'].mean()
        })
    except Exception as e:
        continue

# Create results dataframe
results_df = pd.DataFrame(results)
print(f"\nSuccessfully modeled {len(results_df)} species with logistic regression")

# ============================================================================
# MULTIPLE TESTING CORRECTION
# ============================================================================
print("\n" + "="*70)
print("MULTIPLE TESTING CORRECTION")
print("="*70)

# Apply FDR correction
reject_fdr, p_adjusted_fdr, _, _ = multipletests(
    results_df['p_value'],
    alpha=FDR_THRESHOLD,
    method='fdr_bh'
)

results_df['p_adjusted'] = p_adjusted_fdr
results_df['sig_fdr'] = reject_fdr

# Calculate -log10 for visualization
results_df['log10p'] = -np.log10(results_df['p_value'])
results_df['log10p_adjusted'] = -np.log10(results_df['p_adjusted'])

# Sort by adjusted p-value
results_df = results_df.sort_values('p_adjusted')

print(f"\nSignificant associations (FDR < {FDR_THRESHOLD}): {results_df['sig_fdr'].sum()}")
print(f"Total species analyzed: {len(results_df)}")

print(f"\nEffect size statistics (prevalence difference):")
print(f"  Mean: {results_df['effect_size'].mean():.4f}")
print(f"  Median: {results_df['effect_size'].median():.4f}")
print(f"  Range: [{results_df['effect_size'].min():.4f}, {results_df['effect_size'].max():.4f}]")

# Show top 10 associated species
print(f"\nTop 10 associated species (FDR-corrected):")
print(results_df[['species', 'p_value', 'p_adjusted', 'effect_size', 'max_prev_group', 'min_prev_group', 'sig_fdr']].head(10).to_string(index=False))

# ============================================================================
# MERGE WITH LINEAR REGRESSION RESULTS
# ============================================================================
print("\n" + "="*70)
print("COMPARING LINEAR vs LOGISTIC REGRESSION")
print("="*70)

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(HISTORY_DIR, exist_ok=True)

# Merge results
comparison_df = results_df[['species', 'p_value', 'p_adjusted', 'effect_size']].copy()
comparison_df.rename(columns={
    'p_value': 'log_p_value',
    'p_adjusted': 'log_p_adjusted',
    'effect_size': 'log_effect_size'
}, inplace=True)

# Merge with linear results
linear_subset = linear_results[['species', 'p_value', 'p_adjusted', 'effect_size']].copy()
linear_subset.rename(columns={
    'p_value': 'lin_p_value',
    'p_adjusted': 'lin_p_adjusted',
    'effect_size': 'lin_effect_size'
}, inplace=True)

comparison_df = comparison_df.merge(linear_subset, on='species', how='inner')

# Calculate correlation between methods
corr_p = np.corrcoef(comparison_df['log_p_value'], comparison_df['lin_p_value'])[0, 1]
corr_effect = np.corrcoef(comparison_df['log_effect_size'], comparison_df['lin_effect_size'])[0, 1]

print(f"\nCorrelation between methods:")
print(f"  P-values: r = {corr_p:.4f}")
print(f"  Effect sizes: r = {corr_effect:.4f}")

# Identify concordance
comparison_df['both_sig_fdr'] = (comparison_df['log_p_adjusted'] < 0.05) & (comparison_df['lin_p_adjusted'] < 0.05)
comparison_df['log_only'] = (comparison_df['log_p_adjusted'] < 0.05) & (comparison_df['lin_p_adjusted'] >= 0.05)
comparison_df['lin_only'] = (comparison_df['log_p_adjusted'] >= 0.05) & (comparison_df['lin_p_adjusted'] < 0.05)

print(f"\nSignificance concordance:")
print(f"  Both significant: {comparison_df['both_sig_fdr'].sum()}")
print(f"  Logistic only: {comparison_df['log_only'].sum()}")
print(f"  Linear only: {comparison_df['lin_only'].sum()}")

# Save merged comparison
comparison_df.to_csv(f'{OUTPUT_DIR}/comparison_linear_vs_logistic.csv', index=False)
print(f"\n✓ Saved comparison: {OUTPUT_DIR}/comparison_linear_vs_logistic.csv")

# ============================================================================
# SAVE RESULTS
# ============================================================================
print("\n" + "="*70)
print("SAVING RESULTS")
print("="*70)

# Save full logistic results
results_df.to_csv(f'{OUTPUT_DIR}/logistic_regression_results_all.csv', index=False)
print(f"✓ Saved all results: {OUTPUT_DIR}/logistic_regression_results_all.csv")

# Save significant results
sig_results = results_df[results_df['sig_fdr']].copy()
sig_results.to_csv(f'{OUTPUT_DIR}/logistic_regression_results_significant.csv', index=False)
print(f"✓ Saved significant results: {OUTPUT_DIR}/logistic_regression_results_significant.csv ({len(sig_results)} species)")

# ============================================================================
# VISUALIZATION 1: VOLCANO PLOT (BINARY)
# ============================================================================
print("\n" + "="*70)
print("GENERATING VISUALIZATIONS")
print("="*70)

fig, ax = plt.subplots(figsize=(12, 8))

# Color by significance
colors = ['#d62728' if sig else '#1f77b4' for sig in results_df['sig_fdr']]

# Volcano plot: effect size vs p-value
ax.scatter(results_df['effect_size'], results_df['log10p_adjusted'],
          c=colors, alpha=0.6, s=50, edgecolors='black', linewidth=0.5)

# Add significance threshold line
sig_threshold = -np.log10(FDR_THRESHOLD)
ax.axhline(y=sig_threshold, color='red', linestyle='--', linewidth=2,
          label=f'FDR threshold (p={FDR_THRESHOLD})')

# Highlight top significant species
top_sig = results_df[results_df['sig_fdr']].head(10)
for idx, row in top_sig.iterrows():
    ax.annotate(row['species'],
               xy=(row['effect_size'], row['log10p_adjusted']),
               xytext=(0, 10),
               textcoords='offset points',
               fontsize=8,
               ha='center',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3),
               arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', lw=0.5))

ax.set_xlabel('Effect Size (Prevalence Difference)', fontsize=12, fontweight='bold')
ax.set_ylabel('-log10(Adjusted P-value)', fontsize=12, fontweight='bold')
ax.set_title('Volcano Plot: Binary Presence/Absence Analysis\n(Logistic Regression with FDR Correction)',
            fontsize=14, fontweight='bold', pad=20)
ax.grid(True, alpha=0.3, axis='y')
ax.legend(fontsize=11, loc='upper right')

# Add statistics box
stats_text = f'Species tested: {len(results_df)}\nSignificant (FDR): {results_df["sig_fdr"].sum()}\nMethod: Logistic Regression'
stats_text = (
    f'Species tested: {len(results_df)}\n'
    f'Significant (FDR): {results_df["sig_fdr"].sum()}\n'
    f'Method: Logistic Regression\n'
    f'{normalization_note}'
)
ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
       fontsize=10, verticalalignment='top',
       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/01_volcano_plot_binary.png', dpi=300, bbox_inches='tight')
print(f"✓ Saved: {OUTPUT_DIR}/01_volcano_plot_binary.png")
plt.close()

# ============================================================================
# VISUALIZATION 2: COMPARISON SCATTER PLOTS
# ============================================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# P-value comparison
ax = axes[0]
ax.scatter(comparison_df['lin_p_value'], comparison_df['log_p_value'],
          alpha=0.5, s=30, edgecolors='black', linewidth=0.5)
ax.plot([0, 1], [0, 1], 'r--', lw=2, label='Perfect agreement')
ax.set_xlabel('Linear Regression -log10(p)', fontsize=11, fontweight='bold')
ax.set_ylabel('Logistic Regression -log10(p)', fontsize=11, fontweight='bold')
ax.set_title('P-value Comparison\n(r={:.3f})'.format(corr_p), fontsize=12, fontweight='bold')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True, alpha=0.3)
ax.legend()

# Effect size comparison
ax = axes[1]
ax.scatter(comparison_df['lin_effect_size'], comparison_df['log_effect_size'],
          alpha=0.5, s=30, edgecolors='black', linewidth=0.5)
max_val = max(comparison_df['lin_effect_size'].max(), comparison_df['log_effect_size'].max())
ax.plot([0, max_val], [0, max_val], 'r--', lw=2, label='Perfect agreement')
ax.set_xlabel('Linear: CLR Effect Size', fontsize=11, fontweight='bold')
ax.set_ylabel('Logistic: Prevalence Difference', fontsize=11, fontweight='bold')
ax.set_title('Effect Size Comparison\n(r={:.3f})'.format(corr_effect), fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/02_method_comparison.png', dpi=300, bbox_inches='tight')
print(f"✓ Saved: {OUTPUT_DIR}/02_method_comparison.png")
plt.close()

# ============================================================================
# VISUALIZATION 3: CONCORDANCE/DISCORDANCE
# ============================================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Venn diagram equivalent
ax = axes[0]
categories = ['Both\nSignificant', 'Logistic\nOnly', 'Linear\nOnly', 'Neither']
values = [
    comparison_df['both_sig_fdr'].sum(),
    comparison_df['log_only'].sum(),
    comparison_df['lin_only'].sum(),
    len(comparison_df) - comparison_df['both_sig_fdr'].sum() - comparison_df['log_only'].sum() - comparison_df['lin_only'].sum()
]
colors_venn = ['#2ecc71', '#f39c12', '#e74c3c', '#95a5a6']
bars = ax.bar(categories, values, color=colors_venn, edgecolor='black', linewidth=1.5)
ax.set_ylabel('Number of Species', fontsize=11, fontweight='bold')
ax.set_title('Significance Concordance (FDR < 0.05)', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')

# Add value labels on bars
for bar, val in zip(bars, values):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
           f'{int(val)}', ha='center', va='bottom', fontweight='bold')

# Manhattan plot comparison
ax = axes[1]
plot_data = comparison_df.sort_values('lin_p_value').reset_index(drop=True)
plot_data['position'] = range(len(plot_data))

# Plot both methods
ax.scatter(plot_data['position'], -np.log10(plot_data['lin_p_value']),
          color='#3498db', alpha=0.6, s=30, label='Linear', edgecolors='black', linewidth=0.3)
ax.scatter(plot_data['position'], -np.log10(plot_data['log_p_value']),
          color='#e74c3c', alpha=0.6, s=20, marker='^', label='Logistic', edgecolors='black', linewidth=0.3)

sig_threshold = -np.log10(0.05)
ax.axhline(y=sig_threshold, color='gray', linestyle='--', linewidth=1.5, alpha=0.7, label=f'FDR threshold')

ax.set_xlabel('Species (ranked by linear p-value)', fontsize=11, fontweight='bold')
ax.set_ylabel('-log10(Adjusted P-value)', fontsize=11, fontweight='bold')
ax.set_title('Method Comparison: P-value Rankings', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')
ax.legend()

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/03_concordance_analysis.png', dpi=300, bbox_inches='tight')
print(f"✓ Saved: {OUTPUT_DIR}/03_concordance_analysis.png")
plt.close()

# ============================================================================
# VISUALIZATION 4: PREVALENCE BY GROUP (TOP SPECIES)
# ============================================================================

top_species_list = results_df.nsmallest(TOP_N_SPECIES, 'p_value')['species'].tolist()

if len(top_species_list) > 0:
    n_plots = len(top_species_list)
    n_cols = 3
    n_rows = (n_plots + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
    if n_plots == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for idx, species_name in enumerate(top_species_list):
        ax = axes[idx]

        # Get prevalence by group
        group_data = []
        for group in sorted(metadata[GROUP_VAR].unique()):
            mask = metadata[GROUP_VAR] == group
            prev = binary_abundance[mask][species_name].mean() * 100
            group_data.append({'group': group, 'prevalence': prev})

        group_df = pd.DataFrame(group_data)

        # Create bar plot
        bars = ax.bar(group_df['group'], group_df['prevalence'],
                     color='steelblue', edgecolor='black', linewidth=1)

        # Get statistics
        met_row = results_df[results_df['species'] == species_name].iloc[0]
        sig_mark = '**' if met_row['sig_fdr'] else ''

        title = f"{species_name}{sig_mark}\np={met_row['p_value']:.2e}"
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.set_ylabel('Prevalence (%)', fontsize=9)
        ax.set_xlabel('Wastewater Source', fontsize=9)
        ax.tick_params(axis='x', rotation=45)
        ax.set_ylim([0, 100])
        ax.grid(True, alpha=0.3, axis='y')

        # Add percentage labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.1f}%', ha='center', va='bottom', fontsize=8)

    # Hide unused subplots
    for idx in range(len(top_species_list), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle(f'Top {len(top_species_list)} Binary Associations: Prevalence by Group',
                fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/04_prevalence_by_group.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {OUTPUT_DIR}/04_prevalence_by_group.png")
    plt.close()

# ============================================================================
# SUMMARY REPORT
# ============================================================================
print("\n" + "="*70)
print("ANALYSIS SUMMARY")
print("="*70)

report = f"""
BINARY PRESENCE/ABSENCE ANALYSIS: LOGISTIC REGRESSION
{'='*70}

OBJECTIVE:
Identify which microbial species are significantly associated with wastewater
source using binary presence/absence (logistic) models, accounting for sparse
taxa with high zero-inflation.

DATA:
- Total samples analyzed: {len(binary_abundance)}
- Species tested: {len(results_df)}
- Binary conversion threshold: CLR > -1 (any detectable presence)
- Groups ({GROUP_VAR}): {', '.join(sorted(metadata[GROUP_VAR].unique()))}
- Sample-size normalization: {normalization_note}

BINARY DATA CHARACTERISTICS:
- Mean species prevalence: {prevalence.mean():.1f}%
- Median species prevalence: {prevalence.median():.1f}%
- Range: {prevalence.min():.1f}% - {prevalence.max():.1f}%
- Sparse species (< 10% prevalence): {len(sparse_species)} species

STATISTICAL METHOD:
- Model: Logistic regression (binomial GLM)
- Formula: presence ~ {GROUP_VAR}
- Multiple testing correction: FDR (Benjamini-Hochberg)
- Significance threshold: FDR < {FDR_THRESHOLD}

RESULTS (LOGISTIC REGRESSION):
- Total species with logistic models: {len(results_df)}
- Significant associations (FDR < {FDR_THRESHOLD}): {results_df['sig_fdr'].sum()}
- Effect size range (prevalence difference): [{results_df['effect_size'].min():.4f}, {results_df['effect_size'].max():.4f}]
- Mean effect size: {results_df['effect_size'].mean():.4f}

COMPARISON WITH LINEAR REGRESSION:
- Correlation (p-values): r = {corr_p:.4f}
- Correlation (effect sizes): r = {corr_effect:.4f}
- Both methods significant: {comparison_df['both_sig_fdr'].sum()} species
- Logistic only: {comparison_df['log_only'].sum()} species
- Linear only: {comparison_df['lin_only'].sum()} species

TOP 10 ASSOCIATED SPECIES (LOGISTIC):

"""

for i, row in results_df.head(10).iterrows():
    report += f"\n{i+1}. {row['species']}"
    report += f"\n   p-value: {row['p_value']:.2e}"
    report += f"\n   FDR-adjusted p: {row['p_adjusted']:.2e}"
    report += f"\n   Prevalence effect size: {row['effect_size']:.4f}"
    report += f"\n   Enriched in: {row['max_prev_group']} ({row['max_prevalence']:.1%} prevalence)"
    report += f"\n   Depleted in: {row['min_prev_group']} ({row['min_prevalence']:.1%} prevalence)\n"

report += f"""

OUTPUT FILES:
✓ logistic_regression_results_all.csv - All {len(results_df)} species with statistics
✓ logistic_regression_results_significant.csv - FDR-significant species ({len(sig_results)})
✓ comparison_linear_vs_logistic.csv - Side-by-side comparison of methods
✓ 01_volcano_plot_binary.png - Volcano plot for binary analysis
✓ 02_method_comparison.png - Scatter plots comparing methods
✓ 03_concordance_analysis.png - Significance concordance and ranking comparison
✓ 04_prevalence_by_group.png - Prevalence of top species by group

INTERPRETATION:
Binary models identify species that show consistent presence/absence patterns
across groups, potentially capturing ecological niches better than abundance
alone. Sparse species may show stronger associations in binary models.

{'='*70}
"""

print(report)

# Save report
report_path = f'{OUTPUT_DIR}/logistic_regression_report.txt'
with open(report_path, 'w') as f:
    f.write(report)
print(f"✓ Saved report: {OUTPUT_DIR}/logistic_regression_report.txt")
snapshot_report(report_path)

print("\n" + "="*70)
print("BINARY ANALYSIS COMPLETE!")
print("="*70)
