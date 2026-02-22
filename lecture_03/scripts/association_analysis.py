#!/usr/bin/env python3
"""
Association Analysis: Linear Regression for Species x Group Associations
=========================================================================
Analysis identifies which CLR-transformed species are significantly associated
with predefined groups (wastewater studies), adjusting for confounders and
correcting for multiple testing.

Method:
- Linear regression: CLR_abundance ~ study_code + confounders
- Confounder adjustment: pH, temperature, dissolved_oxygen
- Multiple testing correction: FDR (Benjamini-Hochberg)
- P-value threshold: FDR < 0.05
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import shutil
from scipy import stats
from scipy.stats import beta, kruskal
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
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
OUTPUT_DIR = os.getenv('OUTPUT_DIR', 'output/association_analysis')

# Grouping variable (primary factor of interest)
GROUP_VAR = 'study_code'

# Sample-size normalization (lecture_02 style)
APPLY_SAMPLE_SIZE_NORMALIZATION = _env_bool('APPLY_SAMPLE_SIZE_NORMALIZATION', True)
NORMALIZE_STUDY = os.getenv('NORMALIZE_STUDY', 'Schulz_2017_wastewater')
NORMALIZED_N = int(os.getenv('NORMALIZED_N', '20'))
RANDOM_SEED = int(os.getenv('RANDOM_SEED', '42'))

# Potential confounders (will use those with sufficient non-null data)
CONFOUNDERS = ['pH', 'temperature', 'dissolved_oxygen_uM']

# Multiple testing correction
FDR_THRESHOLD = 0.05

# Visualization parameters
TOP_N_SPECIES = 12  # Top species for boxplots
MANHATTAN_THRESHOLD = 0.05  # FDR threshold for Manhattan plot

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

print(f"Abundance shape: {abundance.shape}")
print(f"Metadata shape: {metadata.shape}")

# ============================================================================
# DATA PREPARATION
# ============================================================================
print("\n" + "="*70)
print("DATA PREPARATION")
print("="*70)

# Ensure samples are aligned - remove duplicates by keeping first occurrence
metadata_clean = metadata[~metadata.index.duplicated(keep='first')].copy()
abundance_clean = abundance.copy()

# Ensure samples are aligned
common_samples = abundance_clean.index.intersection(metadata_clean.index)
abundance = abundance_clean.loc[common_samples]
metadata = metadata_clean.loc[common_samples]

print(f"Aligned samples: {len(common_samples)}")
print(f"\nGrouping variable (study_code) distribution:")
print(metadata[GROUP_VAR].value_counts())

# Sample-size normalization: downsample oversized study to match lecture_02 approach
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

print(f"\nSample-size normalization: {normalization_note}")
print("Updated group distribution:")
print(metadata[GROUP_VAR].value_counts())

# For this analysis, we'll use study_code as the primary grouping variable
# without adjusting for confounders (environmental data not available in aligned metadata)
available_confounders = []

print(f"\nNote: Environmental confounder data not available in aligned metadata.")
print(f"Analysis will test for associations with '{GROUP_VAR}' without adjustment.")

# Create analysis dataset  
analysis_data = pd.DataFrame({
    GROUP_VAR: metadata[GROUP_VAR]
})
analysis_data = pd.concat([analysis_data, abundance], axis=1)

# Drop rows with missing group variable
analysis_data = analysis_data[analysis_data[GROUP_VAR].notna()].copy()
print(f"\nFinal analysis dataset: {analysis_data.shape[0]} samples × {len(abundance.columns)} species")

# ============================================================================
# LINEAR REGRESSION ANALYSIS
# ============================================================================
print("\n" + "="*70)
print("LINEAR REGRESSION ANALYSIS FOR EACH SPECIES")
print("="*70)

species = abundance.columns.tolist()
n_species = len(species)

# Store results
results = []

print(f"\nFitting models for {n_species} species...")

for i, sp in enumerate(species):
    if (i + 1) % 50 == 0:
        print(f"  {i+1}/{n_species} species completed")
    
    try:
        # Build formula: species ~ study_code 
        # (no confounders available in aligned metadata)
        formula = f'{sp} ~ C({GROUP_VAR})'
        
        # Prepare data for this species (remove NA values)
        model_data = analysis_data[[sp, GROUP_VAR]].dropna()
        
        if len(model_data) < 10:  # Need minimum samples
            continue
        
        # Fit model
        model = ols(formula, data=model_data).fit()
        anova_table = anova_lm(model, typ=2)
        
        # Extract p-value for GROUP_VAR
        group_var_key = f'C({GROUP_VAR})'
        if group_var_key in anova_table.index:
            p_value = anova_table.loc[group_var_key, 'PR(>F)']
            f_statistic = anova_table.loc[group_var_key, 'F']

            # Non-parametric robustness check (same grouping variable)
            grouped_values = [group_df[sp].values for _, group_df in model_data.groupby(GROUP_VAR)]
            kw_statistic, p_value_kw = kruskal(*grouped_values)
            
            # Get effect sizes (mean abundance by group)
            group_means = model_data.groupby(GROUP_VAR)[sp].mean()
            effect_size = group_means.max() - group_means.min()
            
            results.append({
                'species': sp,
                'p_value': p_value,
                'f_statistic': f_statistic,
                'p_value_kw': p_value_kw,
                'kw_statistic': kw_statistic,
                'effect_size': effect_size,
                'n_samples': len(model_data),
                'r_squared': model.rsquared
            })
    except Exception as e:
        continue

# Create results dataframe
results_df = pd.DataFrame(results)
print(f"\nSuccessfully modeled {len(results_df)} species")

# ============================================================================
# MULTIPLE TESTING CORRECTION
# ============================================================================
print("\n" + "="*70)
print("MULTIPLE TESTING CORRECTION")
print("="*70)

# Apply FDR correction (Benjamini-Hochberg)
reject_fdr, p_adjusted_fdr, _, _ = multipletests(
    results_df['p_value'], 
    alpha=FDR_THRESHOLD,
    method='fdr_bh'
)

results_df['p_adjusted'] = p_adjusted_fdr
results_df['sig_fdr'] = reject_fdr

# Apply FDR correction for Kruskal-Wallis p-values
reject_kw_fdr, p_adjusted_kw_fdr, _, _ = multipletests(
    results_df['p_value_kw'],
    alpha=FDR_THRESHOLD,
    method='fdr_bh'
)
results_df['p_adjusted_kw'] = p_adjusted_kw_fdr
results_df['sig_fdr_kw'] = reject_kw_fdr

# Calculate -log10(p) for visualization
results_df['log10p'] = -np.log10(results_df['p_value'])
results_df['log10p_adjusted'] = -np.log10(results_df['p_adjusted'])

# Apply Bonferroni correction as well
reject_bonf, p_adjusted_bonf, _, _ = multipletests(
    results_df['p_value'],
    alpha=FDR_THRESHOLD,
    method='bonferroni'
)
results_df['p_adjusted_bonf'] = p_adjusted_bonf
results_df['sig_bonf'] = reject_bonf

# Sort by adjusted p-value
results_df = results_df.sort_values('p_adjusted')

# Print summary statistics
print(f"\nSignificant associations (FDR < {FDR_THRESHOLD}): {results_df['sig_fdr'].sum()}")
print(f"Significant associations (Kruskal FDR < {FDR_THRESHOLD}): {results_df['sig_fdr_kw'].sum()}")
print(f"Significant associations (Bonferroni < {FDR_THRESHOLD}): {results_df['sig_bonf'].sum()}")
print(f"\nEffect size range: [{results_df['effect_size'].min():.4f}, {results_df['effect_size'].max():.4f}]")
print(f"F-statistic range: [{results_df['f_statistic'].min():.2f}, {results_df['f_statistic'].max():.2f}]")

# Show top 10 associated species
print(f"\nTop 10 associated species (FDR-corrected):")
print(results_df[['species', 'p_value', 'p_adjusted', 'p_value_kw', 'p_adjusted_kw', 'effect_size', 'f_statistic', 'sig_fdr', 'sig_fdr_kw']].head(10).to_string(index=False))

# ============================================================================
# SAVE RESULTS
# ============================================================================
print("\n" + "="*70)
print("SAVING RESULTS")
print("="*70)

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(HISTORY_DIR, exist_ok=True)

# Save full results
results_df.to_csv(f'{OUTPUT_DIR}/association_results_all.csv', index=False)
print(f"✓ Saved full results: {OUTPUT_DIR}/association_results_all.csv")

# Save significant results only
sig_results = results_df[results_df['sig_fdr']].copy()
sig_results.to_csv(f'{OUTPUT_DIR}/association_results_significant.csv', index=False)
print(f"✓ Saved significant results: {OUTPUT_DIR}/association_results_significant.csv ({len(sig_results)} species)")

# ============================================================================
# VISUALIZATION 1: MANHATTAN PLOT
# ============================================================================
print("\n" + "="*70)
print("GENERATING VISUALIZATIONS")
print("="*70)

fig, ax = plt.subplots(figsize=(16, 8))

# Prepare data for Manhattan plot
plot_data = results_df.sort_values('p_value').reset_index(drop=True)
plot_data['position'] = range(len(plot_data))

# Color by significance
colors = ['#d62728' if sig else '#1f77b4' for sig in plot_data['sig_fdr']]

# Manhattan plot
ax.scatter(plot_data['position'], plot_data['log10p_adjusted'], 
          c=colors, alpha=0.6, s=50, edgecolors='black', linewidth=0.5)

# Add significance threshold line
sig_threshold = -np.log10(FDR_THRESHOLD)
ax.axhline(y=sig_threshold, color='red', linestyle='--', linewidth=2, 
          label=f'FDR threshold (p={FDR_THRESHOLD})')

# Highlight top significant species
top_sig = plot_data[plot_data['sig_fdr']].head(15)
for idx, row in top_sig.iterrows():
    ax.annotate(row['species'], 
               xy=(row['position'], row['log10p_adjusted']),
               xytext=(0, 10),
               textcoords='offset points',
               fontsize=8,
               ha='center',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3),
               arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', lw=0.5))

ax.set_xlabel('Species (ranked by adjusted p-value)', fontsize=12, fontweight='bold')
ax.set_ylabel('-log10(Adjusted P-value)', fontsize=12, fontweight='bold')
ax.set_title('Manhattan Plot: Species-Group Associations\n(after FDR correction, adjusting for confounders)', 
            fontsize=14, fontweight='bold', pad=20)
ax.grid(True, alpha=0.3, axis='y')
ax.legend(fontsize=11, loc='upper right')

# Add statistics box
stats_text = (
    f'Species tested: {len(results_df)}\n'
    f'Significant ANOVA FDR: {results_df["sig_fdr"].sum()}\n'
    f'Significant Kruskal FDR: {results_df["sig_fdr_kw"].sum()}\n'
    f'Confounders: {len(available_confounders)}'
)
ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
       fontsize=10, verticalalignment='top',
       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/01_manhattan_plot.png', dpi=300, bbox_inches='tight')
print(f"✓ Saved: {OUTPUT_DIR}/01_manhattan_plot.png")
plt.close()

# ============================================================================
# VISUALIZATION 2: BOXPLOTS FOR TOP SPECIES
# ============================================================================

# Get top species by smallest p-value (not just FDR-sig)
top_species = results_df.nsmallest(min(TOP_N_SPECIES, len(results_df[results_df['sig_fdr']])), 'p_value')['species'].tolist()

if len(top_species) > 0:
    n_plots = len(top_species)
    n_cols = 3
    n_rows = (n_plots + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
    if n_plots == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Define color palette for study codes
    study_codes = sorted(metadata[GROUP_VAR].unique())
    colors = sns.color_palette("husl", len(study_codes))
    color_dict = dict(zip(study_codes, colors))
    
    for idx, species_name in enumerate(top_species):
        ax = axes[idx]
        
        # Prepare data
        plot_df = pd.DataFrame({
            GROUP_VAR: metadata[GROUP_VAR],
            'CLR_abundance': abundance[species_name]
        }).dropna()
        
        # Create boxplot
        sns.boxplot(data=plot_df, x=GROUP_VAR, y='CLR_abundance', 
                   palette=color_dict, ax=ax)
        
        # Add individual points
        sns.stripplot(data=plot_df, x=GROUP_VAR, y='CLR_abundance',
                     color='black', alpha=0.4, size=4, ax=ax)
        
        # Get statistics for title
        met_row = results_df[results_df['species'] == species_name].iloc[0]
        sig_mark = '**' if met_row['sig_fdr'] else ''
        
        title = f"{species_name}{sig_mark}\np={met_row['p_value']:.2e}, FDR={met_row['p_adjusted']:.2e}"
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.set_xlabel('Study Code', fontsize=9)
        ax.set_ylabel('CLR-transformed Abundance', fontsize=9)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3, axis='y')
    
    # Hide unused subplots
    for idx in range(len(top_species), len(axes)):
        axes[idx].set_visible(False)
    
    fig.suptitle(f'Top {len(top_species)} Associated Species\n(Boxplots by Study)', 
                fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/02_boxplots_top_species.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {OUTPUT_DIR}/02_boxplots_top_species.png")
    plt.close()

# ============================================================================
# VISUALIZATION 3: Q-Q PLOT AND VOLCANO PLOT
# ============================================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Q-Q plot: observed vs expected -log10(p)
ax = axes[0]
sorted_p = np.sort(results_df['p_value'])
expected_p = (np.arange(1, len(sorted_p) + 1) - 0.5) / len(sorted_p)
expected_log10p = -np.log10(expected_p)
observed_log10p = -np.log10(sorted_p)

# 95% confidence envelope under null (uniform p-values)
n_points = len(sorted_p)
lower_q = beta.ppf(0.025, np.arange(1, n_points + 1), np.arange(n_points, 0, -1))
upper_q = beta.ppf(0.975, np.arange(1, n_points + 1), np.arange(n_points, 0, -1))
lower_log10p = -np.log10(np.clip(lower_q, 1e-320, 1))
upper_log10p = -np.log10(np.clip(upper_q, 1e-320, 1))

ax.scatter(expected_log10p, observed_log10p, alpha=0.6, s=30, edgecolors='black', linewidth=0.5)
ax.plot([0, max(expected_log10p)], [0, max(expected_log10p)], 'r--', lw=2, label='Expected')
ax.fill_between(expected_log10p, lower_log10p, upper_log10p, color='gray', alpha=0.2, label='95% null envelope')
ax.set_xlabel('Expected -log10(p)', fontsize=11, fontweight='bold')
ax.set_ylabel('Observed -log10(p)', fontsize=11, fontweight='bold')
ax.set_title('Q-Q Plot: ANOVA P-value Distribution', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

# Volcano plot: effect size vs p-value
ax = axes[1]
colors = ['#d62728' if sig else '#1f77b4' for sig in results_df['sig_fdr']]
ax.scatter(results_df['effect_size'], results_df['log10p_adjusted'], 
          c=colors, alpha=0.6, s=50, edgecolors='black', linewidth=0.5)

# Add threshold lines
ax.axhline(y=sig_threshold, color='red', linestyle='--', linewidth=1.5, alpha=0.7)

ax.set_xlabel('Effect Size (Max - Min Abundance)', fontsize=11, fontweight='bold')
ax.set_ylabel('-log10(Adjusted P-value)', fontsize=11, fontweight='bold')
ax.set_title('Volcano Plot: Effect Size vs P-value', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)

# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#d62728', alpha=0.6, label='Significant (FDR)'),
                   Patch(facecolor='#1f77b4', alpha=0.6, label='Not significant')]
ax.legend(handles=legend_elements, loc='upper right')

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/03_qq_volcano_plots.png', dpi=300, bbox_inches='tight')
print(f"✓ Saved: {OUTPUT_DIR}/03_qq_volcano_plots.png")
plt.close()

# ============================================================================
# VISUALIZATION 4: EFFECT SIZES AND STATISTICS
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Histogram of p-values
ax = axes[0, 0]
ax.hist(results_df['p_value'], bins=50, alpha=0.7, color='steelblue', edgecolor='black')
ax.axvline(x=FDR_THRESHOLD, color='red', linestyle='--', linewidth=2, label=f'Threshold = {FDR_THRESHOLD}')
ax.set_xlabel('P-value', fontsize=11, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
ax.set_title('Distribution of Raw P-values', fontsize=12, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3, axis='y')

# Histogram of effect sizes
ax = axes[0, 1]
ax.hist(results_df['effect_size'], bins=50, alpha=0.7, color='coral', edgecolor='black')
ax.set_xlabel('Effect Size (Max - Min Abundance)', fontsize=11, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
ax.set_title('Distribution of Effect Sizes', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')

# P-value vs Effect Size scatter
ax = axes[1, 0]
sig_df = results_df[results_df['sig_fdr']]
ax.scatter(results_df['effect_size'], -np.log10(results_df['p_value']), 
          alpha=0.5, s=30, c='lightblue', edgecolors='black', linewidth=0.5, label='Not significant')
ax.scatter(sig_df['effect_size'], -np.log10(sig_df['p_value']), 
          alpha=0.8, s=60, c='red', edgecolors='darkred', linewidth=1, label='Significant (FDR)')
ax.set_xlabel('Effect Size', fontsize=11, fontweight='bold')
ax.set_ylabel('-log10(P-value)', fontsize=11, fontweight='bold')
ax.set_title('Effect Size vs Raw P-value', fontsize=12, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

# R-squared distribution
ax = axes[1, 1]
ax.hist(results_df['r_squared'], bins=50, alpha=0.7, color='mediumseagreen', edgecolor='black')
ax.set_xlabel('R-squared', fontsize=11, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
ax.set_title('Distribution of Model R-squared', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')

fig.suptitle('Statistical Summary Plots', fontsize=14, fontweight='bold', y=1.00)
plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/04_statistical_summary.png', dpi=300, bbox_inches='tight')
print(f"✓ Saved: {OUTPUT_DIR}/04_statistical_summary.png")
plt.close()

# ============================================================================
# SUMMARY REPORT
# ============================================================================
print("\n" + "="*70)
print("ANALYSIS SUMMARY")
print("="*70)

report = f"""
ASSOCIATION ANALYSIS: SPECIES-GROUP ASSOCIATIONS
{'='*70}

OBJECTIVE:
Identify which CLR-transformed microbial species are significantly associated
with predefined groups (wastewater studies), adjusting for confounders.

DATA:
- Total samples analyzed: {len(analysis_data)}
- Species tested: {len(results_df)}
- Groups ({GROUP_VAR}): {', '.join(sorted(metadata[GROUP_VAR].unique()))}
- Sample-size normalization: {normalization_note}
- Sample sizes by group:
{metadata[GROUP_VAR].value_counts().to_string().replace(chr(10), chr(10) + '  ')}

CONFOUNDERS ADJUSTED:
(None - environmental data not available in aligned metadata)

STATISTICAL METHOD:
- Model: Linear regression (Type II ANOVA)
- Formula: CLR_abundance ~ {GROUP_VAR}
- Robustness method: Kruskal-Wallis test (same grouping)
- Multiple testing correction: FDR (Benjamini-Hochberg)
- Significance threshold: FDR < {FDR_THRESHOLD}

RESULTS:
- Total species associated with group (p < 0.05): {(results_df['p_value'] < 0.05).sum()}
- Significant after FDR correction (FDR < {FDR_THRESHOLD}): {results_df['sig_fdr'].sum()}
- Significant after Kruskal-Wallis + FDR: {results_df['sig_fdr_kw'].sum()}
- Significant after Bonferroni (p < {FDR_THRESHOLD}): {results_df['sig_bonf'].sum()}

EFFECT SIZES:
- Mean effect size (Max-Min): {results_df['effect_size'].mean():.4f}
- Median effect size: {results_df['effect_size'].median():.4f}
- Range: [{results_df['effect_size'].min():.4f}, {results_df['effect_size'].max():.4f}]

MODEL FIT:
- Mean R-squared: {results_df['r_squared'].mean():.4f}
- Median R-squared: {results_df['r_squared'].median():.4f}

TOP 10 ASSOCIATED SPECIES (FDR-CORRECTED):

"""

for rank, (_, row) in enumerate(results_df.head(10).iterrows(), start=1):
    report += f"\n{rank}. {row['species']}"
    report += f"\n   p-value: {row['p_value']:.2e}"
    report += f"\n   FDR-adjusted p: {row['p_adjusted']:.2e}"
    report += f"\n   Kruskal p-value: {row['p_value_kw']:.2e}"
    report += f"\n   Kruskal FDR-adjusted p: {row['p_adjusted_kw']:.2e}"
    report += f"\n   Effect size: {row['effect_size']:.4f}"
    report += f"\n   F-statistic: {row['f_statistic']:.2f}"
    report += f"\n   R-squared: {row['r_squared']:.4f}\n"

report += f"""

OUTPUT FILES:
✓ association_results_all.csv - All {len(results_df)} species with statistics
✓ association_results_significant.csv - FDR-significant species only ({len(sig_results)} species)
✓ 01_manhattan_plot.png - Manhattan plot showing associations
✓ 02_boxplots_top_species.png - Boxplots for top associated species
✓ 03_qq_volcano_plots.png - Q-Q and volcano plots
✓ 04_statistical_summary.png - Histograms and scatter plots

INTERPRETATION:
Species with small FDR-adjusted p-values show significant associations with
the grouping variable (wastewater source) after adjusting for confounders.
Effect size indicates the magnitude of abundance difference between groups.

{'='*70}
"""

print(report)

# Save report
report_path = f'{OUTPUT_DIR}/association_analysis_report.txt'
with open(report_path, 'w') as f:
    f.write(report)
print(f"✓ Saved report: {OUTPUT_DIR}/association_analysis_report.txt")
snapshot_report(report_path)

print("\n" + "="*70)
print("ANALYSIS COMPLETE!")
print("="*70)
