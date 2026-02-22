#!/usr/bin/env python3
"""
Pairwise Comparisons for MWAS Analysis
Performs post-hoc pairwise tests between all group pairs
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import shutil
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests
from itertools import combinations
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Configuration
def _env_bool(name, default):
    value = os.getenv(name)
    if value is None:
        return default
    return value.strip().lower() in {'1', 'true', 'yes', 'on'}


INPUT_ABUNDANCE = os.getenv('INPUT_ABUNDANCE', 'output/filtering_clr_analysis/abundance_clr_aligned.csv')
INPUT_METADATA = os.getenv('INPUT_METADATA', 'output/filtering_clr_analysis/metadata_aligned.csv')
OUTPUT_DIR = os.getenv('OUTPUT_DIR', 'output/pairwise_comparisons')
FDR_THRESHOLD = 0.05
TOP_N_SPECIES = 20

# Sample-size normalization (lecture_02 style)
APPLY_SAMPLE_SIZE_NORMALIZATION = _env_bool('APPLY_SAMPLE_SIZE_NORMALIZATION', True)
NORMALIZE_STUDY = os.getenv('NORMALIZE_STUDY', 'Schulz_2017_wastewater')
NORMALIZED_N = int(os.getenv('NORMALIZED_N', '20'))
RANDOM_SEED = int(os.getenv('RANDOM_SEED', '42'))

RUN_STAMP = datetime.now().strftime('%Y%m%d_%H%M%S')
HISTORY_DIR = f'{OUTPUT_DIR}/report_history'


def snapshot_report(file_path):
    if os.path.exists(file_path):
        base_name = os.path.basename(file_path)
        stem, ext = os.path.splitext(base_name)
        archived = f'{HISTORY_DIR}/{stem}_{RUN_STAMP}{ext}'
        shutil.copy2(file_path, archived)
        print(f"  ✓ Archived report: {archived}")

print("="*80)
print("PAIRWISE COMPARISONS ANALYSIS")
print("="*80)

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(HISTORY_DIR, exist_ok=True)

# ============================================================================
# 1. LOAD PREPARED DATA
# ============================================================================
print("\n[1/6] Loading prepared data...")

# Load CLR-transformed abundance data (already filtered and aligned)
abundance_clr = pd.read_csv(INPUT_ABUNDANCE, index_col=0)

# Transpose if needed (species should be rows, samples should be columns)
if abundance_clr.shape[0] < abundance_clr.shape[1]:
    print(f"  Transposing abundance matrix...")
    abundance_clr = abundance_clr.T

print(f"  Loaded CLR abundance: {abundance_clr.shape[0]} species × {abundance_clr.shape[1]} samples")

# Load aligned metadata
metadata = pd.read_csv(INPUT_METADATA, index_col=0)
print(f"  Loaded metadata: {metadata.shape[0]} samples")

# Remove duplicate indices (keep first occurrence)
if metadata.index.duplicated().any():
    n_dup = metadata.index.duplicated().sum()
    print(f"  Removing {n_dup} duplicate samples...")
    metadata = metadata[~metadata.index.duplicated(keep='first')]
    print(f"  After deduplication: {metadata.shape[0]} samples")

# Align samples
common_samples = list(abundance_clr.columns.intersection(metadata.index))
print(f"  Common samples: {len(common_samples)}")

# Subset both
abundance_clr = abundance_clr[common_samples]
metadata = metadata.loc[common_samples]
print(f"  ✓ Aligned to {len(common_samples)} samples")
print(f"  DEBUG: metadata shape after alignment: {metadata.shape}")
print(f"  DEBUG: abundance_clr shape after alignment: {abundance_clr.shape}")

# Sample-size normalization to reduce cohort imbalance
normalization_note = "No sample-size normalization applied"
if APPLY_SAMPLE_SIZE_NORMALIZATION:
    target_mask = metadata['study_code'] == NORMALIZE_STUDY
    target_n = int(target_mask.sum())
    if target_n > NORMALIZED_N:
        normalized_subset = metadata[target_mask].sample(n=NORMALIZED_N, random_state=RANDOM_SEED)
        metadata = pd.concat([metadata[~target_mask], normalized_subset]).sort_index()
        abundance_clr = abundance_clr[metadata.index]
        normalization_note = f"Downsampled {NORMALIZE_STUDY}: {target_n} -> {NORMALIZED_N}"
    else:
        normalization_note = (
            f"Normalization skipped for {NORMALIZE_STUDY} "
            f"(available: {target_n}, threshold: {NORMALIZED_N})"
        )

print(f"  Sample-size normalization: {normalization_note}")

# Get groups
groups = metadata['study_code'].values
print(f"  DEBUG: groups.shape = {groups.shape}")
group_names = sorted(metadata['study_code'].unique())
n_groups = len(group_names)
print(f"  Groups: {group_names}")
print(f"  Group sizes: {metadata['study_code'].value_counts().to_dict()}")

# ============================================================================
# 2. PAIRWISE COMPARISONS
# ============================================================================
print(f"\n[2/6] Performing pairwise comparisons...")

# Generate all pairs
pairs = list(combinations(group_names, 2))
n_pairs = len(pairs)
print(f"  Testing {n_pairs} pairwise comparisons per species")
print(f"  Pairs: {pairs}")

# Store results
pairwise_results = []

for species_idx, species in enumerate(abundance_clr.index):
    if (species_idx + 1) % 100 == 0:
        print(f"  Progress: {species_idx + 1}/{len(abundance_clr.index)} species...")
    
    # Get CLR values for this species
    clr_values = abundance_clr.loc[species].values
    
    # Perform all pairwise t-tests
    for group1, group2 in pairs:
        # Get values for each group
        mask1 = metadata['study_code'] == group1
        mask2 = metadata['study_code'] == group2
        
        values1 = clr_values[mask1]
        values2 = clr_values[mask2]
        
        # Calculate statistics
        mean1 = np.mean(values1)
        mean2 = np.mean(values2)
        effect_size = mean1 - mean2  # Difference in CLR units
        
        # Welch t-test (reference)
        t_stat, p_value_t = stats.ttest_ind(values1, values2, equal_var=False)

        # Mann-Whitney U test (robust primary test)
        u_stat, p_value_mwu = stats.mannwhitneyu(values1, values2, alternative='two-sided')
        
        # Cohen's d (standardized effect size)
        pooled_std = np.sqrt((np.var(values1) + np.var(values2)) / 2)
        cohens_d = effect_size / pooled_std if pooled_std > 0 else 0
        
        # Store result
        pairwise_results.append({
            'species': species,
            'group1': group1,
            'group2': group2,
            'pair': f"{group1} vs {group2}",
            'mean1': mean1,
            'mean2': mean2,
            'effect_size': effect_size,
            'cohens_d': cohens_d,
            't_statistic': t_stat,
            'p_value_ttest': p_value_t,
            'u_statistic': u_stat,
            'p_value_mwu': p_value_mwu,
            'n1': len(values1),
            'n2': len(values2)
        })

# Convert to DataFrame
results_df = pd.DataFrame(pairwise_results)
print(f"  Completed {len(results_df)} total comparisons ({len(abundance_clr.index)} species × {n_pairs} pairs)")

# ============================================================================
# 3. FDR CORRECTION
# ============================================================================
print(f"\n[3/6] Applying FDR correction...")

# FDR correction across ALL comparisons (Mann-Whitney as primary)
_, p_adjusted_mwu, _, _ = multipletests(results_df['p_value_mwu'], method='fdr_bh')
results_df['p_adjusted'] = p_adjusted_mwu
results_df['significant'] = p_adjusted_mwu < FDR_THRESHOLD

# Keep Welch t-test FDR for reference
_, p_adjusted_ttest, _, _ = multipletests(results_df['p_value_ttest'], method='fdr_bh')
results_df['p_adjusted_ttest'] = p_adjusted_ttest
results_df['significant_ttest'] = p_adjusted_ttest < FDR_THRESHOLD

n_significant = results_df['significant'].sum()
print(f"  Significant comparisons: {n_significant}/{len(results_df)} ({100*n_significant/len(results_df):.1f}%)")

# Count significant species per pair
print("\n  Significant species by pair:")
for pair in pairs:
    pair_str = f"{pair[0]} vs {pair[1]}"
    n_sig = results_df[(results_df['pair'] == pair_str) & (results_df['significant'])].shape[0]
    print(f"    {pair_str}: {n_sig}/{len(abundance_clr.index)} species")

# Save results
results_df.to_csv(f'{OUTPUT_DIR}/pairwise_comparisons_all.csv', index=False)
results_significant = results_df[results_df['significant']].sort_values('p_adjusted')
results_significant.to_csv(f'{OUTPUT_DIR}/pairwise_comparisons_significant.csv', index=False)
print(f"\n  ✓ Saved pairwise_comparisons_all.csv ({len(results_df)} comparisons)")
print(f"  ✓ Saved pairwise_comparisons_significant.csv ({len(results_significant)} comparisons)")

# ============================================================================
# 4. SUMMARY STATISTICS
# ============================================================================
print(f"\n[4/6] Generating summary statistics...")

# Species-level summary: count significant pairs per species
species_summary = []
for species in abundance_clr.index:
    species_data = results_df[results_df['species'] == species]
    n_sig_pairs = species_data['significant'].sum()
    
    # Get strongest association
    sig_data = species_data[species_data['significant']]
    if len(sig_data) > 0:
        strongest = sig_data.loc[sig_data['p_adjusted'].idxmin()]
        strongest_pair = strongest['pair']
        strongest_p = strongest['p_adjusted']
        strongest_effect = strongest['effect_size']
    else:
        strongest_pair = 'None'
        strongest_p = 1.0
        strongest_effect = 0.0
    
    species_summary.append({
        'species': species,
        'n_significant_pairs': n_sig_pairs,
        'n_total_pairs': n_pairs,
        'pct_significant': 100 * n_sig_pairs / n_pairs,
        'strongest_pair': strongest_pair,
        'strongest_p_adjusted': strongest_p,
        'strongest_effect_size': strongest_effect
    })

species_summary_df = pd.DataFrame(species_summary)
species_summary_df = species_summary_df.sort_values('n_significant_pairs', ascending=False)
species_summary_df.to_csv(f'{OUTPUT_DIR}/species_pairwise_summary.csv', index=False)
print(f"  ✓ Saved species_pairwise_summary.csv")

# Print top species
print("\n  Top 10 species by number of significant pairwise differences:")
for idx, row in species_summary_df.head(10).iterrows():
    print(f"    {row['species'][:50]:50s} {int(row['n_significant_pairs'])}/{n_pairs} pairs ({row['pct_significant']:.0f}%)")

# ============================================================================
# 5. VISUALIZATIONS
# ============================================================================
print(f"\n[5/6] Creating visualizations...")

# --- Visualization 1: Heatmap of significant comparisons (top species) ---
print("  Creating pairwise significance heatmap...")

top_species = species_summary_df.head(TOP_N_SPECIES)['species'].values

# Create matrix: rows=species, cols=pairs
heatmap_data = []
for species in top_species:
    row = []
    for pair in pairs:
        pair_str = f"{pair[0]} vs {pair[1]}"
        comparison = results_df[(results_df['species'] == species) & (results_df['pair'] == pair_str)]
        if len(comparison) > 0:
            p_adj = comparison['p_adjusted'].values[0]
            # Use -log10(p) for visualization
            row.append(-np.log10(p_adj) if comparison['significant'].values[0] else 0)
        else:
            row.append(0)
    heatmap_data.append(row)

heatmap_df = pd.DataFrame(
    heatmap_data,
    index=[s.split('|')[-1][:40] for s in top_species],
    columns=[f"{p[0][:4]}\nvs\n{p[1][:4]}" for p in pairs]
)

fig, ax = plt.subplots(figsize=(10, 12))
sns.heatmap(heatmap_df, cmap='YlOrRd', annot=False, fmt='.1f', 
            cbar_kws={'label': '-log10(FDR p-value)'}, ax=ax)
ax.set_xlabel('Pairwise Comparison', fontsize=11, fontweight='bold')
ax.set_ylabel('Species', fontsize=11, fontweight='bold')
ax.set_title(f'Top {TOP_N_SPECIES} Species: Pairwise Comparison Significance', 
             fontsize=12, fontweight='bold', pad=20)
plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/01_pairwise_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved 01_pairwise_heatmap.png")

# --- Visualization 2: Effect size heatmap (top species) ---
print("  Creating effect size heatmap...")

effect_data = []
for species in top_species:
    row = []
    for pair in pairs:
        pair_str = f"{pair[0]} vs {pair[1]}"
        comparison = results_df[(results_df['species'] == species) & (results_df['pair'] == pair_str)]
        if len(comparison) > 0 and comparison['significant'].values[0]:
            row.append(comparison['effect_size'].values[0])
        else:
            row.append(0)
    effect_data.append(row)

effect_df = pd.DataFrame(
    effect_data,
    index=[s.split('|')[-1][:40] for s in top_species],
    columns=[f"{p[0][:4]}\nvs\n{p[1][:4]}" for p in pairs]
)

fig, ax = plt.subplots(figsize=(10, 12))
sns.heatmap(effect_df, cmap='RdBu_r', center=0, annot=False, fmt='.1f',
            cbar_kws={'label': 'Effect Size (CLR difference)'}, ax=ax)
ax.set_xlabel('Pairwise Comparison', fontsize=11, fontweight='bold')
ax.set_ylabel('Species', fontsize=11, fontweight='bold')
ax.set_title(f'Top {TOP_N_SPECIES} Species: Pairwise Effect Sizes (FDR<0.05 only)', 
             fontsize=12, fontweight='bold', pad=20)
plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/02_effect_size_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved 02_effect_size_heatmap.png")

# --- Visualization 3: Summary barplot ---
print("  Creating summary barplot...")

# Count species by number of significant pairs
sig_counts = species_summary_df['n_significant_pairs'].value_counts().sort_index()

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: Distribution of significant pairs per species
ax = axes[0, 0]
sig_counts.plot(kind='bar', color='steelblue', ax=ax)
ax.set_xlabel('Number of Significant Pairs', fontsize=10, fontweight='bold')
ax.set_ylabel('Number of Species', fontsize=10, fontweight='bold')
ax.set_title('A) Distribution: How Many Pairs Differ per Species?', 
             fontsize=11, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Panel B: Pairwise comparison counts
ax = axes[0, 1]
pair_counts = results_df.groupby('pair')['significant'].sum().sort_values(ascending=False)
pair_counts.plot(kind='barh', color='coral', ax=ax)
ax.set_xlabel('Number of Significant Species', fontsize=10, fontweight='bold')
ax.set_ylabel('Pairwise Comparison', fontsize=10, fontweight='bold')
ax.set_title('B) Pairwise Comparison: Number of Significant Species', 
             fontsize=11, fontweight='bold')
ax.grid(axis='x', alpha=0.3)

# Panel C: Effect size distribution
ax = axes[1, 0]
sig_effects = results_df[results_df['significant']]['effect_size'].abs()
ax.hist(sig_effects, bins=50, color='forestgreen', edgecolor='black', alpha=0.7)
ax.axvline(sig_effects.median(), color='red', linestyle='--', linewidth=2, 
           label=f'Median: {sig_effects.median():.2f}')
ax.set_xlabel('Absolute Effect Size (CLR units)', fontsize=10, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=10, fontweight='bold')
ax.set_title('C) Distribution of Effect Sizes (Significant Only)', 
             fontsize=11, fontweight='bold')
ax.legend()
ax.grid(axis='y', alpha=0.3)

# Panel D: P-value distribution
ax = axes[1, 1]
ax.hist(-np.log10(results_df[results_df['significant']]['p_adjusted']), 
        bins=50, color='mediumpurple', edgecolor='black', alpha=0.7)
ax.axvline(-np.log10(FDR_THRESHOLD), color='red', linestyle='--', linewidth=2,
           label=f'FDR threshold (0.05)')
ax.set_xlabel('-log10(FDR p-value)', fontsize=10, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=10, fontweight='bold')
ax.set_title('D) Distribution of Adjusted P-values (Significant Only)', 
             fontsize=11, fontweight='bold')
ax.legend()
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/03_summary_statistics.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved 03_summary_statistics.png")

# --- Visualization 4: Network-style visualization ---
print("  Creating pairwise network diagram...")

fig, ax = plt.subplots(figsize=(12, 10))

# Position groups in a circle
n_groups = len(group_names)
angles = np.linspace(0, 2*np.pi, n_groups, endpoint=False)
positions = {group: (np.cos(angle), np.sin(angle)) 
             for group, angle in zip(group_names, angles)}

# Draw nodes (groups)
for group, (x, y) in positions.items():
    n_samples = (metadata['study_code'] == group).sum()
    ax.scatter(x, y, s=2000, c='lightblue', edgecolor='navy', linewidth=3, zorder=3)
    ax.text(x, y, f"{group}\n(n={n_samples})", ha='center', va='center', 
            fontsize=11, fontweight='bold', zorder=4)

# Draw edges (connections) - thickness by number of significant species
for pair in pairs:
    pair_str = f"{pair[0]} vs {pair[1]}"
    n_sig = results_df[(results_df['pair'] == pair_str) & (results_df['significant'])].shape[0]
    
    x1, y1 = positions[pair[0]]
    x2, y2 = positions[pair[1]]
    
    # Line thickness proportional to number of significant species
    thickness = 0.5 + (n_sig / len(abundance_clr.index)) * 15
    color = 'red' if n_sig > 0.9 * len(abundance_clr.index) else 'orange' if n_sig > 0.5 * len(abundance_clr.index) else 'gray'
    
    ax.plot([x1, x2], [y1, y2], linewidth=thickness, color=color, alpha=0.6, zorder=1)
    
    # Add label at midpoint
    mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
    ax.text(mid_x, mid_y, f"{n_sig}", ha='center', va='center',
            fontsize=9, bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8), zorder=2)

ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title('Pairwise Comparisons Network\n(Numbers = significant species; thickness ∝ significance)', 
             fontsize=13, fontweight='bold', pad=20)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='red', lw=4, label=f'>90% species differ'),
    Line2D([0], [0], color='orange', lw=4, label='50-90% species differ'),
    Line2D([0], [0], color='gray', lw=4, label='<50% species differ')
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/04_pairwise_network.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved 04_pairwise_network.png")

# ============================================================================
# 6. GENERATE SUMMARY STATEMENTS
# ============================================================================
print(f"\n[6/6] Generating pairwise comparison statements...")

statements = []
statements.append("=" * 80)
statements.append("PAIRWISE COMPARISON STATEMENTS")
statements.append("=" * 80)
statements.append("")
statements.append(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
statements.append(f"Species Analyzed: {len(abundance_clr.index)}")
statements.append(f"Pairwise Comparisons: {n_pairs} per species ({n_pairs * len(abundance_clr.index)} total)")
statements.append(f"Sample-size normalization: {normalization_note}")
statements.append("Primary test: Mann-Whitney U (Welch t-test retained as reference)")
statements.append(f"FDR Threshold: {FDR_THRESHOLD}")
statements.append(f"Significant Comparisons: {n_significant}/{len(results_df)} ({100*n_significant/len(results_df):.1f}%)")
statements.append("")

# Overall summary
statements.append("=" * 80)
statements.append("OVERALL SUMMARY")
statements.append("=" * 80)
statements.append("")

for pair in pairs:
    pair_str = f"{pair[0]} vs {pair[1]}"
    pair_data = results_df[results_df['pair'] == pair_str]
    n_sig = pair_data['significant'].sum()
    pct_sig = 100 * n_sig / len(pair_data)
    
    statements.append(f"{pair_str}:")
    statements.append(f"  - {n_sig}/{len(pair_data)} species significantly different ({pct_sig:.1f}%)")
    
    if n_sig > 0:
        sig_data = pair_data[pair_data['significant']]
        mean_effect = sig_data['effect_size'].abs().mean()
        median_effect = sig_data['effect_size'].abs().median()
        statements.append(f"  - Mean effect size: {mean_effect:.2f} CLR units")
        statements.append(f"  - Median effect size: {median_effect:.2f} CLR units")
    statements.append("")

# Top species statements
statements.append("=" * 80)
statements.append(f"TOP {TOP_N_SPECIES} SPECIES WITH MOST PAIRWISE DIFFERENCES")
statements.append("=" * 80)
statements.append("")

for idx, row in species_summary_df.head(TOP_N_SPECIES).iterrows():
    species_name = row['species'].split('|')[-1]
    n_sig = int(row['n_significant_pairs'])
    
    statements.append(f"\n{idx+1}. {species_name}")
    statements.append(f"   {'-' * 60}")
    statements.append(f"   Significant in {n_sig}/{n_pairs} pairwise comparisons ({row['pct_significant']:.0f}%)")
    
    # Get details for all pairs
    species_data = results_df[results_df['species'] == row['species']].sort_values('p_adjusted')
    
    for _, comp in species_data.iterrows():
        if comp['significant']:
            direction = "higher" if comp['effect_size'] > 0 else "lower"
            statements.append(f"   • {comp['pair']}: {direction} in {comp['group1']} "
                            f"(effect: {comp['effect_size']:+.2f} CLR, "
                            f"Cohen's d: {comp['cohens_d']:+.2f}, "
                            f"p={comp['p_adjusted']:.2e})")

statements.append("\n" + "=" * 80)
statements.append("END OF REPORT")
statements.append("=" * 80)

# Save statements
statement_text = '\n'.join(statements)
report_path = f'{OUTPUT_DIR}/pairwise_comparison_statements.txt'
with open(report_path, 'w') as f:
    f.write(statement_text)
print(f"  ✓ Saved pairwise_comparison_statements.txt")
snapshot_report(report_path)

# Display summary
print("\n" + "="*80)
print("PAIRWISE ANALYSIS COMPLETE")
print("="*80)
print(f"\nGenerated outputs in: {OUTPUT_DIR}/")
print(f"  • CSV files: 3")
print(f"  • Visualizations: 4 PNG files")
print(f"  • Statement report: 1 TXT file")
print(f"\nKey Finding: {n_significant}/{len(results_df)} pairwise comparisons significant (FDR < {FDR_THRESHOLD})")
print("\nMost discriminative pair:")
most_discriminative = results_df.groupby('pair')['significant'].sum().sort_values(ascending=False).iloc[0]
most_discriminative_pair = results_df.groupby('pair')['significant'].sum().sort_values(ascending=False).index[0]
print(f"  {most_discriminative_pair}: {int(most_discriminative)}/{len(abundance_clr.index)} species differ")
print("\n" + "="*80)
