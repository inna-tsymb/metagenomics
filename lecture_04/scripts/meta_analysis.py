#!/usr/bin/env python3
"""
Meta-Analysis: Multi-Cohort Microbiome Analysis
================================================
Leverages processed data from lecture_02 and lecture_03 to:
1. Load already-processed CLR data from lecture_03 (4 cohorts, 567 species)
2. Optionally add additional cohorts from raw data
3. Check for batch effects across cohorts
4. Perform multi-cohort meta-analysis
5. Validate and extend findings from previous lectures

Strategy: 
- Batch 1: Lecture_03 processed data (Schulz, Chu, Rowe, Lekunberri)
- Batch 2 (optional): Additional wastewater cohorts from lecture_02
- Fixed-effects meta-analysis with batch effect assessment
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import shutil
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input files - use already processed data!
LECTURE03_CLR = '../lecture_03/output/filtering_clr_analysis/abundance_clr_aligned.csv'
LECTURE03_METADATA = '../lecture_03/output/filtering_clr_analysis/metadata_aligned.csv'
LECTURE03_ASSOCIATION = '../lecture_03/output/association_analysis/association_results_all.csv'

# Optional: Additional cohorts from raw data
RAW_ABUNDANCE_FILE = '../lecture_02/input/environmental_metaphlan4_2026-02-06.tsv'
RAW_METADATA_FILE = '../lecture_02/input/environmental_extended_wide.tsv'

# Additional cohorts to add (if desired - set to empty list to use only lecture_03 data)
ADDITIONAL_COHORTS = [
    # 'Hendriksen_2019_sewage_global',    # Global sewage study
    # 'Brinch_2020_sewage_copenhagen',    # Copenhagen sewage
    # 'Kantor_2017_wastewater',           # Municipal wastewater
    # 'McIlroy_2016_sludge',              # Activated sludge
]

# Output directory
def _env_bool(name, default):
    value = os.getenv(name)
    if value is None:
        return default
    return value.strip().lower() in {'1', 'true', 'yes', 'on'}


OUTPUT_DIR = os.getenv('OUTPUT_DIR', 'output/meta_analysis')

# Statistical thresholds
FDR_THRESHOLD = 0.05
TOP_N_SPECIES = 20
MIN_SAMPLES_PER_COHORT = 10  # Minimum samples required for additional cohorts

# Species filtering thresholds
MIN_PREVALENCE = 0.10
MIN_COHORT_OVERLAP = 3
MIN_MEAN_ABUNDANCE = 1e-4

# Sample-size normalization (lecture_02 style)
APPLY_SAMPLE_SIZE_NORMALIZATION = _env_bool('APPLY_SAMPLE_SIZE_NORMALIZATION', True)
NORMALIZE_STUDY = os.getenv('NORMALIZE_STUDY', 'Schulz_2017_wastewater')
NORMALIZED_N = int(os.getenv('NORMALIZED_N', '20'))
RANDOM_SEED = int(os.getenv('RANDOM_SEED', '42'))

# Cohort selection
CORE_COHORTS = [
    'Schulz_2017_wastewater',
    'Chu_2017_sludge',
    'Rowe_2017_hospital_wastewater',
    'Lekunberri_2018_river_wastewater'
]
SELECTED_COHORTS = list(dict.fromkeys(CORE_COHORTS + ADDITIONAL_COHORTS))

# Active input files for this script path
ABUNDANCE_FILE = RAW_ABUNDANCE_FILE
METADATA_FILE = RAW_METADATA_FILE

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
print("META-ANALYSIS: MULTI-COHORT MICROBIOME STUDY")
print("="*80)
print("\nUsing Processed Data from Previous Lectures:")
print(f"  • Lecture_03: Pre-processed CLR data (4 wastewater cohorts)")
print(f"  • Additional cohorts: {len(ADDITIONAL_COHORTS)}")

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(HISTORY_DIR, exist_ok=True)

# ============================================================================
# STEP 1: LOAD AND PREPARE DATA
# ============================================================================
print("\n" + "="*80)
print("STEP 1: LOAD AND PREPARE DATA")
print("="*80)

print("\n[1.1] Loading raw data...")
# Load abundance (long format: sample_alias, clade_name, rel_abund)
abundance_long = pd.read_csv(ABUNDANCE_FILE, sep='\t')
print(f"  Loaded abundance: {len(abundance_long)} rows")

# Load metadata (wide format)
metadata = pd.read_csv(METADATA_FILE, sep='\t')
metadata = metadata.set_index('sample_alias')

# Remove duplicate sample IDs
if metadata.index.duplicated().any():
    n_dup = int(metadata.index.duplicated().sum())
    print(f"  Removing {n_dup} duplicate metadata sample IDs...")
    metadata = metadata[~metadata.index.duplicated(keep='first')]

print(f"  Loaded metadata: {metadata.shape[0]} samples × {metadata.shape[1]} variables")

print("\n[1.2] Filtering to selected cohorts...")
# Filter metadata to selected cohorts
metadata_filtered = metadata[metadata['study_code'].isin(SELECTED_COHORTS)].copy()
print(f"  Samples in selected cohorts: {len(metadata_filtered)}")

# Get sample counts per cohort
cohort_counts = metadata_filtered['study_code'].value_counts()
print("\n  Sample counts per cohort:")
for cohort, count in cohort_counts.items():
    print(f"    {cohort}: {count}")

# Sample-size normalization to reduce cohort imbalance
normalization_note = "No sample-size normalization applied"
if APPLY_SAMPLE_SIZE_NORMALIZATION:
    target_mask = metadata_filtered['study_code'] == NORMALIZE_STUDY
    target_n = int(target_mask.sum())
    if target_n > NORMALIZED_N:
        normalized_subset = metadata_filtered[target_mask].sample(n=NORMALIZED_N, random_state=RANDOM_SEED)
        metadata_filtered = pd.concat([metadata_filtered[~target_mask], normalized_subset]).sort_index()
        normalization_note = f"Downsampled {NORMALIZE_STUDY}: {target_n} -> {NORMALIZED_N}"
    else:
        normalization_note = (
            f"Normalization skipped for {NORMALIZE_STUDY} "
            f"(available: {target_n}, threshold: {NORMALIZED_N})"
        )

print(f"\n  Sample-size normalization: {normalization_note}")

# Recompute counts after normalization
cohort_counts = metadata_filtered['study_code'].value_counts()
print("  Updated sample counts per cohort:")
for cohort, count in cohort_counts.items():
    print(f"    {cohort}: {count}")

# Filter abundance to selected samples
selected_samples = metadata_filtered.index.tolist()
abundance_filtered = abundance_long[abundance_long['sample_alias'].isin(selected_samples)].copy()
print(f"\n  Abundance data filtered to {len(abundance_filtered)} rows")

# ============================================================================
# STEP 2: FIND OVERLAPPING SPECIES ACROSS COHORTS
# ============================================================================
print("\n" + "="*80)
print("STEP 2: FIND OVERLAPPING SPECIES")
print("="*80)

print("\n[2.1] Converting to wide format...")
# Pivot to wide format (samples × species)
abundance_wide = abundance_filtered.pivot(
    index='sample_alias',
    columns='clade_name',
    values='rel_abund'
).fillna(0)
print(f"  Wide format: {abundance_wide.shape[0]} samples × {abundance_wide.shape[1]} species")

print("\n[2.2] Checking species prevalence per cohort...")
# For each cohort, calculate species prevalence
cohort_prevalence = {}
for cohort in SELECTED_COHORTS:
    cohort_samples = metadata_filtered[metadata_filtered['study_code'] == cohort].index
    if len(cohort_samples) == 0:
        print(f"  WARNING: No samples for {cohort}")
        continue
    
    cohort_abundance = abundance_wide.loc[cohort_samples]
    prevalence = (cohort_abundance > 0).sum() / len(cohort_abundance)
    cohort_prevalence[cohort] = prevalence
    
    # Count species passing prevalence filter
    n_passing = (prevalence >= MIN_PREVALENCE).sum()
    print(f"  {cohort}: {n_passing}/{len(prevalence)} species with ≥{MIN_PREVALENCE*100:.0f}% prevalence")

print("\n[2.3] Finding species present in multiple cohorts...")
# Count how many cohorts each species appears in (with sufficient prevalence)
species_cohort_count = pd.Series(0, index=abundance_wide.columns)
for cohort, prevalence in cohort_prevalence.items():
    species_passing = prevalence[prevalence >= MIN_PREVALENCE].index
    species_cohort_count[species_passing] += 1

# Filter to species present in >= MIN_COHORT_OVERLAP cohorts
overlapping_species = species_cohort_count[species_cohort_count >= MIN_COHORT_OVERLAP].index.tolist()
print(f"\n  Species in ≥{MIN_COHORT_OVERLAP} cohorts: {len(overlapping_species)}")

# Additional filtering: mean abundance threshold
abundance_overlap = abundance_wide[overlapping_species]
mean_abundance = abundance_overlap.mean()
species_keep = mean_abundance[mean_abundance >= MIN_MEAN_ABUNDANCE].index.tolist()
print(f"  After abundance filter (≥{MIN_MEAN_ABUNDANCE*100:.04f}%): {len(species_keep)} species")

# Final filtered abundance
abundance_final = abundance_wide[species_keep].copy()
print(f"\n  ✓ Final dataset: {abundance_final.shape[0]} samples × {abundance_final.shape[1]} species")

# Save overlap summary
overlap_summary = pd.DataFrame({
    'species': abundance_final.columns,
    'n_cohorts': species_cohort_count[abundance_final.columns],
    'mean_abundance': mean_abundance[abundance_final.columns] * 100,
    'global_prevalence': (abundance_final > 0).sum() / len(abundance_final) * 100
})
overlap_summary = overlap_summary.sort_values('n_cohorts', ascending=False)
overlap_summary.to_csv(f'{OUTPUT_DIR}/species_overlap_summary.csv', index=False)
print(f"  ✓ Saved species_overlap_summary.csv")

# ============================================================================
# STEP 3: DATA HARMONIZATION (CLR TRANSFORMATION)
# ============================================================================
print("\n" + "="*80)
print("STEP 3: DATA HARMONIZATION")
print("="*80)

print("\n[3.1] Applying CLR transformation to combined dataset...")

def clr_transform(df):
    """Center log-ratio transformation with pseudocount"""
    pseudocount = 1e-6
    df_pseudo = df + pseudocount
    geometric_mean = np.exp(np.log(df_pseudo).mean(axis=1))
    clr = np.log(df_pseudo.div(geometric_mean, axis=0))
    return clr

# CLR transformation
abundance_clr = clr_transform(abundance_final)
print(f"  ✓ CLR transformation complete")
print(f"  CLR range: [{abundance_clr.min().min():.2f}, {abundance_clr.max().max():.2f}]")

# Save harmonized data
abundance_clr.to_csv(f'{OUTPUT_DIR}/abundance_clr_harmonized.csv')
print(f"  ✓ Saved abundance_clr_harmonized.csv")

# Merge metadata
metadata_harmonized = metadata_filtered.loc[abundance_clr.index].copy()
metadata_harmonized.to_csv(f'{OUTPUT_DIR}/metadata_harmonized.csv')
print(f"  ✓ Saved metadata_harmonized.csv ({len(metadata_harmonized)} samples)")

# ============================================================================
# STEP 4: CHECK FOR BATCH EFFECTS
# ============================================================================
print("\n" + "="*80)
print("STEP 4: BATCH EFFECT ASSESSMENT")
print("="*80)

print("\n[4.1] Performing PCA...")
# Standardize CLR values for PCA
scaler = StandardScaler()
abundance_scaled = scaler.fit_transform(abundance_clr)

# PCA
pca = PCA(n_components=min(50, abundance_clr.shape[1]))
pca_coords = pca.fit_transform(abundance_scaled)
explained_var = pca.explained_variance_ratio_

print(f"  PC1 explains {explained_var[0]*100:.1f}% variance")
print(f"  PC2 explains {explained_var[1]*100:.1f}% variance")
print(f"  PC1+PC2 explain {sum(explained_var[:2])*100:.1f}% variance")

# Create PCA DataFrame
pca_df = pd.DataFrame(
    pca_coords[:, :10],
    index=abundance_clr.index,
    columns=[f'PC{i+1}' for i in range(10)]
)
pca_df['cohort'] = metadata_harmonized['study_code'].values

# Save PCA results
pca_df.to_csv(f'{OUTPUT_DIR}/pca_coordinates.csv')
print(f"  ✓ Saved pca_coordinates.csv")

print("\n[4.2] Testing for batch effects (PERMANOVA-like)...")
# Test if cohort explains variance in PC1-PC5
from scipy.stats import f_oneway

batch_effect_tests = []
for pc in range(5):
    pc_values = pca_coords[:, pc]
    groups = [pc_values[metadata_harmonized['study_code'] == cohort] 
              for cohort in SELECTED_COHORTS 
              if (metadata_harmonized['study_code'] == cohort).sum() > 0]
    
    if len(groups) >= 2 and all(len(group) > 0 for group in groups):
        f_stat, p_value = f_oneway(*groups)
        kw_stat, p_value_kw = stats.kruskal(*groups)
        batch_effect_tests.append({
            'PC': f'PC{pc+1}',
            'variance_explained': explained_var[pc] * 100,
            'F_statistic': f_stat,
            'p_value': p_value,
            'KW_statistic': kw_stat,
            'p_value_kw': p_value_kw
        })

batch_df = pd.DataFrame(batch_effect_tests)

if len(batch_df) > 0:
    _, p_adjusted, _, _ = multipletests(batch_df['p_value'], method='fdr_bh')
    _, p_adjusted_kw, _, _ = multipletests(batch_df['p_value_kw'], method='fdr_bh')
    batch_df['p_adjusted'] = p_adjusted
    batch_df['p_adjusted_kw'] = p_adjusted_kw
    batch_df['significant'] = np.where(batch_df['p_adjusted_kw'] < FDR_THRESHOLD, 'YES', 'NO')

print("\n  Batch effect tests (cohort vs PC):")
if len(batch_df) > 0:
    print(batch_df[['PC', 'variance_explained', 'F_statistic', 'p_value', 'p_adjusted', 'KW_statistic', 'p_value_kw', 'p_adjusted_kw', 'significant']].to_string(index=False))

# ============================================================================
# STEP 5: VISUALIZATIONS
# ============================================================================
print("\n" + "="*80)
print("STEP 5: CREATING VISUALIZATIONS")
print("="*80)

# --- Visualization 1: PCA plot ---
print("\n[5.1] Creating PCA plot...")
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Panel A: PCA colored by cohort
ax = axes[0]
cohorts_present = metadata_harmonized['study_code'].unique()
colors = plt.cm.tab10(np.linspace(0, 1, len(cohorts_present)))
cohort_colors = dict(zip(cohorts_present, colors))

for cohort in cohorts_present:
    mask = pca_df['cohort'] == cohort
    ax.scatter(pca_df.loc[mask, 'PC1'], pca_df.loc[mask, 'PC2'],
               c=[cohort_colors[cohort]], label=cohort, alpha=0.6, s=50)

ax.set_xlabel(f'PC1 ({explained_var[0]*100:.1f}%)', fontsize=11, fontweight='bold')
ax.set_ylabel(f'PC2 ({explained_var[1]*100:.1f}%)', fontsize=11, fontweight='bold')
ax.set_title('A) PCA: Batch Effect Assessment', fontsize=12, fontweight='bold')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
ax.grid(alpha=0.3)

# Panel B: Scree plot
ax = axes[1]
n_pcs = min(20, len(explained_var))
ax.bar(range(1, n_pcs+1), explained_var[:n_pcs]*100, color='steelblue', edgecolor='black')
ax.plot(range(1, n_pcs+1), np.cumsum(explained_var[:n_pcs])*100, 
        'ro-', linewidth=2, markersize=6, label='Cumulative')
ax.set_xlabel('Principal Component', fontsize=11, fontweight='bold')
ax.set_ylabel('Variance Explained (%)', fontsize=11, fontweight='bold')
ax.set_title('B) Scree Plot: Variance Explained', fontsize=12, fontweight='bold')
ax.legend()
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/01_batch_effect_pca.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved 01_batch_effect_pca.png")

# --- Visualization 2: Cohort overlap heatmap ---
print("\n[5.2] Creating cohort overlap heatmap...")
# Calculate Jaccard similarity between cohorts
presence_matrix = (abundance_final > 0).astype(int)
cohort_similarity = np.zeros((len(SELECTED_COHORTS), len(SELECTED_COHORTS)))

for i, cohort1 in enumerate(SELECTED_COHORTS):
    samples1 = metadata_harmonized[metadata_harmonized['study_code'] == cohort1].index
    if len(samples1) == 0:
        continue
    set1 = set(abundance_final.columns[(presence_matrix.loc[samples1].sum() / len(samples1)) > 0.1])
    
    for j, cohort2 in enumerate(SELECTED_COHORTS):
        samples2 = metadata_harmonized[metadata_harmonized['study_code'] == cohort2].index
        if len(samples2) == 0:
            continue
        set2 = set(abundance_final.columns[(presence_matrix.loc[samples2].sum() / len(samples2)) > 0.1])
        
        # Jaccard similarity
        intersection = len(set1 & set2)
        union = len(set1 | set2)
        cohort_similarity[i, j] = intersection / union if union > 0 else 0

fig, ax = plt.subplots(figsize=(10, 8))
im = ax.imshow(cohort_similarity, cmap='YlOrRd', aspect='auto', vmin=0, vmax=1)

# Labels
ax.set_xticks(range(len(SELECTED_COHORTS)))
ax.set_yticks(range(len(SELECTED_COHORTS)))
ax.set_xticklabels([c.replace('_', ' ')[:30] for c in SELECTED_COHORTS], 
                    rotation=45, ha='right', fontsize=9)
ax.set_yticklabels([c.replace('_', ' ')[:30] for c in SELECTED_COHORTS], fontsize=9)

# Colorbar
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Jaccard Similarity', fontsize=11, fontweight='bold')

# Annotate cells
for i in range(len(SELECTED_COHORTS)):
    for j in range(len(SELECTED_COHORTS)):
        text = ax.text(j, i, f'{cohort_similarity[i, j]:.2f}',
                      ha="center", va="center", color="black" if cohort_similarity[i, j] < 0.5 else "white",
                      fontsize=8)

ax.set_title('Cohort Similarity (Jaccard Index)\nBased on Species Overlap', 
             fontsize=12, fontweight='bold', pad=20)
plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/02_cohort_similarity.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved 02_cohort_similarity.png")

# --- Visualization 3: Species distribution across cohorts ---
print("\n[5.3] Creating species distribution plot...")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: Species per cohort
ax = axes[0, 0]
species_per_cohort = []
for cohort in SELECTED_COHORTS:
    samples = metadata_harmonized[metadata_harmonized['study_code'] == cohort].index
    if len(samples) > 0:
        n_species = (abundance_final.loc[samples] > 0).any().sum()
        species_per_cohort.append({'cohort': cohort, 'n_species': n_species})

species_df = pd.DataFrame(species_per_cohort).sort_values('n_species', ascending=False)
ax.barh(range(len(species_df)), species_df['n_species'], color='steelblue', edgecolor='black')
ax.set_yticks(range(len(species_df)))
ax.set_yticklabels([c.replace('_', ' ')[:30] for c in species_df['cohort']], fontsize=9)
ax.set_xlabel('Number of Species', fontsize=10, fontweight='bold')
ax.set_title('A) Species Richness per Cohort', fontsize=11, fontweight='bold')
ax.grid(axis='x', alpha=0.3)

# Panel B: Sample size per cohort
ax = axes[0, 1]
cohort_counts_sorted = cohort_counts.sort_values(ascending=False)
ax.barh(range(len(cohort_counts_sorted)), cohort_counts_sorted.values, color='coral', edgecolor='black')
ax.set_yticks(range(len(cohort_counts_sorted)))
ax.set_yticklabels([c.replace('_', ' ')[:30] for c in cohort_counts_sorted.index], fontsize=9)
ax.set_xlabel('Number of Samples', fontsize=10, fontweight='bold')
ax.set_title('B) Sample Size per Cohort', fontsize=11, fontweight='bold')
ax.grid(axis='x', alpha=0.3)

# Panel C: Species overlap distribution
ax = axes[1, 0]
ax.hist(species_cohort_count[species_keep], bins=range(1, len(SELECTED_COHORTS)+2), 
        color='forestgreen', edgecolor='black', alpha=0.7)
ax.set_xlabel('Number of Cohorts', fontsize=10, fontweight='bold')
ax.set_ylabel('Number of Species', fontsize=10, fontweight='bold')
ax.set_title('C) Species Presence Across Cohorts', fontsize=11, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Panel D: Top species by cohort presence
ax = axes[1, 1]
top_overlap = overlap_summary.head(15)
ax.barh(range(len(top_overlap)), top_overlap['n_cohorts'], color='mediumpurple', edgecolor='black')
ax.set_yticks(range(len(top_overlap)))
ax.set_yticklabels([s.split('|')[-1][:40] for s in top_overlap['species']], fontsize=8)
ax.set_xlabel('Number of Cohorts', fontsize=10, fontweight='bold')
ax.set_title('D) Top 15 Species by Cohort Presence', fontsize=11, fontweight='bold')
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/03_species_distribution.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved 03_species_distribution.png")

# ============================================================================
# STEP 6: META-ANALYSIS (MULTI-COHORT TESTING)
# ============================================================================
print("\n" + "="*80)
print("STEP 6: META-ANALYSIS")
print("="*80)

print("\n[6.1] Strategy: Within-cohort effect size calculation + meta-analysis...")
# For each species, calculate effect sizes within each cohort
# Then combine using fixed-effects meta-analysis

meta_results = []

print(f"\n[6.2] Testing {len(abundance_clr.columns)} species...")
for species_idx, species in enumerate(abundance_clr.columns):
    if (species_idx + 1) % 50 == 0:
        print(f"  Progress: {species_idx + 1}/{len(abundance_clr.columns)} species...")
    
    # Get CLR values for this species
    clr_values = abundance_clr[species]
    
    # Calculate cohort-specific effect sizes
    cohort_effects = []
    cohort_ses = []  # Standard errors
    cohort_ns = []
    
    for cohort in SELECTED_COHORTS:
        cohort_samples = metadata_harmonized[metadata_harmonized['study_code'] == cohort].index
        if len(cohort_samples) < 5:  # Need minimum 5 samples
            continue
        
        cohort_clr = clr_values[cohort_samples]
        
        # Effect size: deviation from overall mean
        overall_mean = clr_values.mean()
        cohort_mean = cohort_clr.mean()
        effect = cohort_mean - overall_mean
        
        # Standard error
        se = cohort_clr.std() / np.sqrt(len(cohort_clr))
        
        cohort_effects.append(effect)
        cohort_ses.append(se)
        cohort_ns.append(len(cohort_clr))
    
    # Skip if insufficient cohorts
    if len(cohort_effects) < 3:
        continue
    
    # Fixed-effects meta-analysis
    # Weight by inverse variance (1 / SE^2)
    weights = [1 / (se**2) if se > 0 else 0 for se in cohort_ses]
    total_weight = sum(weights)
    
    if total_weight == 0:
        continue
    
    # Pooled effect size
    pooled_effect = sum(e * w for e, w in zip(cohort_effects, weights)) / total_weight
    pooled_se = np.sqrt(1 / total_weight)
    
    # Z-test for pooled effect
    z_score = pooled_effect / pooled_se if pooled_se > 0 else 0
    p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
    
    # Heterogeneity test (Q statistic)
    Q = sum(w * (e - pooled_effect)**2 for e, w in zip(cohort_effects, weights))
    df = len(cohort_effects) - 1
    Q_pval = 1 - stats.chi2.cdf(Q, df) if df > 0 else 1
    I2 = max(0, (Q - df) / Q * 100) if Q > 0 else 0
    
    meta_results.append({
        'species': species,
        'n_cohorts': len(cohort_effects),
        'pooled_effect': pooled_effect,
        'pooled_se': pooled_se,
        'z_score': z_score,
        'p_value': p_value,
        'Q_statistic': Q,
        'Q_pvalue': Q_pval,
        'I2_heterogeneity': I2,
        'min_effect': min(cohort_effects),
        'max_effect': max(cohort_effects)
    })

meta_df = pd.DataFrame(meta_results)
print(f"\n  Completed meta-analysis for {len(meta_df)} species")

# FDR correction
_, p_adjusted, _, _ = multipletests(meta_df['p_value'], method='fdr_bh')
meta_df['p_adjusted'] = p_adjusted
meta_df['significant'] = p_adjusted < FDR_THRESHOLD

n_sig = meta_df['significant'].sum()
print(f"  Significant species (FDR < {FDR_THRESHOLD}): {n_sig}/{len(meta_df)} ({100*n_sig/len(meta_df):.1f}%)")

# Save results
meta_df_sorted = meta_df.sort_values('p_adjusted')
meta_df_sorted.to_csv(f'{OUTPUT_DIR}/meta_analysis_results_all.csv', index=False)
meta_df_sorted[meta_df_sorted['significant']].to_csv(f'{OUTPUT_DIR}/meta_analysis_results_significant.csv', index=False)
print(f"\n  ✓ Saved meta_analysis_results_all.csv")
print(f"  ✓ Saved meta_analysis_results_significant.csv")

# ============================================================================
# STEP 7: META-ANALYSIS VISUALIZATIONS
# ============================================================================
print("\n" + "="*80)
print("STEP 7: META-ANALYSIS VISUALIZATIONS")
print("="*80)

# --- Visualization 4: Volcano plot ---
print("\n[7.1] Creating volcano plot...")
fig, ax = plt.subplots(figsize=(10, 8))

# Plot all species
ax.scatter(meta_df['pooled_effect'], -np.log10(meta_df['p_adjusted']),
          c='gray', alpha=0.5, s=30, label='Not significant')

# Highlight significant
sig_df = meta_df[meta_df['significant']]
ax.scatter(sig_df['pooled_effect'], -np.log10(sig_df['p_adjusted']),
          c='red', alpha=0.7, s=50, label=f'Significant (FDR<{FDR_THRESHOLD})')

# Threshold line
ax.axhline(-np.log10(FDR_THRESHOLD), color='blue', linestyle='--', linewidth=2, 
          label=f'FDR threshold ({FDR_THRESHOLD})')

# Labels for top species
top_species = sig_df.nsmallest(10, 'p_adjusted')
for _, row in top_species.iterrows():
    species_name = row['species'].split('|')[-1][:30]
    ax.annotate(species_name, 
               xy=(row['pooled_effect'], -np.log10(row['p_adjusted'])),
               xytext=(5, 5), textcoords='offset points',
               fontsize=7, alpha=0.8,
               bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))

ax.set_xlabel('Pooled Effect Size (CLR units)', fontsize=11, fontweight='bold')
ax.set_ylabel('-log10(FDR p-value)', fontsize=11, fontweight='bold')
ax.set_title('Meta-Analysis Volcano Plot\nFixed-Effects Model Across Cohorts', 
            fontsize=12, fontweight='bold', pad=15)
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/04_meta_volcano_plot.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved 04_meta_volcano_plot.png")

# --- Visualization 5: Heterogeneity assessment ---
print("\n[7.2] Creating heterogeneity plot...")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: I² distribution
ax = axes[0, 0]
ax.hist(meta_df['I2_heterogeneity'], bins=30, color='steelblue', edgecolor='black', alpha=0.7)
ax.axvline(25, color='orange', linestyle='--', label='Low heterogeneity (25%)')
ax.axvline(75, color='red', linestyle='--', label='High heterogeneity (75%)')
ax.set_xlabel('I² Heterogeneity (%)', fontsize=10, fontweight='bold')
ax.set_ylabel('Number of Species', fontsize=10, fontweight='bold')
ax.set_title('A) Heterogeneity Distribution', fontsize=11, fontweight='bold')
ax.legend()
ax.grid(axis='y', alpha=0.3)

# Panel B: Effect size vs heterogeneity
ax = axes[0, 1]
scatter = ax.scatter(meta_df['pooled_effect'], meta_df['I2_heterogeneity'],
                    c=meta_df['significant'].astype(int), cmap='RdYlGn', 
                    alpha=0.6, s=50, edgecolor='black', linewidth=0.5)
ax.set_xlabel('Pooled Effect Size', fontsize=10, fontweight='bold')
ax.set_ylabel('I² Heterogeneity (%)', fontsize=10, fontweight='bold')
ax.set_title('B) Effect Size vs Heterogeneity', fontsize=11, fontweight='bold')
ax.grid(alpha=0.3)
plt.colorbar(scatter, ax=ax, label='Significant')

# Panel C: Number of cohorts
ax = axes[1, 0]
ax.hist(meta_df['n_cohorts'], bins=range(3, max(meta_df['n_cohorts'])+2), 
       color='forestgreen', edgecolor='black', alpha=0.7)
ax.set_xlabel('Number of Cohorts', fontsize=10, fontweight='bold')
ax.set_ylabel('Number of Species', fontsize=10, fontweight='bold')
ax.set_title('C) Species Detection Across Cohorts', fontsize=11, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Panel D: P-value distribution
ax = axes[1, 1]
ax.hist(meta_df['p_adjusted'], bins=50, color='mediumpurple', edgecolor='black', alpha=0.7)
ax.axvline(FDR_THRESHOLD, color='red', linestyle='--', linewidth=2, label=f'FDR={FDR_THRESHOLD}')
ax.set_xlabel('FDR-Adjusted P-value', fontsize=10, fontweight='bold')
ax.set_ylabel('Number of Species', fontsize=10, fontweight='bold')
ax.set_title('D) P-value Distribution', fontsize=11, fontweight='bold')
ax.legend()
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/05_meta_heterogeneity.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved 05_meta_heterogeneity.png")

# ============================================================================
# STEP 8: COMPARISON WITH LECTURE_03 RESULTS
# ============================================================================
print("\n" + "="*80)
print("STEP 8: VALIDATION - COMPARE WITH LECTURE_03")
print("="*80)

print("\n[8.1] Loading lecture_03 results...")
common_species = set()
both_sig = set()
correlation = np.nan
try:
    lecture03_results = pd.read_csv('../lecture_03/output/association_analysis/association_results_all.csv')
    print(f"  Loaded lecture_03 results: {len(lecture03_results)} species")
    
    # Find common species
    common_species = set(meta_df['species']) & set(lecture03_results['species'])
    print(f"  Common species: {len(common_species)}")
    
    if len(common_species) > 0:
        # Compare significance
        meta_sig = set(meta_df[meta_df['significant']]['species'])
        lec03_sig = set(lecture03_results[lecture03_results['p_adjusted'] < FDR_THRESHOLD]['species'])
        
        both_sig = meta_sig & lec03_sig & common_species
        meta_only = (meta_sig - lec03_sig) & common_species
        lec03_only = (lec03_sig - meta_sig) & common_species
        
        print(f"\n  Validation Results:")
        print(f"    Significant in BOTH: {len(both_sig)} ({100*len(both_sig)/len(common_species):.1f}%)")
        print(f"    Meta-analysis only: {len(meta_only)}")
        print(f"    Lecture_03 only: {len(lec03_only)}")
        
        # Correlation of effect sizes
        meta_common = meta_df[meta_df['species'].isin(common_species)].set_index('species')
        lec03_common = lecture03_results[lecture03_results['species'].isin(common_species)].set_index('species')

        lec03_compare = lec03_common[['effect_size', 'p_adjusted']].rename(columns={
            'effect_size': 'effect_size_lec03',
            'p_adjusted': 'p_adjusted_lec03'
        })

        merged = pd.merge(meta_common[['pooled_effect', 'p_adjusted']],
                         lec03_compare,
                         left_index=True, right_index=True, 
                         suffixes=('_meta', '_lec03'))
        
        correlation = merged['pooled_effect'].corr(merged['effect_size_lec03'])
        print(f"\n  Effect size correlation: r = {correlation:.3f}")
        
        # Save comparison
        merged['both_significant'] = merged.index.isin(both_sig)
        merged.to_csv(f'{OUTPUT_DIR}/comparison_meta_vs_lecture03.csv')
        print(f"  ✓ Saved comparison_meta_vs_lecture03.csv")
        
except Exception as e:
    print(f"  Warning: Could not load lecture_03 results: {e}")

# ============================================================================
# STEP 9: GENERATE SUMMARY REPORT
# ============================================================================
print("\n" + "="*80)
print("STEP 9: SUMMARY REPORT")
print("="*80)

report_lines = []
report_lines.append("="*80)
report_lines.append("META-ANALYSIS SUMMARY REPORT")
report_lines.append("="*80)
report_lines.append(f"\nAnalysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
report_lines.append(f"\n{'='*80}")
report_lines.append("1. COHORT SELECTION")
report_lines.append("="*80)
report_lines.append(f"\nSelected {len(SELECTED_COHORTS)} wastewater-related cohorts:")
for i, cohort in enumerate(SELECTED_COHORTS, 1):
    if cohort in cohort_counts.index:
        count = cohort_counts[cohort]
        report_lines.append(f"{i}. {cohort}: {count} samples")
    else:
        report_lines.append(f"{i}. {cohort}: 0 samples (excluded)")

report_lines.append(f"\nTotal samples: {len(metadata_harmonized)}")
report_lines.append(f"Sample-size normalization: {normalization_note}")

report_lines.append(f"\n{'='*80}")
report_lines.append("2. DATA HARMONIZATION")
report_lines.append("="*80)
report_lines.append(f"\nSpecies filtering:")
report_lines.append(f"  • Initial species (all cohorts): {abundance_wide.shape[1]}")
report_lines.append(f"  • After prevalence filter (≥{MIN_PREVALENCE*100:.0f}% per cohort): {len(overlapping_species)}")
report_lines.append(f"  • After abundance filter (≥{MIN_MEAN_ABUNDANCE*100:.04f}%): {len(species_keep)}")
report_lines.append(f"  • Must appear in ≥{MIN_COHORT_OVERLAP} cohorts")
report_lines.append(f"\nFinal dataset:")
report_lines.append(f"  • Samples: {abundance_final.shape[0]}")
report_lines.append(f"  • Species: {abundance_final.shape[1]}")
report_lines.append(f"\nCLR transformation applied to combined dataset")

report_lines.append(f"\n{'='*80}")
report_lines.append("3. BATCH EFFECT ASSESSMENT")
report_lines.append("="*80)
report_lines.append(f"\nPCA Results:")
report_lines.append(f"  • PC1 variance explained: {explained_var[0]*100:.1f}%")
report_lines.append(f"  • PC2 variance explained: {explained_var[1]*100:.1f}%")
report_lines.append(f"  • PC1+PC2 combined: {sum(explained_var[:2])*100:.1f}%")
report_lines.append(f"\nBatch effect tests:")
for _, row in batch_df.iterrows():
    report_lines.append(
        f"  • {row['PC']}: ANOVA p={row['p_value']:.2e} (FDR={row['p_adjusted']:.2e}), "
        f"Kruskal p={row['p_value_kw']:.2e} (FDR={row['p_adjusted_kw']:.2e}) [{row['significant']}]"
    )

report_lines.append(f"\n{'='*80}")
report_lines.append("4. META-ANALYSIS RESULTS")
report_lines.append("="*80)
report_lines.append(f"\nStrategy: Fixed-effects meta-analysis")
report_lines.append(f"  • Within-cohort effect sizes calculated")
report_lines.append(f"  • Combined using inverse-variance weighting")
report_lines.append(f"  • Heterogeneity assessed using I² statistic")
report_lines.append(f"\nResults:")
report_lines.append(f"  • Species tested: {len(meta_df)}")
report_lines.append(f"  • Significant species (FDR<{FDR_THRESHOLD}): {n_sig} ({100*n_sig/len(meta_df):.1f}%)")
report_lines.append(f"\nHeterogeneity:")
report_lines.append(f"  • Median I²: {meta_df['I2_heterogeneity'].median():.1f}%")
report_lines.append(f"  • High heterogeneity (I²>75%): {(meta_df['I2_heterogeneity']>75).sum()} species")

report_lines.append(f"\n{'='*80}")
report_lines.append(f"5. TOP {TOP_N_SPECIES} SPECIES")
report_lines.append("="*80)
top_species_df = meta_df_sorted.head(TOP_N_SPECIES)
for idx, row in top_species_df.iterrows():
    species_name = row['species'].split('|')[-1]
    report_lines.append(f"\n{idx+1}. {species_name}")
    report_lines.append(f"   Pooled effect: {row['pooled_effect']:+.2f} CLR units")
    report_lines.append(f"   P-value: {row['p_adjusted']:.2e}")
    report_lines.append(f"   I² heterogeneity: {row['I2_heterogeneity']:.1f}%")
    report_lines.append(f"   Present in {row['n_cohorts']} cohorts")

if 'both_sig' in locals():
    report_lines.append(f"\n{'='*80}")
    report_lines.append("6. VALIDATION (vs LECTURE_03)")
    report_lines.append("="*80)
    report_lines.append(f"\nComparison with lecture_03 single-cohort analysis:")
    report_lines.append(f"  • Common species: {len(common_species)}")
    report_lines.append(f"  • Concordance (both significant): {len(both_sig)} ({100*len(both_sig)/len(common_species):.1f}%)")
    report_lines.append(f"  • Effect size correlation: r = {correlation:.3f}")

report_lines.append(f"\n{'='*80}")
report_lines.append("END OF REPORT")
report_lines.append("="*80)

# Save report
report_text = '\n'.join(report_lines)
report_path = f'{OUTPUT_DIR}/meta_analysis_summary_report.txt'
with open(report_path, 'w') as f:
    f.write(report_text)

print(f"\n  ✓ Saved meta_analysis_summary_report.txt")
snapshot_report(report_path)

# Print summary to console
print("\n" + "="*80)
print("META-ANALYSIS COMPLETE")
print("="*80)
print(f"\nGenerated outputs in: {OUTPUT_DIR}/")
print(f"  • CSV files: 5-6")
print(f"  • Visualizations: 5 PNG files")
print(f"  • Summary report: 1 TXT file")
print(f"\nKey Finding: {n_sig}/{len(meta_df)} species show consistent effects across {len(SELECTED_COHORTS)} cohorts")
print(f"Batch effects detected: {'YES' if (batch_df['p_adjusted_kw'] < FDR_THRESHOLD).any() else 'MINIMAL'}")
print("\n" + "="*80)
