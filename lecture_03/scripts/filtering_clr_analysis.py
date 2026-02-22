"""
================================================================================
MICROBIOME DATA FILTERING & CLR TRANSFORMATION ANALYSIS
================================================================================
This script performs:
1. Species filtering (prevalence & abundance)
2. Zero handling with pseudocount
3. CLR (Centered Log-Ratio) transformation for compositionality
4. Metadata alignment
5. Wastewater sample filtering
================================================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
import os
import json
from datetime import datetime

# ==========================================
# 1. SETUP & CONFIGURATION
# ==========================================
base_dir = '/Users/user/Documents/metagenomics'
input_dir = os.path.join(base_dir, 'lecture_02', 'input')
output_dir = os.path.join(base_dir, 'lecture_03', 'output', 'filtering_clr_analysis')
output_dir = os.getenv('OUTPUT_DIR', output_dir)
os.makedirs(output_dir, exist_ok=True)

# Wastewater studies to include
WASTEWATER_STUDIES = [
    'Schulz_2017_wastewater', 
    'Rowe_2017_hospital_wastewater', 
    'Chu_2017_sludge',
    'Lekunberri_2018_river_wastewater'
]

# Filtering thresholds
PREVALENCE_THRESHOLD = 0.10  # Keep species present in at least 10% of samples
ABUNDANCE_THRESHOLD = 0.01   # Keep species with mean abundance >= 0.01%
PSEUDOCOUNT = 1e-6           # Small value added before log transformation
RANDOM_SEED = 42

# Sample-size normalization (lecture_02 style)
def _env_bool(name, default):
    value = os.getenv(name)
    if value is None:
        return default
    return value.strip().lower() in {'1', 'true', 'yes', 'on'}


APPLY_SAMPLE_SIZE_NORMALIZATION = _env_bool('APPLY_SAMPLE_SIZE_NORMALIZATION', True)
NORMALIZE_STUDY = os.getenv('NORMALIZE_STUDY', 'Schulz_2017_wastewater')
NORMALIZED_N = int(os.getenv('NORMALIZED_N', '20'))

print("=" * 80)
print("MICROBIOME DATA FILTERING & CLR TRANSFORMATION")
print("=" * 80)
print(f"Output directory: {output_dir}\n")


def apply_sample_size_normalization(abund_df, meta_df):
    """Downsample oversized study cohort to reduce imbalance (lecture_02 style)."""
    common_samples = abund_df.index.intersection(meta_df.index)
    abund_df = abund_df.loc[common_samples].copy()
    meta_df = meta_df.loc[common_samples].copy()

    if 'study_code' not in meta_df.columns:
        return abund_df, meta_df, "Normalization skipped: 'study_code' not found", None

    counts_before = meta_df['study_code'].value_counts().sort_index()
    normalization_note = "No sample-size normalization applied"

    if not APPLY_SAMPLE_SIZE_NORMALIZATION:
        comparison_df = pd.DataFrame({
            'study_code': counts_before.index,
            'n_before': counts_before.values,
            'n_after': counts_before.values
        })
        return abund_df, meta_df, normalization_note, comparison_df

    target_mask = meta_df['study_code'] == NORMALIZE_STUDY
    target_n = int(target_mask.sum())

    if target_n > NORMALIZED_N:
        keep_target = meta_df[target_mask].sample(n=NORMALIZED_N, random_state=RANDOM_SEED)
        meta_after = pd.concat([meta_df[~target_mask], keep_target]).sort_index()
        common_after = abund_df.index.intersection(meta_after.index)
        abund_after = abund_df.loc[common_after].copy()
        meta_after = meta_after.loc[common_after].copy()
        normalization_note = f"Downsampled {NORMALIZE_STUDY}: {target_n} -> {NORMALIZED_N}"
    else:
        meta_after = meta_df.copy()
        abund_after = abund_df.copy()
        normalization_note = (
            f"Normalization skipped for {NORMALIZE_STUDY} "
            f"(available: {target_n}, threshold: {NORMALIZED_N})"
        )

    counts_after = meta_after['study_code'].value_counts().sort_index()
    all_studies = sorted(set(counts_before.index).union(set(counts_after.index)))
    comparison_df = pd.DataFrame({
        'study_code': all_studies,
        'n_before': [int(counts_before.get(study, 0)) for study in all_studies],
        'n_after': [int(counts_after.get(study, 0)) for study in all_studies]
    })

    return abund_after, meta_after, normalization_note, comparison_df

# ==========================================
# 2. HELPER FUNCTIONS
# ==========================================

def transform_long_to_wide(df_long):
    """Convert long format metaphlan data to wide format (samples x species)."""
    print("Converting long format to wide format...")
    
    # Filter for species level only (s__)
    df_species = df_long[
        df_long['clade_name'].str.contains('s__', na=False) & 
        ~df_long['clade_name'].str.contains('t__', na=False)
    ].copy()

    n_before = len(df_species)
    df_species = df_species.drop_duplicates()
    n_removed = n_before - len(df_species)
    if n_removed > 0:
        print(f"  Removed {n_removed} duplicate long-format abundance rows")
    
    # Pivot table (samples as rows, species as columns)
    abund_wide = df_species.pivot_table(
        index='sample_alias', 
        columns='clade_name', 
        values='rel_abund', 
        aggfunc='mean'
    )
    
    # Fill NaN with 0
    abund_wide = abund_wide.fillna(0)
    
    print(f"  ✓ Converted to wide format: {abund_wide.shape[0]} samples × {abund_wide.shape[1]} species")
    
    return abund_wide


def calculate_prevalence_abundance(abund_df):
    """Calculate prevalence and mean abundance for each species."""
    n_samples = abund_df.shape[0]
    
    # Prevalence: proportion of samples where species is detected (>0)
    prevalence = (abund_df > 0).sum() / n_samples
    
    # Mean abundance
    mean_abundance = abund_df.mean()
    
    stats_df = pd.DataFrame({
        'species': abund_df.columns,
        'prevalence': prevalence.values,
        'mean_abundance': mean_abundance.values,
        'median_abundance': abund_df.median().values,
        'max_abundance': abund_df.max().values,
        'n_samples_detected': (abund_df > 0).sum().values
    }).set_index('species')
    
    return stats_df


def filter_species(abund_df, prevalence_thresh, abundance_thresh):
    """Filter species based on prevalence and abundance thresholds."""
    print(f"\nFiltering species (prevalence ≥ {prevalence_thresh:.1%}, mean abundance ≥ {abundance_thresh}%)...")
    
    # Calculate statistics
    stats = calculate_prevalence_abundance(abund_df)
    
    # Apply filters
    mask = (stats['prevalence'] >= prevalence_thresh) | (stats['mean_abundance'] >= abundance_thresh)
    filtered_species = stats[mask].index.tolist()
    removed_species = stats[~mask].index.tolist()
    
    # Filter abundance table
    abund_filtered = abund_df[filtered_species].copy()
    
    print(f"  ✓ Removed {len(removed_species)} species, kept {len(filtered_species)} species")
    print(f"    Before: {abund_df.shape[1]} species")
    print(f"    After:  {abund_filtered.shape[1]} species (retained: {len(filtered_species)/abund_df.shape[1]*100:.1f}%)")
    
    return abund_filtered, filtered_species, removed_species, stats


def add_pseudocount_and_clr(abund_df, pseudocount=1e-6):
    """
    Apply pseudocount and CLR transformation.
    
    CLR (Centered Log-Ratio) is defined as:
    clr(x_i) = log(x_i / geom_mean(x))
    
    This handles compositionality in microbiome data.
    """
    print(f"\nApplying pseudocount ({pseudocount}) and CLR transformation...")
    
    # Add pseudocount
    abund_pseudo = abund_df + pseudocount
    
    # Calculate geometric mean for each sample
    # geom_mean = (prod(x))^(1/n) = exp(mean(log(x)))
    log_abund = np.log(abund_pseudo)
    geom_mean = np.exp(log_abund.mean(axis=1))
    
    # CLR transformation
    clr_transformed = np.zeros_like(abund_pseudo.values, dtype=float)
    for i, (idx, row) in enumerate(abund_pseudo.iterrows()):
        clr_transformed[i, :] = np.log(row.values / geom_mean.iloc[i])
    
    clr_df = pd.DataFrame(
        clr_transformed,
        index=abund_df.index,
        columns=abund_df.columns
    )
    
    print("  ✓ CLR transformation complete")
    print(f"    Mean of CLR values: {clr_df.values.mean():.6f} (should be ~0)")
    print(f"    Std of CLR values: {clr_df.values.std():.4f}")
    
    return clr_df


def create_metadata_mapping(meta_long):
    """Convert long format metadata to wide format for each sample."""
    print("\nProcessing metadata...")

    n_before = len(meta_long)
    meta_long = meta_long.drop_duplicates()
    n_removed = n_before - len(meta_long)
    if n_removed > 0:
        print(f"  Removed {n_removed} duplicate metadata rows")
    
    # Pivot metadata to wide format
    meta_wide = meta_long.pivot_table(
        index='sample_alias',
        columns='metadata_item',
        values='value',
        aggfunc='first'
    )
    
    print(f"  ✓ Metadata wide format: {meta_wide.shape[0]} samples × {meta_wide.shape[1]} metadata fields")
    
    return meta_wide


def align_abundance_metadata(abund_df, meta_df):
    """Align abundance and metadata tables."""
    print("\nAligning abundance table with metadata...")

    if abund_df.index.duplicated().any():
        n_dup = abund_df.index.duplicated().sum()
        print(f"  Removing {n_dup} duplicate abundance sample IDs...")
        abund_df = abund_df[~abund_df.index.duplicated(keep='first')]

    if meta_df.index.duplicated().any():
        n_dup = meta_df.index.duplicated().sum()
        print(f"  Removing {n_dup} duplicate metadata sample IDs...")
        meta_df = meta_df[~meta_df.index.duplicated(keep='first')]
    
    # Find common samples
    common_samples = abund_df.index.intersection(meta_df.index)
    
    print(f"  Samples in abundance table: {abund_df.shape[0]}")
    print(f"  Samples in metadata: {meta_df.shape[0]}")
    print(f"  Common samples: {len(common_samples)}")
    
    # Filter both tables
    abund_aligned = abund_df.loc[common_samples].copy()
    meta_aligned = meta_df.loc[common_samples].copy()
    
    print(f"  ✓ Aligned to {abund_aligned.shape[0]} samples")
    
    return abund_aligned, meta_aligned


def filter_wastewater_samples(abund_df, meta_df, wastewater_studies):
    """Filter for wastewater studies only."""
    print(f"\nFiltering for wastewater studies: {wastewater_studies}...")
    
    # Get study information from metadata
    if 'study_code' in meta_df.columns:
        studies_in_data = meta_df['study_code'].unique()
        print(f"  Studies found in metadata: {studies_in_data}")
        
        # Filter for wastewater studies only
        wastewater_mask = meta_df['study_code'].isin(wastewater_studies)
        meta_filtered = meta_df[wastewater_mask].copy()
        
        print(f"  ✓ Wastewater metadata rows: {len(meta_filtered)}")
        print(f"    Study breakdown:")
        for study in wastewater_studies:
            count = (meta_filtered['study_code'] == study).sum()
            if count > 0:
                print(f"      - {study}: {count}")
        
        # Find common samples between abundance and filtered metadata
        common_samples = abund_df.index.intersection(meta_filtered.index)
        print(f"  Common samples between abundance and metadata: {len(common_samples)}")
        
        # Filter both tables
        abund_filtered = abund_df.loc[common_samples].copy()
        meta_filtered = meta_filtered.loc[common_samples].copy()
        
        print(f"  ✓ Final wastewater samples: {abund_filtered.shape[0]}")
        
        return abund_filtered, meta_filtered
    else:
        print("  ✗ 'study_code' not found in metadata. Returning unfiltered data.")
        return abund_df, meta_df


def save_statistics(stats_df, filtered_species, removed_species, output_dir):
    """Save filtering statistics."""
    report_path = os.path.join(output_dir, 'filtering_statistics.txt')
    
    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("SPECIES FILTERING STATISTICS\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Total species analyzed: {len(stats_df)}\n")
        f.write(f"Species retained: {len(filtered_species)}\n")
        f.write(f"Species removed: {len(removed_species)}\n")
        f.write(f"Retention rate: {len(filtered_species)/len(stats_df)*100:.1f}%\n\n")
        
        f.write("RETAINED SPECIES (Top 20 by mean abundance):\n")
        f.write("-" * 80 + "\n")
        retained_stats = stats_df.loc[filtered_species].sort_values('mean_abundance', ascending=False)
        f.write(retained_stats.head(20).to_string())
        f.write("\n\n")
        
        f.write("REMOVED SPECIES (Reasons):\n")
        f.write("-" * 80 + "\n")
        removed_stats = stats_df.loc[removed_species]
        f.write(f"Low prevalence species: {(removed_stats['prevalence'] < 0.10).sum()}\n")
        f.write(f"Low abundance species: {(removed_stats['mean_abundance'] < 0.01).sum()}\n\n")
        
        f.write("Top 10 removed species (by max abundance):\n")
        f.write(removed_stats.sort_values('max_abundance', ascending=False).head(10).to_string())
    
    print(f"  ✓ Statistics saved to {report_path}")


# ==========================================
# 3. MAIN ANALYSIS PIPELINE
# ==========================================

try:
    # Load data - combine both Schulz-specific and environmental datasets
    print("\nLoading data...")
    
    # Load Schulz data
    meta_schulz_long = pd.read_csv(
        os.path.join(input_dir, 'metadata_Schulz_2017_wastewater_environmental_core_long.tsv'),
        sep='\t'
    )
    abund_schulz_long = pd.read_csv(
        os.path.join(input_dir, 'metaphlan4_Schulz_2017_wastewater_2026-02-06.tsv'),
        sep='\t'
    )
    
    print(f"  ✓ Schulz metadata: {meta_schulz_long.shape[0]} rows")
    print(f"  ✓ Schulz abundance: {abund_schulz_long.shape[0]} rows")
    
    # Load environmental metadata and filter for wastewater studies
    print("\n  Loading environmental metadata to filter for wastewater studies...")
    meta_env_wide = pd.read_csv(
        os.path.join(input_dir, 'environmental_extended_wide.tsv'),
        sep='\t'
    )
    
    # Filter for wastewater studies
    wastewater_mask = meta_env_wide['study_code'].isin(WASTEWATER_STUDIES)
    meta_env_wastewater = meta_env_wide[wastewater_mask].copy()
    if meta_env_wastewater['sample_alias'].duplicated().any():
        n_dup = meta_env_wastewater['sample_alias'].duplicated().sum()
        print(f"  Removing {n_dup} duplicate environmental metadata sample rows...")
        meta_env_wastewater = meta_env_wastewater.drop_duplicates(subset='sample_alias', keep='first')
    
    print(f"  ✓ Environmental metadata loaded: {meta_env_wide.shape[0]} total rows")
    print(f"  ✓ Wastewater studies filtered: {meta_env_wastewater.shape[0]} rows")
    print(f"    Study breakdown:")
    for study in WASTEWATER_STUDIES:
        count = (meta_env_wastewater['study_code'] == study).sum()
        if count > 0:
            print(f"      - {study}: {count}")
    
    # Get sample IDs for wastewater studies
    wastewater_sample_ids = set(meta_env_wastewater['sample_alias'].unique())
    print(f"  ✓ Wastewater sample IDs: {len(wastewater_sample_ids)}")
    
    # Load environmental abundance in chunks and filter on-the-fly
    print("\n  Loading environmental abundance (reading in chunks)...")
    abund_env_chunks = []
    chunk_count = 0
    for chunk in pd.read_csv(
        os.path.join(input_dir, 'environmental_metaphlan4_2026-02-06.tsv'),
        sep='\t',
        chunksize=50000
    ):
        # Filter this chunk for wastewater samples
        chunk_filtered = chunk[chunk['sample_alias'].isin(wastewater_sample_ids)]
        if len(chunk_filtered) > 0:
            abund_env_chunks.append(chunk_filtered)
        chunk_count += 1
        if chunk_count % 10 == 0:
            print(f"    ...processed {chunk_count} chunks")
    
    abund_env_wastewater = pd.concat(abund_env_chunks, ignore_index=True) if abund_env_chunks else pd.DataFrame()
    print(f"  ✓ Environmental abundance filtered: {abund_env_wastewater.shape[0]} rows")
    
    # Combine abundance data
    print("\n  Combining abundance tables...")
    abund_combined = pd.concat([abund_schulz_long, abund_env_wastewater], ignore_index=True)
    n_before = len(abund_combined)
    abund_combined = abund_combined.drop_duplicates()
    n_removed = n_before - len(abund_combined)
    if n_removed > 0:
        print(f"  Removed {n_removed} duplicate combined abundance rows")
    print(f"  ✓ Combined abundance: {abund_combined.shape[0]} rows")
    print(f"    Schulz: {abund_schulz_long.shape[0]} rows")
    print(f"    Environmental (wastewater): {abund_env_wastewater.shape[0]} rows")
    
    # Convert to wide format
    abund_wide = transform_long_to_wide(abund_combined)
    
    # Prepare metadata
    print("\n  Preparing metadata...")
    meta_schulz_wide = create_metadata_mapping(meta_schulz_long)
    
    # Combine metadata tables - union by index
    meta_combined = pd.concat([
        meta_schulz_wide,
        meta_env_wastewater.set_index('sample_alias')
    ], ignore_index=False)
    if meta_combined.index.duplicated().any():
        n_dup = meta_combined.index.duplicated().sum()
        print(f"  Removing {n_dup} duplicate metadata sample IDs after merge...")
        meta_combined = meta_combined[~meta_combined.index.duplicated(keep='first')]
    
    meta_wide = meta_combined
    print(f"  ✓ Combined metadata: {meta_wide.shape[0]} samples × {meta_wide.shape[1]} fields")
    
    # Filter for wastewater samples only (already filtered during loading)
    abund_wastewater = abund_wide.copy()
    meta_wastewater = meta_wide.copy()

    # Sample-size normalization (optional, lecture_02 style)
    abund_wastewater, meta_wastewater, normalization_note, normalization_comparison = apply_sample_size_normalization(
        abund_wastewater,
        meta_wastewater
    )
    print(f"\nSample-size normalization: {normalization_note}")
    if normalization_comparison is not None:
        print("  Cohort sizes before/after normalization:")
        for _, row in normalization_comparison.iterrows():
            print(f"    - {row['study_code']}: {row['n_before']} -> {row['n_after']}")
    
    print(f"\nFinal wastewater dataset:")
    print(f"  Samples: {abund_wastewater.shape[0]}")
    print(f"  Species: {abund_wastewater.shape[1]}")
    
    # Calculate initial statistics
    stats_initial = calculate_prevalence_abundance(abund_wastewater)
    
    # Filter species
    abund_filtered, retained_species, removed_species, stats_all = filter_species(
        abund_wastewater, 
        PREVALENCE_THRESHOLD, 
        ABUNDANCE_THRESHOLD
    )
    
    # Apply CLR transformation
    clr_transformed = add_pseudocount_and_clr(abund_filtered, PSEUDOCOUNT)
    
    # Align with metadata
    abund_aligned, meta_aligned = align_abundance_metadata(clr_transformed, meta_wastewater)
    
    # Save statistics
    save_statistics(stats_all, retained_species, removed_species, output_dir)
    
    # Save processed data
    print("\nSaving processed data...")
    abund_filtered.to_csv(os.path.join(output_dir, 'abundance_filtered_raw.csv'))
    clr_transformed.to_csv(os.path.join(output_dir, 'abundance_clr_transformed.csv'))
    abund_aligned.to_csv(os.path.join(output_dir, 'abundance_clr_aligned.csv'))
    meta_aligned.to_csv(os.path.join(output_dir, 'metadata_aligned.csv'))
    if normalization_comparison is not None:
        normalization_comparison.to_csv(
            os.path.join(output_dir, 'sample_size_normalization_comparison.csv'),
            index=False
        )
    print("  ✓ Data files saved")
    
    # ==========================================
    # 4. VISUALIZATIONS
    # ==========================================
    print("\nGenerating visualizations...")
    sns.set_style("whitegrid")
    
    # === VIZ 1: Pie chart of retained vs removed species ===
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Pie chart
    sizes = [len(retained_species), len(removed_species)]
    labels = [f'Retained\n({len(retained_species)})', f'Removed\n({len(removed_species)})']
    colors = ['#2ecc71', '#e74c3c']
    explode = (0.05, 0)
    
    ax1.pie(sizes, labels=labels, colors=colors, explode=explode,
            autopct='%1.1f%%', startangle=90, textprops={'fontsize': 12, 'weight': 'bold'})
    ax1.set_title('Species Retention After Filtering', fontsize=14, weight='bold', pad=20)
    
    # Filtering criteria breakdown
    low_prev = (stats_all['prevalence'] < PREVALENCE_THRESHOLD).sum()
    low_abund = (stats_all['mean_abundance'] < ABUNDANCE_THRESHOLD).sum()
    both = ((stats_all['prevalence'] < PREVALENCE_THRESHOLD) & 
            (stats_all['mean_abundance'] < ABUNDANCE_THRESHOLD)).sum()
    
    removal_reasons = [
        f"Low prevalence\n(<{PREVALENCE_THRESHOLD:.0%})\nonly",
        f"Low abundance\n(<{ABUNDANCE_THRESHOLD}%)\nonly",
        "Both criteria"
    ]
    removal_counts = [
        low_prev - both,
        low_abund - both,
        both
    ]
    
    colors_reason = ['#3498db', '#f39c12', '#9b59b6']
    bars = ax2.bar(removal_reasons, removal_counts, color=colors_reason, edgecolor='black', linewidth=1.5)
    ax2.set_ylabel('Number of Species', fontsize=12, weight='bold')
    ax2.set_title('Reasons for Species Removal', fontsize=14, weight='bold', pad=20)
    ax2.set_ylim(0, max(removal_counts) * 1.15)
    
    # Add count labels on bars
    for bar in bars:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom', fontsize=11, weight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '01_species_retention_pie.png'), dpi=300, bbox_inches='tight')
    print("  ✓ Saved: 01_species_retention_pie.png")
    plt.close()
    
    # === VIZ 2: Scatter plot of prevalence vs abundance (before/after filtering) ===
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Before filtering
    ax1.scatter(stats_all['prevalence'] * 100, stats_all['mean_abundance'], 
               alpha=0.6, s=60, c='#e74c3c', edgecolors='black', linewidth=0.5)
    ax1.axvline(PREVALENCE_THRESHOLD * 100, color='blue', linestyle='--', linewidth=2, label=f'Prevalence threshold ({PREVALENCE_THRESHOLD:.0%})')
    ax1.axhline(ABUNDANCE_THRESHOLD, color='orange', linestyle='--', linewidth=2, label=f'Abundance threshold ({ABUNDANCE_THRESHOLD}%)')
    ax1.set_xlabel('Prevalence (%)', fontsize=12, weight='bold')
    ax1.set_ylabel('Mean Relative Abundance (%)', fontsize=12, weight='bold')
    ax1.set_title(f'Before Filtering ({len(stats_all)} species)', fontsize=13, weight='bold')
    ax1.set_yscale('log')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # After filtering
    stats_retained = stats_all.loc[retained_species]
    colors_after = ['#2ecc71' if (p >= PREVALENCE_THRESHOLD or m >= ABUNDANCE_THRESHOLD) else '#e74c3c' 
                   for p, m in zip(stats_retained['prevalence']*100, stats_retained['mean_abundance'])]
    ax2.scatter(stats_retained['prevalence'] * 100, stats_retained['mean_abundance'], 
               alpha=0.6, s=60, c='#2ecc71', edgecolors='black', linewidth=0.5)
    ax2.axvline(PREVALENCE_THRESHOLD * 100, color='blue', linestyle='--', linewidth=2, label=f'Prevalence threshold')
    ax2.axhline(ABUNDANCE_THRESHOLD, color='orange', linestyle='--', linewidth=2, label=f'Abundance threshold')
    ax2.set_xlabel('Prevalence (%)', fontsize=12, weight='bold')
    ax2.set_ylabel('Mean Relative Abundance (%)', fontsize=12, weight='bold')
    ax2.set_title(f'After Filtering ({len(stats_retained)} species)', fontsize=13, weight='bold')
    ax2.set_yscale('log')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '02_prevalence_abundance_filtering.png'), dpi=300, bbox_inches='tight')
    print("  ✓ Saved: 02_prevalence_abundance_filtering.png")
    plt.close()
    
    # === VIZ 3: Boxplots of CLR-transformed abundances (Top 12 species) ===
    top_species = stats_retained.nlargest(12, 'mean_abundance').index.tolist()
    
    fig, axes = plt.subplots(3, 4, figsize=(18, 12))
    axes = axes.flatten()
    
    for idx, species in enumerate(top_species):
        species_name = species.replace('s__', '')
        clr_values = clr_transformed[species]
        
        bp = axes[idx].boxplot([clr_values], vert=True, patch_artist=True, widths=0.6)
        bp['boxes'][0].set_facecolor('#3498db')
        bp['boxes'][0].set_alpha(0.7)
        
        # Add individual points
        y = clr_values
        x = np.random.normal(1, 0.04, size=len(y))
        axes[idx].scatter(x, y, alpha=0.4, s=40, color='black')
        
        axes[idx].set_ylabel('CLR-transformed Abundance', fontsize=10, weight='bold')
        axes[idx].set_title(species_name[:40], fontsize=10, weight='bold')
        axes[idx].grid(True, alpha=0.3, axis='y')
        axes[idx].set_xticklabels([])
        
        # Add statistics
        mean_val = clr_values.mean()
        median_val = clr_values.median()
        axes[idx].axhline(mean_val, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label=f'Mean: {mean_val:.2f}')
        axes[idx].axhline(median_val, color='green', linestyle='--', linewidth=1.5, alpha=0.7, label=f'Median: {median_val:.2f}')
        axes[idx].legend(fontsize=8)
    
    plt.suptitle('CLR-Transformed Abundance Distributions (Top 12 Species)', 
                fontsize=16, weight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '03_clr_boxplots_top12.png'), dpi=300, bbox_inches='tight')
    print("  ✓ Saved: 03_clr_boxplots_top12.png")
    plt.close()
    
    # ===  VIZ 4: Histogram and density plots of CLR values ===
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # All CLR values distribution
    all_clr_values = clr_transformed.values.flatten()
    axes[0, 0].hist(all_clr_values, bins=50, color='#3498db', alpha=0.7, edgecolor='black')
    axes[0, 0].axvline(all_clr_values.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {all_clr_values.mean():.4f}')
    axes[0, 0].set_xlabel('CLR Value', fontsize=11, weight='bold')
    axes[0, 0].set_ylabel('Frequency', fontsize=11, weight='bold')
    axes[0, 0].set_title('Distribution of All CLR-Transformed Values', fontsize=12, weight='bold')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Density plot
    axes[0, 1].hist(all_clr_values, bins=50, density=True, alpha=0.5, color='#3498db', edgecolor='black')
    try:
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(all_clr_values)
        x_range = np.linspace(all_clr_values.min(), all_clr_values.max(), 200)
        axes[0, 1].plot(x_range, kde(x_range), 'r-', linewidth=2.5, label='KDE')
    except:
        pass
    axes[0, 1].set_xlabel('CLR Value', fontsize=11, weight='bold')
    axes[0, 1].set_ylabel('Density', fontsize=11, weight='bold')
    axes[0, 1].set_title('CLR Distribution with Kernel Density Estimate', fontsize=12, weight='bold')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Q-Q plot
    from scipy import stats as sp_stats
    sp_stats.probplot(all_clr_values, dist="norm", plot=axes[1, 0])
    axes[1, 0].set_title('Q-Q Plot (Normality Assessment)', fontsize=12, weight='bold')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Sample-wise statistics
    sample_means = clr_transformed.mean(axis=1)
    sample_stds = clr_transformed.std(axis=1)
    
    axes[1, 1].scatter(sample_means, sample_stds, alpha=0.6, s=80, color='#2ecc71', edgecolors='black', linewidth=1)
    axes[1, 1].set_xlabel('Mean CLR per Sample', fontsize=11, weight='bold')
    axes[1, 1].set_ylabel('Std Dev CLR per Sample', fontsize=11, weight='bold')
    axes[1, 1].set_title('Sample-wise CLR Statistics', fontsize=12, weight='bold')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '04_clr_distribution_analysis.png'), dpi=300, bbox_inches='tight')
    print("  ✓ Saved: 04_clr_distribution_analysis.png")
    plt.close()
    
    # === VIZ 5: Heatmap of top species in filtered dataset ===
    top_n_species = 20
    top_species_for_heatmap = stats_retained.nlargest(top_n_species, 'mean_abundance').index.tolist()
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    heatmap_data = clr_transformed[top_species_for_heatmap].T
    
    # Clean species names for display
    display_names = [s.replace('s__', '').replace('_', ' ')[:50] for s in top_species_for_heatmap]
    
    sns.heatmap(heatmap_data, 
               cmap='RdBu_r', 
               center=0,
               cbar_kws={'label': 'CLR-Transformed Abundance'},
               ax=ax,
               yticklabels=display_names,
               xticklabels=False)
    
    ax.set_xlabel('Samples', fontsize=12, weight='bold')
    ax.set_ylabel('Species', fontsize=12, weight='bold')
    ax.set_title(f'Top {top_n_species} Species - CLR-Transformed Abundances', fontsize=14, weight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '05_clr_heatmap_top20.png'), dpi=300, bbox_inches='tight')
    print("  ✓ Saved: 05_clr_heatmap_top20.png")
    plt.close()
    
    # === VIZ 6: Summary statistics comparison (before vs after) ===
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Number of species per sample (before/after filtering)
    species_per_sample_before = (abund_wide > 0).sum(axis=1)
    species_per_sample_after = (abund_filtered > 0).sum(axis=1)
    
    axes[0, 0].hist([species_per_sample_before, species_per_sample_after], 
                   bins=15, label=['Before', 'After'], color=['#e74c3c', '#2ecc71'], alpha=0.7, edgecolor='black')
    axes[0, 0].set_xlabel('Number of Detected Species per Sample', fontsize=11, weight='bold')
    axes[0, 0].set_ylabel('Frequency', fontsize=11, weight='bold')
    axes[0, 0].set_title('Species Richness: Before vs After Filtering', fontsize=12, weight='bold')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3, axis='y')
    
    # Total relative abundance per sample (distribution before/after)
    total_abund_before = abund_wide.sum(axis=1).values
    total_abund_after = abund_filtered.sum(axis=1).values

    bins_total = np.linspace(
        min(total_abund_before.min(), total_abund_after.min()),
        max(total_abund_before.max(), total_abund_after.max()),
        20
    )
    axes[0, 1].hist(total_abund_before, bins=bins_total, alpha=0.55,
                    color='#3498db', edgecolor='black', label='Before')
    axes[0, 1].hist(total_abund_after, bins=bins_total, alpha=0.55,
                    color='#9b59b6', edgecolor='black', label='After')
    axes[0, 1].axvline(total_abund_before.mean(), color='#1f77b4', linestyle='--', linewidth=2,
                      label=f"Before mean={total_abund_before.mean():.1f}")
    axes[0, 1].axvline(total_abund_after.mean(), color='#6c3483', linestyle='--', linewidth=2,
                      label=f"After mean={total_abund_after.mean():.1f}")
    axes[0, 1].set_xlabel('Total Relative Abundance (%)', fontsize=11, weight='bold')
    axes[0, 1].set_ylabel('Frequency', fontsize=11, weight='bold')
    axes[0, 1].set_title('Total Abundance per Sample', fontsize=12, weight='bold')
    axes[0, 1].legend(fontsize=9)
    axes[0, 1].grid(True, alpha=0.3, axis='y')
    
    # Species abundance distribution (log10 transformed non-zero values)
    before_nonzero = abund_wide.values[abund_wide.values > 0]
    after_nonzero = abund_filtered.values[abund_filtered.values > 0]
    log_before = np.log10(before_nonzero)
    log_after = np.log10(after_nonzero)

    bins_log = np.linspace(min(log_before.min(), log_after.min()),
                           max(log_before.max(), log_after.max()), 35)
    axes[1, 0].hist(log_before, bins=bins_log, alpha=0.55,
                    color='#f39c12', edgecolor='black', label='Before')
    axes[1, 0].hist(log_after, bins=bins_log, alpha=0.55,
                    color='#27ae60', edgecolor='black', label='After')
    axes[1, 0].axvline(np.median(log_before), color='#d35400', linestyle='--', linewidth=2,
                      label=f"Before median={np.median(log_before):.2f}")
    axes[1, 0].axvline(np.median(log_after), color='#1e8449', linestyle='--', linewidth=2,
                      label=f"After median={np.median(log_after):.2f}")
    axes[1, 0].set_xlabel('log10(Relative Abundance %)', fontsize=11, weight='bold')
    axes[1, 0].set_ylabel('Frequency', fontsize=11, weight='bold')
    axes[1, 0].set_title('Distribution of Non-zero Abundances', fontsize=12, weight='bold')
    axes[1, 0].legend(fontsize=9)
    axes[1, 0].grid(True, alpha=0.3, axis='y')
    
    # Sparsity comparison
    sparsity_before = (abund_wide == 0).sum().sum() / (abund_wide.shape[0] * abund_wide.shape[1]) * 100
    sparsity_after = (abund_filtered == 0).sum().sum() / (abund_filtered.shape[0] * abund_filtered.shape[1]) * 100
    
    categories = ['Before\nFiltering', 'After\nFiltering']
    sparsity_values = [sparsity_before, sparsity_after]
    bars = axes[1, 1].bar(categories, sparsity_values, color=['#e74c3c', '#2ecc71'], alpha=0.7, edgecolor='black', linewidth=2)
    axes[1, 1].set_ylabel('Sparsity (%)', fontsize=11, weight='bold')
    axes[1, 1].set_title('Data Sparsity: Zeros in Matrix', fontsize=12, weight='bold')
    axes[1, 1].set_ylim(0, 100)
    
    for bar, val in zip(bars, sparsity_values):
        height = bar.get_height()
        axes[1, 1].text(bar.get_x() + bar.get_width()/2., height,
                       f'{val:.1f}%',
                       ha='center', va='bottom', fontsize=12, weight='bold')
    
    axes[1, 1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '06_summary_statistics_comparison.png'), dpi=300, bbox_inches='tight')
    print("  ✓ Saved: 06_summary_statistics_comparison.png")
    plt.close()
    
    # ==========================================
    # 5. GENERATE SUMMARY REPORT
    # ==========================================
    
    report_path = os.path.join(output_dir, 'analysis_summary_report.txt')
    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("MICROBIOME DATA FILTERING & CLR TRANSFORMATION - ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("1. FILTERING PARAMETERS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Prevalence threshold: {PREVALENCE_THRESHOLD:.1%} (species in ≥ x% of samples)\n")
        f.write(f"Abundance threshold: {ABUNDANCE_THRESHOLD}% (mean relative abundance ≥ x%)\n")
        f.write(f"Pseudocount: {PSEUDOCOUNT}\n\n")

        f.write("Sample-size normalization\n")
        f.write("-" * 40 + "\n")
        f.write(f"Applied: {APPLY_SAMPLE_SIZE_NORMALIZATION}\n")
        f.write(f"Rule: {NORMALIZE_STUDY} -> {NORMALIZED_N} max samples\n")
        f.write(f"Status: {normalization_note}\n")
        if normalization_comparison is not None:
            f.write("Cohort sizes (before -> after):\n")
            for _, row in normalization_comparison.iterrows():
                f.write(f"  - {row['study_code']}: {row['n_before']} -> {row['n_after']}\n")
        f.write("\n")
        
        f.write("2. FILTERING RESULTS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Initial samples: {abund_wide.shape[0]}\n")
        f.write(f"Initial species: {abund_wide.shape[1]}\n")
        f.write(f"Filtered samples: {abund_filtered.shape[0]} (Δ: {abund_wide.shape[0] - abund_filtered.shape[0]} removed)\n")
        f.write(f"Filtered species: {abund_filtered.shape[1]} (Δ: {abund_wide.shape[1] - abund_filtered.shape[1]} removed)\n")
        f.write(f"Species retention rate: {len(retained_species)/len(stats_all)*100:.1f}%\n\n")
        
        f.write("3. SPARSITY ANALYSIS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Before filtering:\n")
        f.write(f"  - Sparsity: {sparsity_before:.2f}%\n")
        f.write(f"  - Min detected species per sample: {species_per_sample_before.min()}\n")
        f.write(f"  - Max detected species per sample: {species_per_sample_before.max()}\n")
        f.write(f"  - Mean detected species per sample: {species_per_sample_before.mean():.1f}\n")
        
        f.write(f"\nAfter filtering:\n")
        f.write(f"  - Sparsity: {sparsity_after:.2f}%\n")
        f.write(f"  - Min detected species per sample: {species_per_sample_after.min()}\n")
        f.write(f"  - Max detected species per sample: {species_per_sample_after.max()}\n")
        f.write(f"  - Mean detected species per sample: {species_per_sample_after.mean():.1f}\n\n")
        
        f.write("4. CLR TRANSFORMATION VALIDATION\n")
        f.write("-" * 80 + "\n")
        f.write(f"CLR mean: {all_clr_values.mean():.8f} (should be ≈ 0)\n")
        f.write(f"CLR std: {all_clr_values.std():.4f}\n")
        f.write(f"CLR min: {all_clr_values.min():.4f}\n")
        f.write(f"CLR max: {all_clr_values.max():.4f}\n")
        f.write(f"CLR median: {np.median(all_clr_values):.4f}\n\n")
        
        f.write("5. METADATA ALIGNMENT\n")
        f.write("-" * 80 + "\n")
        f.write(f"Aligned samples: {abund_aligned.shape[0]}\n")
        f.write(f"Aligned species: {abund_aligned.shape[1]}\n")
        f.write(f"Metadata fields: {meta_aligned.shape[1]}\n\n")
        
        f.write("6. TOP 15 RETAINED SPECIES (by mean abundance)\n")
        f.write("-" * 80 + "\n")
        top_15 = stats_retained.nlargest(15, 'mean_abundance')
        for i, (species, row) in enumerate(top_15.iterrows(), 1):
            f.write(f"{i:2d}. {species}\n")
            f.write(f"    Mean abundance: {row['mean_abundance']:.4f}%\n")
            f.write(f"    Prevalence: {row['prevalence']:.1%}\n")
            f.write(f"    Detected in {row['n_samples_detected']} samples\n\n")
        
        f.write("7. OUTPUT FILES GENERATED\n")
        f.write("-" * 80 + "\n")
        f.write("Data tables:\n")
        f.write("  - abundance_filtered_raw.csv: Filtered abundance table (raw %)\n")
        f.write("  - abundance_clr_transformed.csv: CLR-transformed abundances\n")
        f.write("  - abundance_clr_aligned.csv: CLR abundances aligned with metadata\n")
        f.write("  - metadata_aligned.csv: Aligned metadata\n\n")
        
        f.write("Visualizations:\n")
        f.write("  1. 01_species_retention_pie.png - Species retention & removal reasons\n")
        f.write("  2. 02_prevalence_abundance_filtering.png - Scatter plot of filtering criteria\n")
        f.write("  3. 03_clr_boxplots_top12.png - Boxplots of top 12 species distributions\n")
        f.write("  4. 04_clr_distribution_analysis.png - CLR value distributions and normality\n")
        f.write("  5. 05_clr_heatmap_top20.png - Heatmap of top 20 species\n")
        f.write("  6. 06_summary_statistics_comparison.png - Before/after filtering comparison\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("Analysis complete!\n")
        f.write("=" * 80 + "\n")
    
    print(f"\n✓ Summary report saved to {report_path}")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)
    print(f"\nAll outputs saved to: {output_dir}")
    print("\nKey findings:")
    print(f"  • Retained {len(retained_species)}/{len(stats_all)} species ({len(retained_species)/len(stats_all)*100:.1f}%)")
    print(f"  • CLR transformation centered at {all_clr_values.mean():.6f} (±{all_clr_values.std():.4f})")
    print(f"  • Aligned {abund_aligned.shape[0]} samples across abundance and metadata tables")
    
except Exception as e:
    print(f"\n✗ ERROR: {e}")
    import traceback
    traceback.print_exc()
