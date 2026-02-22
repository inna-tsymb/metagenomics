#!/usr/bin/env python3
"""
Association Statements: Formulate Biological Interpretations
================================================================
Generate human-readable statements about species-group associations,
indicating positive/negative associations with specific wastewater sources.
"""

import pandas as pd
import numpy as np
import os

# ============================================================================
# CONFIGURATION
# ============================================================================
ABUNDANCE_FILE = os.getenv('ABUNDANCE_FILE', 'output/filtering_clr_analysis/abundance_clr_aligned.csv')
METADATA_FILE = os.getenv('METADATA_FILE', 'output/filtering_clr_analysis/metadata_aligned.csv')
RESULTS_FILE = os.getenv('RESULTS_FILE', 'output/association_analysis/association_results_all.csv')
OUTPUT_DIR = os.getenv('OUTPUT_DIR', 'output/association_analysis')

GROUP_VAR = 'study_code'
N_TOP_SPECIES = 30  # Generate statements for top N species
FDR_THRESHOLD = 0.05

# ============================================================================
# LOAD DATA
# ============================================================================
print("Loading data...")
abundance = pd.read_csv(ABUNDANCE_FILE, index_col=0)
metadata = pd.read_csv(METADATA_FILE, index_col=0)
results = pd.read_csv(RESULTS_FILE)

# Clean abundance (remove duplicate sample IDs)
if abundance.index.duplicated().any():
    n_dup = abundance.index.duplicated().sum()
    print(f"Removing {n_dup} duplicate samples from abundance...")
    abundance = abundance[~abundance.index.duplicated(keep='first')].copy()

# Clean metadata (remove duplicates)
metadata_clean = metadata[~metadata.index.duplicated(keep='first')].copy()

# Clean results (keep strongest evidence per species)
if results.duplicated(subset=['species']).any():
    n_dup = results.duplicated(subset=['species']).sum()
    print(f"Removing {n_dup} duplicate species rows from association results...")
    results = results.sort_values('p_value').drop_duplicates(subset=['species'], keep='first')

# Align samples
common_samples = abundance.index.intersection(metadata_clean.index)
abundance = abundance.loc[common_samples]
metadata = metadata_clean.loc[common_samples]

print(f"Samples aligned: {len(common_samples)}")

# ============================================================================
# FUNCTION: GENERATE ASSOCIATION STATEMENTS
# ============================================================================
def get_group_means(species_name, group_var, abundance_df, metadata_df):
    """
    Calculate mean abundance for each group.
    Returns dict with group-level statistics.
    """
    species_data = abundance_df[species_name]
    groups = metadata_df[group_var]
    
    group_stats = {}
    for group in groups.unique():
        mask = groups == group
        group_values = species_data[mask]
        
        group_stats[group] = {
            'mean': group_values.mean(),
            'median': group_values.median(),
            'std': group_values.std(),
            'n': mask.sum(),
            'values': group_values
        }
    
    return group_stats

def calculate_beta(species_name, group_var, abundance_df, metadata_df, reference_group=None):
    """
    Calculate effect sizes (beta) as deviations from mean abundance.
    Beta > 0 = enriched in this group
    Beta < 0 = depleted in this group
    """
    group_stats = get_group_means(species_name, group_var, abundance_df, metadata_df)
    overall_mean = abundance_df[species_name].mean()
    
    betas = {}
    for group, stats in group_stats.items():
        betas[group] = {
            'beta': stats['mean'] - overall_mean,
            'mean_abundance': stats['mean'],
            'n_samples': stats['n'],
            'std': stats['std']
        }
    
    return betas, group_stats, overall_mean

def determine_direction(species_name, group_var, abundance_df, metadata_df):
    """
    Determine if species is positively or negatively associated with each group.
    """
    betas, group_stats, overall_mean = calculate_beta(
        species_name, group_var, abundance_df, metadata_df
    )
    
    # Find max and min enriched groups
    sorted_betas = sorted(betas.items(), key=lambda x: x[1]['beta'], reverse=True)
    enriched_group = sorted_betas[0]
    depleted_group = sorted_betas[-1]
    
    return betas, enriched_group, depleted_group, group_stats

# ============================================================================
# GENERATE STATEMENTS FOR TOP SPECIES
# ============================================================================
print(f"\n{'='*80}")
print("GENERATING ASSOCIATION STATEMENTS FOR TOP SPECIES")
print(f"{'='*80}\n")

# Sort by p-value and get top species
top_species_results = results.nsmallest(N_TOP_SPECIES, 'p_value')

statements = []

for idx, (_, row) in enumerate(top_species_results.iterrows(), 1):
    species = row['species']
    p_value = row['p_value']
    p_adjusted = row['p_adjusted']
    effect_size = row['effect_size']
    f_stat = row['f_statistic']
    r_squared = row['r_squared']
    
    # Get direction of association
    betas, enriched_group, depleted_group, group_stats = determine_direction(
        species, GROUP_VAR, abundance, metadata
    )
    
    enriched_name, enriched_beta = enriched_group
    depleted_name, depleted_beta = depleted_group
    
    # Calculate effect size magnitude (percentage change)
    overall_mean = abundance[species].mean()
    if overall_mean != 0:
        pct_enriched = (enriched_beta['beta'] / overall_mean) * 100
        pct_depleted = (depleted_beta['beta'] / overall_mean) * 100
    else:
        pct_enriched = 0
        pct_depleted = 0
    
    # Formulate statement
    statement = f"""
{idx}. {species.replace('s__', '').replace('_', ' ')}

   Association with wastewater source (overall):
   - p-value: {p_value:.2e}
   - FDR-corrected p-value: {p_adjusted:.2e}
   - Effect size (max-min abundance): {effect_size:.4f} (CLR units)
   - Model R²: {r_squared:.4f} (Study source explains {r_squared*100:.1f}% of variation)
   - F-statistic: {f_stat:.2f}

   Group-specific associations:
   ✓ POSITIVELY associated with {enriched_name}:
     - Mean CLR abundance: {enriched_beta['mean_abundance']:.4f}
     - Deviation from global mean: +{enriched_beta['beta']:.4f} (≈ +{pct_enriched:.1f}%)
     - Present in {enriched_beta['n_samples']} samples
     
   ✗ NEGATIVELY associated with {depleted_name}:
     - Mean CLR abundance: {depleted_beta['mean_abundance']:.4f}
     - Deviation from global mean: {depleted_beta['beta']:.4f} (≈ {pct_depleted:.1f}%)
     - Present in {depleted_beta['n_samples']} samples

   Interpretation:
   {species.replace('s__', '').split('_')[0]} shows strong association with wastewater source
   (p={p_value:.2e}), with notably higher relative abundance in {enriched_name}
   compared to other sources. The estimated effect size is {'substantial' if effect_size > 5 else 'moderate'} 
   (beta={enriched_beta['beta']:.4f} CLR units), indicating a {('biologically meaningful' if effect_size > 5 else 'subtle')}
   shift in community composition.
"""
    
    statements.append(statement)
    print(statement)

# ============================================================================
# SAVE DETAILED REPORT
# ============================================================================
report_content = f"""
ASSOCIATION STATEMENTS: Species-Group Relationships
{'='*80}

Analysis of top {N_TOP_SPECIES} species most strongly associated with wastewater source.

Interpretation Guide:
- POSITIVE association: Species is significantly enriched in a particular wastewater source
- NEGATIVE association: Species is significantly depleted in a particular wastewater source
- Beta (effect size): CLR-space units of abundance difference from global mean
  • Beta > 0 indicates enrichment
  • Beta < 0 indicates depletion
  • Magnitude indicates strength of deviation

FDR Significance Threshold: {FDR_THRESHOLD}

{'='*80}
"""

for stmt in statements:
    report_content += stmt

# Save report
report_file = f'{OUTPUT_DIR}/association_statements.txt'
with open(report_file, 'w') as f:
    f.write(report_content)

print(f"\n{'='*80}")
print(f"✓ Saved association statements: {report_file}")
print(f"{'='*80}")

# ============================================================================
# GENERATE SUMMARY TABLE
# ============================================================================
print("\n" + "="*80)
print("SUMMARY TABLE: TOP SPECIES ASSOCIATIONS")
print("="*80 + "\n")

summary_rows = []
for _, row in top_species_results.iterrows():
    species = row['species']
    
    betas, enriched_group, depleted_group, _ = determine_direction(
        species, GROUP_VAR, abundance, metadata
    )
    
    enriched_name, enriched_stats = enriched_group
    
    summary_rows.append({
        'Species': species.replace('s__', '').replace('_', ' '),
        'Enriched in': enriched_name,
        'p-value': f"{row['p_value']:.2e}",
        'Effect Size': f"{row['effect_size']:.3f}",
        'β (enriched)': f"{enriched_stats['beta']:.3f}",
        'R²': f"{row['r_squared']:.3f}"
    })

summary_df = pd.DataFrame(summary_rows)
summary_file = f'{OUTPUT_DIR}/association_summary_table.csv'
summary_df.to_csv(summary_file, index=False)
print(summary_df.to_string(index=False))
print(f"\n✓ Saved summary table: {summary_file}")

# ============================================================================
# GENERATE MARKDOWN VERSION FOR EASY READING
# ============================================================================
markdown_content = f"""# Species-Group Association Statements

## Summary Statistics

- **Total species analyzed:** {len(results)}
- **Significant associations (FDR < {FDR_THRESHOLD}):** {(results['sig_fdr']).sum()}
- **Statements generated for:** Top {N_TOP_SPECIES} species

## Interpretation Guide

- **Positive Association (✓):** Species is enriched (higher mean abundance) in this group
- **Negative Association (✗):** Species is depleted (lower mean abundance) in this group
- **β (Beta):** Effect size in CLR units; deviation from global mean abundance
  - β > 0 = enriched relative to average
  - β < 0 = depleted relative to average
- **Effect size range:** {results['effect_size'].min():.2f} to {results['effect_size'].max():.2f} CLR units
- **Mean R²:** {results['r_squared'].mean():.3f} (average variance explained by group membership)

---

## Top Associated Species

"""

for idx, (_, row) in enumerate(top_species_results.iterrows(), 1):
    species = row['species']
    p_value = row['p_value']
    p_adjusted = row['p_adjusted']
    effect_size = row['effect_size']
    r_squared = row['r_squared']
    
    betas, enriched_group, depleted_group, _ = determine_direction(
        species, GROUP_VAR, abundance, metadata
    )
    
    enriched_name, enriched_stats = enriched_group
    depleted_name, depleted_stats = depleted_group
    
    sig_mark = "**" if row['sig_fdr'] else ""
    
    markdown_content += f"""
### {idx}. {sig_mark}{species.replace('s__', '').replace('_', ' ')}{sig_mark}

**Statistical Significance:** p = {p_value:.2e} | FDR = {p_adjusted:.2e} | Effect size = {effect_size:.3f} CLR units | R² = {r_squared:.3f}

| Metric | Value |
|--------|-------|
| **Enriched in** | {enriched_name} |
| **β (enriched)** | {enriched_stats['beta']:+.4f} |
| **Mean abundance (enriched)** | {enriched_stats['mean_abundance']:.4f} |
| **Depleted in** | {depleted_name} |
| **β (depleted)** | {depleted_stats['beta']:+.4f} |
| **Mean abundance (depleted)** | {depleted_stats['mean_abundance']:.4f} |

**Interpretation:** This species is strongly associated with wastewater source (p={p_value:.2e}). 
It is {'substantially' if effect_size > 5 else 'moderately'} enriched in **{enriched_name}** 
(effect size: {enriched_stats['beta']:+.4f} CLR units) and depleted in **{depleted_name}** 
(effect size: {depleted_stats['beta']:+.4f} CLR units). This species may serve as a biomarker 
for {enriched_name} wastewater.

"""

# Save markdown
markdown_file = f'{OUTPUT_DIR}/association_statements.md'
with open(markdown_file, 'w') as f:
    f.write(markdown_content)

print(f"✓ Saved markdown version: {markdown_file}")

print(f"\n{'='*80}")
print("ALL OUTPUTS COMPLETE!")
print(f"{'='*80}")
