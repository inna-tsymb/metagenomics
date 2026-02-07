import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import entropy, mannwhitneyu
from sklearn.decomposition import PCA

# ==========================================
# 1. SETUP PATHS
# ==========================================
base_dir = '/Users/user/Documents/metagenomics/lecture_02'
input_dir = os.path.join(base_dir, 'input')
output_dir = os.path.join(base_dir, 'output', 'comparative_analysis')
os.makedirs(output_dir, exist_ok=True)

# ==========================================
# 2. LOAD DATA
# ==========================================
print("Loading data...")
try:
    # 2a. Load Metadata
    meta_df = pd.read_csv(os.path.join(input_dir, 'environmental.tsv'), sep='\t')
    print(f"Metadata loaded: {len(meta_df)} rows.")

    # 2b. Load MetaPhlAn table
    mp_path = os.path.join(input_dir, 'environmental_metaphlan4_2026-02-06.tsv')
    mp_df = pd.read_csv(mp_path, sep='\t')
    
    # Clean Column Names (remove #)
    mp_df.columns = [c.replace('#', '').strip() for c in mp_df.columns]
    print(f"MetaPhlAn raw shape: {mp_df.shape}")

    # 2c. CHECK FORMAT & PIVOT
    # Case 1: Long Format (sample_alias, clade_name, rel_abund)
    if 'rel_abund' in mp_df.columns and 'sample_alias' in mp_df.columns:
        print("Detected LONG format. Pivoting to WIDE format...")
        
        # Identify Taxonomy Column
        if 'clade_name' in mp_df.columns:
            tax_col = 'clade_name'
        else:
            tax_col = mp_df.columns[1] # Guess

        # Filter for Species (s__) and drop Strains (t__)
        mp_df = mp_df[mp_df[tax_col].str.contains('s__') & ~mp_df[tax_col].str.contains('t__')]
        
        # PIVOT with Aggregation (Fixes "Index contains duplicate entries" error)
        # We use mean() to handle any accidental duplicates
        abund_df = mp_df.pivot_table(index='sample_alias', columns=tax_col, values='rel_abund', aggfunc='mean')
        abund_df = abund_df.fillna(0)
        
        print(f"Pivoted successfully. Samples: {abund_df.shape[0]}, Species: {abund_df.shape[1]}")

    # Case 2: Wide Format (Standard)
    else:
        print("Detected WIDE format. Processing...")
        if 'clade_name' not in mp_df.columns:
             mp_df.rename(columns={mp_df.columns[0]: 'clade_name'}, inplace=True)
        
        mp_df = mp_df[mp_df['clade_name'].str.contains('s__') & ~mp_df['clade_name'].str.contains('t__')]
        mp_df.set_index('clade_name', inplace=True)
        
        # Drop non-sample columns
        cols_to_drop = [c for c in mp_df.columns if 'tax' in c.lower() or 'ncbi' in c.lower()]
        mp_df.drop(columns=cols_to_drop, inplace=True)
        
        abund_df = mp_df.T
        abund_df = abund_df.apply(pd.to_numeric, errors='coerce').fillna(0)

    # 2d. MERGE WITH METADATA
    # Clean up index names for better matching
    abund_df.index = abund_df.index.str.replace('_metaphlan', '').str.replace('.fastq', '').str.replace('.fq', '')

    print("Merging with metadata...")
    # Try merging on sample_alias
    full_df = pd.merge(abund_df, meta_df, left_index=True, right_on='sample_alias', how='inner')
    
    # Fallback to ID if alias fails
    if len(full_df) == 0:
        print("Merge on 'sample_alias' failed. Trying 'ena_ers_sample_id'...")
        full_df = pd.merge(abund_df, meta_df, left_index=True, right_on='ena_ers_sample_id', how='inner')

    if len(full_df) == 0:
        print("\nCRITICAL ERROR: Merge failed. No common sample names found.")
        print(f"Abundance Samples: {abund_df.index[:3].tolist()}")
        print(f"Metadata Aliases: {meta_df['sample_alias'].head(3).tolist()}")
        exit()

    print(f"Successfully merged {len(full_df)} samples.")

except Exception as e:
    print(f"Error processing data: {e}")
    exit()

# ==========================================
# 3. DEFINE GROUPS & SUBSET
# ==========================================
print("Defining Case/Control Groups...")

anthropogenic_studies = ['Rowe_2017_hospital_wastewater', 'Chu_2017_sludge']
municipal_study = 'Schulz_2017_wastewater'
environmental_studies = ['Chopyk_2020_pond', 'Lekunberri_2018_river_wastewater']

# Filter dataset
target_mask = full_df['study_code'].isin(anthropogenic_studies + environmental_studies + [municipal_study])
subset_df = full_df[target_mask].copy()

if len(subset_df) == 0:
    print("Error: No samples found from the target studies.")
    exit()

# Downsample Schulz
schulz_samples = subset_df[subset_df['study_code'] == municipal_study]
others = subset_df[subset_df['study_code'] != municipal_study]

if len(schulz_samples) > 20:
    print(f"Balancing: Downsampling {municipal_study} to 20 samples.")
    schulz_subset = schulz_samples.sample(n=20, random_state=42)
else:
    schulz_subset = schulz_samples

final_df = pd.concat([others, schulz_subset])

def assign_group(study):
    if study in anthropogenic_studies or study == municipal_study:
        return 'Anthropogenic'
    elif study in environmental_studies:
        return 'Environmental'
    return 'Other'

final_df['Source_Type'] = final_df['study_code'].apply(assign_group)
print(f"Final Analysis Set: {len(final_df)} samples")
print(final_df['Source_Type'].value_counts())

# ==========================================
# 4. STATISTICAL ANALYSIS & PLOTTING
# ==========================================
# Identify species columns (columns present in abundance table AND final_df)
# Note: after merge, species columns are in final_df. We need to exclude metadata cols.
metadata_cols = meta_df.columns.tolist() + ['Source_Type', 'Shannon', 'Richness', 'PC1', 'PC2']
species_cols = [c for c in abund_df.columns if c in final_df.columns]

# A. Alpha Diversity
final_df['Shannon'] = final_df[species_cols].apply(lambda x: entropy(x[x>0]), axis=1)
final_df['Richness'] = final_df[species_cols].apply(lambda x: (x > 0).sum(), axis=1)

# Mann-Whitney U
group_anthro = final_df[final_df['Source_Type'] == 'Anthropogenic']['Shannon']
group_env = final_df[final_df['Source_Type'] == 'Environmental']['Shannon']

p_val = 1.0
if len(group_anthro) > 0 and len(group_env) > 0:
    stat, p_val = mannwhitneyu(group_anthro, group_env)
    print(f"Mann-Whitney U Result: p-value = {p_val:.5e}")

# Save Stats
with open(os.path.join(output_dir, 'stats_report.txt'), 'w') as f:
    f.write("COMPARATIVE ANALYSIS REPORT\n")
    f.write(f"Hypothesis: Environmental diversity > Anthropogenic.\n")
    f.write(f"Mann-Whitney p-value: {p_val}\n\n")
    f.write(final_df.groupby('Source_Type')[['Shannon', 'Richness']].describe().to_string())

# B. Plots
# Alpha Boxplot
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
sns.boxplot(data=final_df, x='Source_Type', y='Shannon', palette='Set2')
plt.title(f'Shannon Diversity\n(p={p_val:.1e})')
plt.subplot(1, 2, 2)
sns.boxplot(data=final_df, x='Source_Type', y='Richness', palette='Set2')
plt.title('Observed Richness')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, '1_alpha_diversity_comparison.png'))
plt.close()

# Beta PCA
pca = PCA(n_components=2)
if len(final_df) > 2:
    coords = pca.fit_transform(final_df[species_cols])
    final_df['PC1'] = coords[:, 0]
    final_df['PC2'] = coords[:, 1]
    var_exp = pca.explained_variance_ratio_

    plt.figure(figsize=(10, 7))
    sns.scatterplot(data=final_df, x='PC1', y='PC2', hue='Source_Type', style='study_code', s=100, alpha=0.9)
    plt.title(f'PCA: Microbiome Structure\nPC1 ({var_exp[0]:.1%}) - PC2 ({var_exp[1]:.1%})')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '2_beta_diversity_pca.png'))
    plt.close()

print(f"Analysis Complete. Results saved to: {output_dir}")