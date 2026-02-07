import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

# ==========================================
# 1. SETUP & LOAD DATA
# ==========================================
base_dir = '/Users/user/Documents/metagenomics/lecture_02'
input_dir = os.path.join(base_dir, 'input')
output_dir = os.path.join(base_dir, 'output', 'advanced_viz')
os.makedirs(output_dir, exist_ok=True)

print("Loading data...")
try:
    meta_df = pd.read_csv(os.path.join(input_dir, 'environmental.tsv'), sep='\t')
    mp_path = os.path.join(input_dir, 'environmental_metaphlan4_2026-02-06.tsv')
    mp_df = pd.read_csv(mp_path, sep='\t')
    
    # Clean Columns
    mp_df.columns = [c.replace('#', '').strip() for c in mp_df.columns]
    
    # Standardize Format (Pivot Long to Wide if needed)
    if 'rel_abund' in mp_df.columns:
        if 'clade_name' in mp_df.columns: tax_col = 'clade_name'
        else: tax_col = mp_df.columns[1]
        mp_df = mp_df[mp_df[tax_col].str.contains('s__') & ~mp_df[tax_col].str.contains('t__')]
        abund_df = mp_df.pivot_table(index='sample_alias', columns=tax_col, values='rel_abund', aggfunc='mean').fillna(0)
    else:
        mp_df.rename(columns={mp_df.columns[0]: 'clade_name'}, inplace=True)
        mp_df = mp_df[mp_df['clade_name'].str.contains('s__') & ~mp_df['clade_name'].str.contains('t__')]
        mp_df.set_index('clade_name', inplace=True)
        cols_to_drop = [c for c in mp_df.columns if 'tax' in c.lower() or 'ncbi' in c.lower()]
        mp_df.drop(columns=cols_to_drop, inplace=True)
        abund_df = mp_df.T.apply(pd.to_numeric, errors='coerce').fillna(0)

    # Clean Index
    abund_df.index = abund_df.index.str.replace('_metaphlan', '').str.replace('.fastq', '').str.replace('.fq', '')

    # Merge
    full_df = pd.merge(abund_df, meta_df, left_index=True, right_on='sample_alias', how='inner')
    if len(full_df) == 0:
        full_df = pd.merge(abund_df, meta_df, left_index=True, right_on='ena_ers_sample_id', how='inner')

    # Define Groups
    anthro = ['Rowe_2017_hospital_wastewater', 'Chu_2017_sludge', 'Schulz_2017_wastewater']
    env = ['Chopyk_2020_pond', 'Lekunberri_2018_river_wastewater']
    
    def assign_group(study):
        if study in anthro: return 'Anthropogenic'
        elif study in env: return 'Environmental'
        return 'Other'

    full_df['Source_Type'] = full_df['study_code'].apply(assign_group)
    
    # Filter for target groups only
    final_df = full_df[full_df['Source_Type'].isin(['Anthropogenic', 'Environmental'])].copy()
    print(f"Data Loaded: {len(final_df)} samples.")

except Exception as e:
    print(f"Error: {e}")
    exit()

# ==========================================
# 2. TAXONOMIC STACKED BAR PLOT (GENUS LEVEL)
# ==========================================
print("Generating Stacked Bar Plot...")

# Extract Genus names from species columns
# Species name format: s__Genus_species -> extract "Genus"
species_cols = [c for c in abund_df.columns if c in final_df.columns]
genus_data = {}

for col in species_cols:
    # Try to extract genus: "s__Escherichia_coli" -> "Escherichia"
    try:
        if 's__' in col:
            genus_name = col.split('s__')[1].split('_')[0]
        else:
            genus_name = col.split('_')[0]
            
        # Sum abundance
        if genus_name in genus_data:
            genus_data[genus_name] += final_df[col]
        else:
            genus_data[genus_name] = final_df[col]
    except:
        continue

genus_df = pd.DataFrame(genus_data)
genus_df['Source_Type'] = final_df['Source_Type']

# Calculate Mean Abundance per Group
grouped_genus = genus_df.groupby('Source_Type').mean()

# Keep Top 10 Genera + "Other"
top_10 = grouped_genus.sum().nlargest(10).index
other_abundance = grouped_genus.drop(columns=top_10).sum(axis=1)

plot_data = grouped_genus[top_10].copy()
plot_data['Other'] = other_abundance

# Plot
plot_data.plot(kind='bar', stacked=True, figsize=(10, 7), colormap='tab20')
plt.title('Mean Genus Composition per Group')
plt.ylabel('Relative Abundance (%)')
plt.xticks(rotation=0)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Genus')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, '1_genus_composition.png'))
plt.close()

# ==========================================
# 3. HIERARCHICAL CLUSTERMAP
# ==========================================
print("Generating Clustermap...")

# Select Top 50 Species across the whole dataset
top_50_species = final_df[species_cols].mean().nlargest(50).index
heatmap_data = final_df[top_50_species]

# Add Group Colors
group_colors = final_df['Source_Type'].map({'Anthropogenic': 'red', 'Environmental': 'green'})

# Plot (Standard Scale=None means raw abundance, z_score=None)
# We use metric='braycurtis' if possible, but 'euclidean' is safer for general use without skbio
g = sns.clustermap(heatmap_data.T, 
                   col_colors=group_colors,
                   cmap='viridis', 
                   standard_scale=1, # Scale rows (species) 0-1 to see patterns
                   figsize=(12, 10),
                   dendrogram_ratio=0.15)

g.ax_heatmap.set_title("Top 50 Species Clustering (Columns=Samples)")
plt.savefig(os.path.join(output_dir, '2_hierarchical_clustering.png'))
plt.close()

print(f"Visualizations saved to: {output_dir}")