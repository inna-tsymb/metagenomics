import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.spatial.distance import pdist, squareform

# ==========================================
# 1. SETUP & DATA LOADING
# ==========================================
base_dir = '/Users/user/Documents/metagenomics/lecture_02'
input_dir = os.path.join(base_dir, 'input')
output_dir = os.path.join(base_dir, 'output', 'final_comprehensive_analysis')
os.makedirs(output_dir, exist_ok=True)

print("--- STARTING COMPREHENSIVE ANALYSIS ---")

try:
    # Load Metadata and Abundance Table
    meta_df = pd.read_csv(os.path.join(input_dir, 'environmental.tsv'), sep='\t')
    mp_path = os.path.join(input_dir, 'environmental_metaphlan4_2026-02-06.tsv')
    mp_df = pd.read_csv(mp_path, sep='\t')
    
    # Clean Column Names
    mp_df.columns = [c.replace('#', '').strip() for c in mp_df.columns]

    # --- DATA PROCESSING (Long to Wide conversion if necessary) ---
    if 'rel_abund' in mp_df.columns and 'sample_alias' in mp_df.columns:
        print("Detected LONG format. Pivoting...")
        tax_col = 'clade_name' if 'clade_name' in mp_df.columns else mp_df.columns[1]
        
        # Filter for Species level only (s__)
        mp_df = mp_df[mp_df[tax_col].str.contains('s__') & ~mp_df[tax_col].str.contains('t__')]
        
        # Pivot table (mean aggregation handles duplicates)
        abund_df = mp_df.pivot_table(index='sample_alias', columns=tax_col, values='rel_abund', aggfunc='mean')
        abund_df = abund_df.fillna(0)
    else:
        print("Detected WIDE format...")
        if 'clade_name' not in mp_df.columns: mp_df.rename(columns={mp_df.columns[0]: 'clade_name'}, inplace=True)
        
        mp_df = mp_df[mp_df['clade_name'].str.contains('s__') & ~mp_df['clade_name'].str.contains('t__')]
        mp_df.set_index('clade_name', inplace=True)
        
        # Drop non-numeric taxonomy columns
        cols_to_drop = [c for c in mp_df.columns if 'tax' in c.lower() or 'ncbi' in c.lower()]
        mp_df.drop(columns=cols_to_drop, inplace=True)
        
        abund_df = mp_df.T.apply(pd.to_numeric, errors='coerce').fillna(0)

    # Clean Index Names
    abund_df.index = abund_df.index.str.replace('_metaphlan', '').str.replace('.fastq', '').str.replace('.fq', '')

    # Merge with Metadata
    full_df = pd.merge(abund_df, meta_df, left_index=True, right_on='sample_alias', how='inner')
    if len(full_df) == 0:
        # Fallback merge on ID
        full_df = pd.merge(abund_df, meta_df, left_index=True, right_on='ena_ers_sample_id', how='inner')

    print(f"Merge successful. Total samples: {len(full_df)}")

except Exception as e:
    print(f"CRITICAL ERROR: {e}")
    exit()

# ==========================================
# 2. BALANCING & GROUPING
# ==========================================
# Define target studies
anthropogenic_studies = ['Rowe_2017_hospital_wastewater', 'Chu_2017_sludge']
municipal_study = 'Schulz_2017_wastewater'
environmental_studies = ['Chopyk_2020_pond', 'Lekunberri_2018_river_wastewater']

# Filter dataset for specific studies
target_mask = full_df['study_code'].isin(anthropogenic_studies + environmental_studies + [municipal_study])
subset_df = full_df[target_mask].copy()

# Downsample Schulz study (Balancing logic)
schulz_samples = subset_df[subset_df['study_code'] == municipal_study]
others = subset_df[subset_df['study_code'] != municipal_study]

if len(schulz_samples) > 20:
    print(f"Balancing: Downsampling {municipal_study} to 20 samples.")
    schulz_subset = schulz_samples.sample(n=20, random_state=42)
else:
    schulz_subset = schulz_samples

# Recombine
final_df = pd.concat([others, schulz_subset])

# Assign "Case/Control" groups
def assign_group(study):
    if study in anthropogenic_studies or study == municipal_study:
        return 'Anthropogenic'
    elif study in environmental_studies:
        return 'Environmental'
    return 'Other'

final_df['Source_Type'] = final_df['study_code'].apply(assign_group)

# ==========================================
# 3. SPECIES FILTERING (Cleaning Zero-Sum Columns)
# ==========================================
# Identify species columns that exist in both abundance and metadata DF
potential_species = [c for c in abund_df.columns if c in final_df.columns]
species_data = final_df[potential_species]

# Drop species that are all zeros in this specific subset (113 samples)
active_species = species_data.loc[:, (species_data != 0).any(axis=0)]
species_cols = active_species.columns.tolist()

# Update final DataFrame
meta_cols = [c for c in final_df.columns if c not in potential_species]
final_df = pd.concat([final_df[meta_cols], active_species], axis=1)

print(f"\n--- FINAL DATASET STATISTICS ---")
print(f"Samples: {len(final_df)} (Expected ~113)")
print(f"Active Species: {len(species_cols)}")
print(final_df['Source_Type'].value_counts())

# ==========================================
# 4. BASIC PLOTS (Core, Ratio, Diversity, Accumulation)
# ==========================================
sns.set_style("whitegrid")

# --- A. Core Microbiome ---
print("Generating Core Microbiome Plot...")
mean_abund = final_df[species_cols].mean()
prevalence = (final_df[species_cols] > 0).sum() / len(final_df)
core_df = pd.DataFrame({'Mean': mean_abund, 'Prevalence': prevalence})
# Clean names: s__Name_name -> Name name
core_df['Name'] = core_df.index.map(lambda x: x.split('s__')[1].replace('_', ' ') if 's__' in x else x)
core_df['Type'] = np.where((core_df['Prevalence'] > 0.5) & (core_df['Mean'] > 1.0), 'Core', 'Transient')

fig, ax = plt.subplots(figsize=(12, 9))

# Create scatter plot colored by Mean abundance
scatter = ax.scatter(core_df['Prevalence'], core_df['Mean'], 
                    c=core_df['Mean'], cmap='YlOrRd',
                    s=core_df['Mean'] * 30, alpha=0.8, edgecolors='black', linewidth=0.5,
                    norm=plt.Normalize(vmin=core_df['Mean'].min(), vmax=core_df['Mean'].max()))

# Add colorbar for abundance scale
cbar = plt.colorbar(scatter, ax=ax, label='Mean Abundance (%)', pad=0.02)

# Add legend for Type (using marker styles or add text annotation)
from matplotlib.patches import Patch
core_mask = core_df['Type'] == 'Core'
transient_mask = core_df['Type'] == 'Transient'

# Create custom legend patches
legend_elements = [Patch(facecolor='white', edgecolor='black', linewidth=2, label='Core Species'),
                   Patch(facecolor='white', edgecolor='black', linewidth=1, label='Transient Species')]

# Add text annotation to distinguish them
ax.text(0.02, 0.98, 'Bold labels = Core Species\nRegular labels = Transient Species', 
        transform=ax.transAxes, fontsize=10, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Label all core species and top transient species
labeled = set()
# First, label all core species
for _, row in core_df[core_df['Type'] == 'Core'].iterrows():
    ax.text(row['Prevalence'] + 0.01, row['Mean'], row['Name'], fontsize=9, fontweight='bold')
    labeled.add(row['Name'])

# Then label top 5 transient species by mean abundance
top_transient = core_df[core_df['Type'] == 'Transient'].nlargest(5, 'Mean')
for _, row in top_transient.iterrows():
    if row['Name'] not in labeled:
        ax.text(row['Prevalence'] + 0.01, row['Mean'], row['Name'], fontsize=8, alpha=0.7)
        labeled.add(row['Name'])

ax.set_yscale('log')
ax.set_xlabel('Prevalence (Fraction of Samples)', fontsize=12)
ax.set_ylabel('Mean Abundance (%)', fontsize=12)
ax.set_title('Core Microbiome Analysis: Species Color-Coded by Abundance (Balanced Subset)', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, '1_Core_Microbiome.png'), dpi=300, bbox_inches='tight')
plt.close()

# --- B. P/B Ratio ---
print("Generating P/B Ratio...")
proteo = [c for c in species_cols if 'p__Proteobacteria' in c]
bacter = [c for c in species_cols if 'p__Bacteroidetes' in c]

if proteo and bacter:
    final_df['Proteo_Sum'] = final_df[proteo].sum(axis=1)
    final_df['Bacter_Sum'] = final_df[bacter].sum(axis=1)
    final_df['PB_Ratio'] = final_df['Proteo_Sum'] / (final_df['Bacter_Sum'] + 0.001)

    plt.figure(figsize=(6, 6))
    sns.boxplot(data=final_df, x='Source_Type', y='PB_Ratio', hue='Source_Type', palette='Set2', legend=False)
    sns.stripplot(data=final_df, x='Source_Type', y='PB_Ratio', color='k', alpha=0.5)
    plt.yscale('log')
    plt.title('Proteobacteria / Bacteroidetes Ratio')
    plt.savefig(os.path.join(output_dir, '2_PB_Ratio.png'))
    plt.close()

# --- C. Diversity Violins ---
print("Generating Violins...")
def calc_shannon(row):
    p = row[row > 0] / row.sum()
    return -np.sum(p * np.log(p))
final_df['Shannon'] = final_df[species_cols].apply(calc_shannon, axis=1)

plt.figure(figsize=(8, 6))
sns.violinplot(data=final_df, x='Source_Type', y='Shannon', hue='Source_Type', palette='muted', inner='quartile', legend=False)
sns.stripplot(data=final_df, x='Source_Type', y='Shannon', color='k', alpha=0.3)
plt.title('Shannon Diversity Distribution')
plt.savefig(os.path.join(output_dir, '3_Shannon_Violin.png'))
plt.close()

# --- D. Overlap Barplot ---
print("Generating Overlap Barplot...")
anthro_sp = final_df[final_df['Source_Type']=='Anthropogenic'][species_cols].sum()
env_sp = final_df[final_df['Source_Type']=='Environmental'][species_cols].sum()
set_a = set(anthro_sp[anthro_sp > 0].index)
set_e = set(env_sp[env_sp > 0].index)

counts = [len(set_a - set_e), len(set_a & set_e), len(set_e - set_a)]
labels = ['Unique Anthro', 'Shared', 'Unique Env']
plt.figure(figsize=(8,6))
bars = plt.bar(labels, counts, color=['#d62728', 'gray', '#2ca02c'], edgecolor='k')
plt.bar_label(bars)
plt.title(f'Species Overlap (Total Active Species: {len(species_cols)})')
plt.savefig(os.path.join(output_dir, '4_Overlap.png'))
plt.close()

# --- E. Accumulation Curve ---
print("Generating Accumulation Curve...")
plt.figure(figsize=(10, 6))
for group in ['Anthropogenic', 'Environmental']:
    grp_data = final_df[final_df['Source_Type'] == group][species_cols]
    n = len(grp_data)
    acc = np.zeros(n)
    for _ in range(10): # 10 permutations
        shuffled = grp_data.sample(frac=1).values
        acc += np.maximum.accumulate(shuffled > 0, axis=0).sum(axis=1)
    plt.plot(range(1, n+1), acc/10, label=group, lw=2)
plt.legend()
plt.title('Species Accumulation Curve')
plt.xlabel('Number of Samples')
plt.ylabel('Species Discovered')
plt.savefig(os.path.join(output_dir, '5_Accumulation.png'))
plt.close()

# ==========================================
# 5. SHARED SPECIES DEEP DIVE
# ==========================================
print("--- STARTING SHARED SPECIES DEEP DIVE ---")

# Calculate means per group
group_means = final_df.groupby('Source_Type')[species_cols].mean()
mean_anthro = group_means.loc['Anthropogenic']
mean_env = group_means.loc['Environmental']

# Identify shared species (Mean > 0 in both groups)
shared_mask = (mean_anthro > 0) & (mean_env > 0)
shared_species = mean_anthro[shared_mask].index.tolist()
print(f"Shared Species Count: {len(shared_species)}")

# Create Report DataFrame
shared_df = pd.DataFrame({
    'Species': shared_species,
    'Mean_Anthro': mean_anthro[shared_species],
    'Mean_Env': mean_env[shared_species]
})

# Add log fold change
shared_df['Log2_FoldChange'] = np.log2((shared_df['Mean_Anthro'] + 1e-6) / (shared_df['Mean_Env'] + 1e-6))
shared_df['Classification'] = np.where(abs(shared_df['Log2_FoldChange']) < 1, 'Generalist (Equal)', 
                                np.where(shared_df['Log2_FoldChange'] > 0, 'Anthro-Biased', 'Env-Biased'))
shared_df['Short_Name'] = shared_df['Species'].apply(lambda x: x.split('s__')[1].replace('_', ' ') if 's__' in x else x)

# Save Report
csv_path = os.path.join(output_dir, 'Shared_Species_Report.csv')
shared_df.sort_values('Mean_Anthro', ascending=False).to_csv(csv_path, index=False)

# --- Plot F: Shared Species Scatter (Generalists vs Specialists) ---
plt.figure(figsize=(10, 8))
sns.scatterplot(data=shared_df, x='Mean_Anthro', y='Mean_Env', 
                hue='Classification', style='Classification', 
                palette={'Anthro-Biased':'#d62728', 'Env-Biased':'#2ca02c', 'Generalist (Equal)':'gray'},
                s=100, alpha=0.8)

# Equality Line
max_val = max(shared_df['Mean_Anthro'].max(), shared_df['Mean_Env'].max())
plt.plot([0.0001, max_val], [0.0001, max_val], ls='--', c='black', alpha=0.3, label='1:1 Line')

# Label Top Shared
shared_df['Total_Mean'] = shared_df['Mean_Anthro'] + shared_df['Mean_Env']
top_shared = shared_df.nlargest(10, 'Total_Mean')
for _, row in top_shared.iterrows():
    plt.text(row['Mean_Anthro'], row['Mean_Env'], row['Short_Name'], fontsize=9, fontweight='bold')

plt.xscale('log')
plt.yscale('log')
plt.title(f'Shared Species Analysis ({len(shared_species)} species)\nGeneralists vs Specialists')
plt.xlabel('Mean Abundance in Anthropogenic (%)')
plt.ylabel('Mean Abundance in Environmental (%)')
plt.legend(title='Preference')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, '6_Shared_Scatter.png'))
plt.close()

# --- Plot G: Heatmap of Top 20 Shared Species ---
# Select Top 20 most abundant shared species
top_20_ids = top_shared.head(20)['Species'].tolist()
heatmap_data = final_df[top_20_ids + ['Source_Type']].copy()

# Sort by group for visual clarity
heatmap_data = heatmap_data.sort_values('Source_Type')
plot_matrix = heatmap_data[top_20_ids].T 
sample_colors = heatmap_data['Source_Type'].map({'Anthropogenic': '#d62728', 'Environmental': '#2ca02c'})

# Clustermap
g = sns.clustermap(plot_matrix, 
                   col_colors=sample_colors, 
                   col_cluster=False, # Keep sorted by group
                   cmap='viridis', 
                   standard_scale=0, # Normalize rows
                   figsize=(12, 10),
                   cbar_kws={'label': 'Normalized Abundance'})

# Clean Y labels
new_labels = [label.get_text().split('s__')[1].replace('_', ' ') if 's__' in label.get_text() else label.get_text() 
              for label in g.ax_heatmap.get_yticklabels()]
g.ax_heatmap.set_yticklabels(new_labels, rotation=0)
g.ax_heatmap.set_title("Heatmap of Top 20 Shared Species")

plt.savefig(os.path.join(output_dir, '7_Shared_Heatmap.png'))
plt.close()

# --- Plot H: IMPROVED Shared Species Preference (Grouped Comparison) ---
print("Generating Improved Preference Distribution Plot...")

# 1. Підготовка даних (Перетворення в Long Format для Seaborn)
# Ми хочемо, щоб для кожного виду було два рядки: один для Anthro, один для Env
plot_data = shared_df.melt(id_vars=['Short_Name', 'Classification'], 
                           value_vars=['Mean_Anthro', 'Mean_Env'],
                           var_name='Source', value_name='Abundance')

# Перейменуємо для краси
plot_data['Source'] = plot_data['Source'].map({'Mean_Anthro': 'Anthropogenic Abundance', 
                                               'Mean_Env': 'Environmental Abundance'})

# Визначаємо порядок категорій для осі X
order = ['Anthro-Biased', 'Generalist (Equal)', 'Env-Biased']

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7), gridspec_kw={'width_ratios': [1, 2]})

# --- PART 1: Barplot (Кількість видів) ---
# Рахуємо кількість
counts = shared_df['Classification'].value_counts().reindex(order)
colors_map = {'Anthro-Biased': '#d62728', 'Generalist (Equal)': 'gray', 'Env-Biased': '#2ca02c'}

bars = ax1.bar(counts.index, counts.values, color=[colors_map[x] for x in counts.index], 
               edgecolor='black', alpha=0.8)

ax1.set_title('A. Number of Shared Species\nby Category', fontsize=14, fontweight='bold')
ax1.set_ylabel('Count of Species', fontsize=12)
ax1.bar_label(bars, fmt='%d', padding=3, fontsize=12, fontweight='bold')
ax1.set_xticks(range(len(counts)))
ax1.set_xticklabels(['Anthro-\nBiased', 'Generalist', 'Env-\nBiased'], fontsize=11)
ax1.grid(axis='y', linestyle='--', alpha=0.3)

# --- PART 2: Grouped Boxplot (Порівняння чисельності) ---
# Це найважливіша частина: порівнюємо рівні в обох середовищах
sns.boxplot(data=plot_data, x='Classification', y='Abundance', hue='Source', 
            order=order, ax=ax2, palette={'Anthropogenic Abundance': '#d62728', 'Environmental Abundance': '#2ca02c'},
            boxprops=dict(alpha=0.7), flierprops=dict(marker='o', markersize=3, alpha=0.5))

# Додаємо логарифмічну шкалу, бо різниця може бути великою
ax2.set_yscale('log')

ax2.set_title('B. Abundance Shifts: The "Leakage" Effect', fontsize=14, fontweight='bold')
ax2.set_ylabel('Mean Relative Abundance (%) - Log Scale', fontsize=12)
ax2.set_xlabel('Ecological Classification', fontsize=12)

# Переміщуємо легенду
ax2.legend(title='Measured In:', loc='upper right', frameon=True)
ax2.grid(True, which="both", ls="--", alpha=0.2)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, '8_Shared_Species_Preference_IMPROVED.png'), dpi=300)
plt.close()

print(f"ANALYSIS COMPLETE. Results saved to: {output_dir}")