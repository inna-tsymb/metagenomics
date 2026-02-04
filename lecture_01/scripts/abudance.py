import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ==========================================
# 0. CONFIGURATION & SETUP
# ==========================================
# Path to input abundance table and output directory
input_file = '/Users/user/Documents/metagenomics/lecture_01/input/abundance_table.csv'
output_dir = '/Users/user/Documents/metagenomics/lecture_01/output/'

# Create output directory if it doesn't exist
if not os.path.exists(output_dir):
    try:
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    except OSError as e:
        print(f"Error creating directory {output_dir}: {e}")
        exit()

# --- NEW: Create 'abundance' subfolder ---
# Text report, CSVs, AND images will go here
abundance_dir = os.path.join(output_dir, 'abundance')
if not os.path.exists(abundance_dir):
    try:
        os.makedirs(abundance_dir)
        print(f"Created abundance directory: {abundance_dir}")
    except OSError as e:
        print(f"Error creating abundance directory {abundance_dir}: {e}")
# -----------------------------------------

# Set visual style
sns.set_style("whitegrid")
plt.rcParams.update({'font.size': 12, 'figure.figsize': (10, 6)})

# ==========================================
# 1. LOAD & PREPROCESS DATA
# ==========================================
print(f"Loading data from: {input_file}")

try:
    raw_df = pd.read_csv(input_file, index_col=0).fillna(0)
except FileNotFoundError:
    print(f"ERROR: File not found at {input_file}. Please check the path.")
    exit()

# --- Orientation Check ---
# We need Rows = Species, Cols = Samples.
first_index = str(raw_df.index[0])
if "MB-" in first_index or "Sample" in first_index:
    print("  -> Detected Samples as Rows. Transposing to (Species x Samples)...")
    df_species = raw_df.T
else:
    print("  -> Detected Species as Rows. Orientation correct.")
    df_species = raw_df

# --- Data Safety ---
# Force numeric types (coerce errors to NaN then fill with 0) to handle accidental strings
df_species = df_species.apply(pd.to_numeric, errors='coerce').fillna(0)

print(f"  -> Final Dimensions: {df_species.shape[0]} Species x {df_species.shape[1]} Samples")

# ==========================================
# 2. CALCULATE STATISTICS
# ==========================================
print("Calculating statistics...")

# --- A. SAMPLE STATS (Alpha Diversity) ---
# Relative Abundance calculation
df_rel = df_species.div(df_species.sum(axis=0), axis=1)

# Richness (Count > 0)
richness = (df_species > 0).sum(axis=0)

# Shannon Diversity: -sum(p * ln(p))
# Using 1e-9 to avoid log(0) error
shannon = - (df_rel * np.log(df_rel + 1e-9)).sum(axis=0)

# Pielou's Evenness: H / ln(S)
# Handle cases where richness is 0 or 1 to avoid division by zero
evenness = shannon / np.log(richness + 1e-9)
evenness[richness <= 1] = 0

sample_stats = pd.DataFrame({
    'Richness': richness,
    'Shannon': shannon,
    'Evenness': evenness,
    'Total_Reads': df_species.sum(axis=0)
})

# CHANGED: Saving CSV to abundance_dir
sample_csv_path = os.path.join(abundance_dir, 'sample_statistics.csv')
sample_stats.to_csv(sample_csv_path)
print(f"  -> Saved {sample_csv_path}")

# --- B. SPECIES STATS (Prevalence) ---
species_stats = pd.DataFrame({
    'Prevalence': (df_species > 0).sum(axis=1) / df_species.shape[1],
    'Mean_Rel_Abund': df_rel.mean(axis=1)
})
# CHANGED: Saving CSV to abundance_dir
species_csv_path = os.path.join(abundance_dir, 'species_statistics.csv')
species_stats.to_csv(species_csv_path)
print(f"  -> Saved {species_csv_path}")

# --- C. CORE MICROBIOME STATS ---
# Prepare data for Plot 6 (Core Analysis)
df_core = pd.DataFrame({
    'Species': df_species.index,
    'Prevalence': (df_species > 0).sum(axis=1) / df_species.shape[1],
    'Mean_Rel_Abund': df_rel.mean(axis=1)
})

# Identify Top 10 species for labeling
top_species = df_core.nlargest(10, 'Mean_Rel_Abund')

# --- D. CONSOLE TEXT REPORT & FILE SAVE ---
# This section prints to console AND saves to a text file in 'abundance' folder
report_file_path = os.path.join(abundance_dir, 'summary_report.txt')

print(f"\nGenerating text report to: {report_file_path}")

with open(report_file_path, 'w') as f:
    # Helper function to print to both screen and file
    def log(text):
        print(text)
        f.write(text + "\n")

    log("\n" + "="*40)
    log("SUMMARY STATISTICS REPORT")
    log("="*40)
    
    # 1. Counts
    log(f"• Total Samples: {df_species.shape[1]}")
    log(f"• Total Bacterial Species: {df_species.shape[0]}")

    # 2. Top 5 Prevalent
    top_5_prev = species_stats.sort_values('Prevalence', ascending=False).head(5)
    log("\n• Top 5 Most Prevalent Species:")
    for name, row in top_5_prev.iterrows():
        log(f"  - {name}: {row['Prevalence']:.2%} of samples")

    # 3. Top 5 Abundant
    top_5_abund = species_stats.sort_values('Mean_Rel_Abund', ascending=False).head(5)
    log("\n• Top 5 Species by Mean Relative Abundance:")
    for name, row in top_5_abund.iterrows():
        log(f"  - {name}: {row['Mean_Rel_Abund']:.4f} (mean fraction)")
    
    log("="*40 + "\n")

# ==========================================
# 3. GENERATE PLOTS
# ==========================================
print("Generating plots...")

def save_plot(name):
    """Helper to save plots consistently to the ABUNDANCE folder"""
    plt.tight_layout()
    # CHANGED: Now saves to abundance_dir
    save_path = os.path.join(abundance_dir, name)
    plt.savefig(save_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"  -> Saved {name} to abundance/")

# --- Plot 1: Violin Plots (Alpha Diversity) ---
fig, axes = plt.subplots(1, 3, figsize=(15, 6))

sns.violinplot(y=sample_stats['Richness'], ax=axes[0], color='forestgreen', inner='quartile')
axes[0].set_title('Species Richness Distribution')
axes[0].set_ylabel('Number of Observed Species')

sns.violinplot(y=sample_stats['Shannon'], ax=axes[1], color='rebeccapurple', inner='quartile')
axes[1].set_title('Shannon Diversity Distribution')
axes[1].set_ylabel('Shannon Index (H)')

sns.violinplot(y=sample_stats['Evenness'], ax=axes[2], color='darkorange', inner='quartile')
axes[2].set_title('Evenness Distribution')
axes[2].set_ylabel('Pielou\'s Evenness (J)')

save_plot('plot_1_diversity_violins.png')

# --- Plot 2: PCA (Beta Diversity) ---
# Log transform (Samples as Rows for PCA)
log_data = np.log(df_rel.T + 1e-6) 
pca = PCA(n_components=2)
coords = pca.fit_transform(log_data)
var_exp = pca.explained_variance_ratio_

plt.figure(figsize=(10, 8))
scatter = plt.scatter(coords[:, 0], coords[:, 1], 
                      c=sample_stats['Shannon'], cmap='viridis', 
                      edgecolor='k', alpha=0.8, s=60)
plt.colorbar(scatter, label='Shannon Diversity')
plt.xlabel(f'PC1 ({var_exp[0]:.1%} Variance)')
plt.ylabel(f'PC2 ({var_exp[1]:.1%} Variance)')
plt.title('PCA of Microbiome Composition (Beta Diversity)')
plt.grid(True, linestyle='--', alpha=0.5)
save_plot('plot_2_pca_samples.png')

# --- Plot 3: Rank-Abundance Curve ---
plt.figure(figsize=(10, 6))
sorted_abund = species_stats['Mean_Rel_Abund'].sort_values(ascending=False).values
ranks = range(1, len(sorted_abund) + 1)

plt.plot(ranks, sorted_abund, marker='o', linestyle='-', color='black', markersize=4, label='Cohort Mean')

# Add "Zone" Legends
plt.axvspan(0, 5, color='red', alpha=0.1, label='Dominant Species')
plt.axvspan(5, len(ranks), color='blue', alpha=0.1, label='Rare Biosphere')

plt.yscale('log')
plt.xlabel('Species Rank')
plt.ylabel('Mean Relative Abundance (Log Scale)')
plt.title('Rank-Abundance Curve')
plt.legend(loc='upper right')
save_plot('plot_3_rank_abundance.png')

# --- Plot 4: QC Correlation (Reads vs Richness) ---
plt.figure(figsize=(8, 6))
sns.scatterplot(data=sample_stats, x='Total_Reads', y='Richness', 
                hue='Shannon', palette='magma', size='Shannon', sizes=(20, 200),
                edgecolor='black', alpha=0.7)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Shannon Diversity')
plt.title('QC Check: Sequencing Depth vs. Richness')
plt.xlabel('Total Reads per Sample')
plt.ylabel('Observed Species Richness')
plt.grid(True, linestyle='--', alpha=0.5)
save_plot('plot_4_qc_depth_vs_richness.png')

# --- Plot 5: Species Accumulation Curve ---
n_samples = df_species.shape[1]
shuffled_indices = np.random.permutation(n_samples)
shuffled_df = df_species.iloc[:, shuffled_indices]

presence = (shuffled_df > 0).values
accumulator = []
seen_species = np.zeros(presence.shape[0], dtype=bool)

for col in range(n_samples):
    seen_species |= presence[:, col]
    accumulator.append(seen_species.sum())

plt.figure(figsize=(10, 6))
plt.plot(range(1, n_samples + 1), accumulator, color='teal', linewidth=2, label='Observed Species')
plt.title('Species Accumulation Curve (Cohort Sufficiency)')
plt.xlabel('Number of Samples Sequenced')
plt.ylabel('Cumulative Number of Unique Species Found')
plt.legend(loc='lower right')
plt.grid(True, linestyle='--')
save_plot('plot_5_accumulation_curve.png')

# --- Plot 6: Core Microbiome (Prevalence vs Abundance) ---
plt.figure(figsize=(10, 8))
sns.scatterplot(data=df_core, x='Mean_Rel_Abund', y='Prevalence', 
                alpha=0.6, color='grey', edgecolor=None, s=30, label='Other Species')
sns.scatterplot(data=top_species, x='Mean_Rel_Abund', y='Prevalence', 
                color='red', s=100, edgecolor='black', label='Dominant Core Species')

# Add Labels
for line in range(0, top_species.shape[0]):
    plt.text(
        top_species.Mean_Rel_Abund.iloc[line], 
        top_species.Prevalence.iloc[line] + 0.02, 
        top_species.Species.iloc[line], 
        horizontalalignment='left', size='small', color='black', weight='semibold'
    )

plt.xscale('log')
plt.xlabel('Mean Relative Abundance (Log Scale)')
plt.ylabel('Prevalence (Fraction of Samples)')
plt.title('Core Microbiome Analysis: Prevalence vs. Abundance')
plt.axhline(0.5, linestyle='--', color='blue', alpha=0.3, label='50% Prevalence')
plt.legend(loc='lower right')
plt.grid(True, which="both", ls="--", alpha=0.3)
save_plot('plot_6_core_prevalence.png')

# --- Plot 7: Clustered Heatmap (Top 20 Species) ---
print("  -> Generating Clustered Heatmap (Top 20 Species)...")
# Filter for Top 20 most abundant species
top_20_indices = species_stats.sort_values('Mean_Rel_Abund', ascending=False).head(20).index
heatmap_data = df_rel.loc[top_20_indices]

# Draw clustermap
try:
    g = sns.clustermap(heatmap_data, cmap='viridis', standard_scale=None, z_score=0, 
                       figsize=(12, 10), method='ward', metric='euclidean',
                       dendrogram_ratio=(.1, .2), cbar_pos=(0, .2, .03, .4))
    g.fig.suptitle('Clustered Heatmap (Top 20 Species)', fontsize=16, y=1.02)
    # Improve tick label readability
    try:
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=8)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=8)
        g.fig.subplots_adjust(bottom=0.25, left=0.2)
    except Exception:
        pass
    # CHANGED: Saving manually to abundance_dir
    g.savefig(os.path.join(abundance_dir, 'plot_7_clustered_heatmap.png'), bbox_inches='tight')
    plt.close()
    print("  -> Saved plot_7_clustered_heatmap.png to abundance/")
except Exception as e:
    print(f"  -> Skipped Heatmap due to error (possibly too few samples/species): {e}")


# --- Plot 8: Bar Plot of Top Prevalent Species (ADDED TASK) ---
print("  -> Generating Bar Plot (Task 1)...")
top_10_prev = species_stats.sort_values('Prevalence', ascending=False).head(10)
# Convert to a tidy DataFrame so we can assign hue and suppress legend (avoids FutureWarning)
df_top = top_10_prev.reset_index()
# Ensure the first column is named 'Species'
df_top.rename(columns={df_top.columns[0]: 'Species'}, inplace=True)

plt.figure(figsize=(10, 6))
ax = sns.barplot(data=df_top, x='Species', y='Prevalence', hue='Species', palette='viridis', dodge=False)
# Remove the legend since hue is identical to x (same visual effect)
if ax.get_legend() is not None:
    ax.get_legend().remove()

plt.title('Top 10 Most Prevalent Bacterial Species')
plt.ylabel('Prevalence (Fraction of Samples)')
plt.xlabel('Species')
plt.xticks(rotation=45, ha='right')
plt.ylim(0, 1.1)
save_plot('plot_8_bar_prevalence.png')

# --- Plot 9: Scatter Prevalence vs Abundance (ADDED TASK - AXES SWAPPED) ---
print("  -> Generating Scatter Plot (Task 2)...")
plt.figure(figsize=(10, 8))

# Background points
sns.scatterplot(data=species_stats, x='Prevalence', y='Mean_Rel_Abund', 
                alpha=0.6, color='grey', edgecolor='k', s=40)

# Highlight Top 5 Abundant
sns.scatterplot(data=top_5_abund, x='Prevalence', y='Mean_Rel_Abund', 
                color='red', s=100, edgecolor='k', label='Top 5 Abundant')

# Label Top 5
for name, row in top_5_abund.iterrows():
    plt.text(row['Prevalence'] + 0.02, row['Mean_Rel_Abund'], name, 
             fontsize=9, fontweight='bold', color='darkred')

plt.title('Species Prevalence vs. Mean Relative Abundance')
plt.xlabel('Species Prevalence (Fraction of Samples)') # X-axis
plt.ylabel('Mean Relative Abundance (Log Scale)') # Y-axis
plt.yscale('log')
plt.grid(True, which="both", ls="--", alpha=0.3)
plt.legend()
save_plot('plot_9_scatter_prevalence_x_abundance_y.png')

# --- Plot 10: Top 10 Samples by Shannon Diversity ---
print("  -> Generating Bar Plot: Top 10 Samples by Shannon Diversity...")
# Select top 10 samples by Shannon index
top10_shannon = sample_stats.sort_values('Shannon', ascending=False).head(10).reset_index()
# Normalize the first column name to 'SampleID' (works even if index had a name)
if top10_shannon.columns[0] != 'SampleID':
    top10_shannon.rename(columns={top10_shannon.columns[0]: 'SampleID'}, inplace=True)

plt.figure(figsize=(12, 6))
palette = sns.color_palette("magma", n_colors=len(top10_shannon))
ax = sns.barplot(data=top10_shannon, x='SampleID', y='Shannon', palette=palette, hue='SampleID', dodge=False)
# Remove legend since hue duplicates x
if ax.get_legend() is not None:
    ax.get_legend().remove()

plt.title('Top 10 Samples by Shannon Diversity')
plt.ylabel('Shannon Diversity Index')
plt.xlabel('Sample ID')
plt.xticks(rotation=45, ha='right')
# Add a small top margin
plt.ylim(0, top10_shannon['Shannon'].max() + 0.2)
save_plot('plot_10_top10_shannon.png')

# --- Plot 11: PCA of Species (Annotated, numbered legend) ---
print("  -> Generating PCA Plot of Species (numbered legend for clarity)...")
from sklearn.decomposition import PCA as sklearnPCA

# PCA on species (rows = species x samples)
X_species = df_rel.fillna(0).values
pca_sp = sklearnPCA(n_components=2)
coords_sp = pca_sp.fit_transform(X_species)
var_sp = pca_sp.explained_variance_ratio_

plt.figure(figsize=(12, 10))
plt.scatter(coords_sp[:, 0], coords_sp[:, 1], s=40, color='navy', alpha=0.8, edgecolor='k')

# Helper to produce short, readable labels
def short_label(s, maxlen=30):
    s = str(s).replace('_', ' ')
    return s if len(s) <= maxlen else s[:maxlen-3] + '...'

# Annotate only the top N species by mean abundance and use numbered markers
top_n = 5
top_species = species_stats.sort_values('Mean_Rel_Abund', ascending=False).head(top_n).index.tolist()
for idx, sp in enumerate(top_species, start=1):
    i = list(df_rel.index).index(sp)
    x, y = coords_sp[i, 0], coords_sp[i, 1]
    # Draw a larger numbered marker
    plt.scatter([x], [y], s=150, color='crimson', edgecolor='k', zorder=4)
    plt.text(x, y, str(idx), fontsize=9, fontweight='bold', color='white', ha='center', va='center', zorder=5)

# Add a boxed legend to the right listing the full (or shortened) species names
legend_lines = [f"{i+1}. {short_label(sp, maxlen=40)}" for i, sp in enumerate(top_species)]
fig = plt.gcf()
fig.subplots_adjust(right=0.78)
fig.text(0.82, 0.5, "\n".join(legend_lines), fontsize=9, ha='left', va='center',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.95, edgecolor='k'))

plt.xlabel(f'Principal Component 1 ({var_sp[0]:.1%})')
plt.ylabel(f'Principal Component 2 ({var_sp[1]:.1%})')
plt.title('PCA of Microbiome Data (Species-level)')
plt.grid(True, linestyle='--', alpha=0.4)
save_plot('plot_11_pca_species_annotated.png')

# Plot 12 removed: annotated QC scatter was removed to reduce label clutter. Consider using the existing Plot 4 (QC Correlation) for depth vs richness checks.

# --- Plot 13: Species Prevalence vs Mean Relative Abundance (Bubble with legend) ---
print("  -> Generating Enhanced Prevalence vs Abundance Bubble Plot (numbered labels)...")
plt.figure(figsize=(12, 8))
# Normalize sizes for visualization
size_norm = (species_stats['Mean_Rel_Abund'] / species_stats['Mean_Rel_Abund'].max()) * 600
sc = plt.scatter(species_stats['Prevalence'], species_stats['Mean_Rel_Abund'],
                 s=size_norm, c=species_stats['Mean_Rel_Abund'], cmap='viridis', alpha=0.8, edgecolor='k')
plt.colorbar(sc, label='Mean Relative Abundance')

# Annotate only the top N species by mean abundance with numbered markers
top_n = 6
top6_abund = species_stats.sort_values('Mean_Rel_Abund', ascending=False).head(top_n)
for idx, (name, row) in enumerate(top6_abund.iterrows(), start=1):
    x = row['Prevalence']
    y = row['Mean_Rel_Abund']
    plt.scatter([x], [y], s=220, color='crimson', edgecolor='k', zorder=4)
    plt.text(x, y, str(idx), fontsize=9, fontweight='bold', color='white', ha='center', va='center', zorder=5)

# Add boxed legend to the right showing which number = species name
# Reuse short_label helper to keep names concise
def short_label(s, maxlen=40):
    s = str(s).replace('_', ' ')
    return s if len(s) <= maxlen else s[:maxlen-3] + '...'

legend_lines = [f"{i+1}. {short_label(name)}" for i, name in enumerate(top6_abund.index)]
fig = plt.gcf()
fig.subplots_adjust(right=0.78)
fig.text(0.82, 0.5, "\n".join(legend_lines), fontsize=9, ha='left', va='center',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.95, edgecolor='k'))

plt.title('Species: Prevalence vs. Mean Relative Abundance')
plt.xlabel('Prevalence')
plt.ylabel('Mean Relative Abundance')
plt.grid(True, linestyle='--', alpha=0.3)
save_plot('plot_13_species_prevalence_vs_abundance.png')

# --- Plot 14: Phylum Core Analysis (requires taxonomy file) ---
print("  -> Generating Phylum Core Analysis (if taxonomy available)...")
try:
    taxonomy_path = os.path.join(os.path.dirname(input_file), 'mag_data_taxa.xlsx')
    if os.path.exists(taxonomy_path):
        taxa_df = pd.read_excel(taxonomy_path, index_col=0)
        taxa_df.columns = [c.capitalize() for c in taxa_df.columns]
        # Find species column
        species_col = None
        for c in taxa_df.columns:
            if c.lower() in ('sp', 'species'):
                species_col = c
                break
        if species_col is None:
            print('  -> Taxonomy file does not contain species column. Skipping phylum plot.')
        else:
            # Build lookup
            def clean_name(name):
                return str(name).lower().replace('_', ' ').strip()
            taxa_by_sp = taxa_df[~taxa_df[species_col].isna()].drop_duplicates(subset=species_col).set_index(species_col)
            taxa_lookup = {clean_name(idx): taxa_by_sp.loc[idx] for idx in taxa_by_sp.index}

            phylum_map = {}
            import difflib
            for sp in df_rel.index:
                csp = clean_name(sp)
                if csp in taxa_lookup:
                    ph = taxa_lookup[csp].get('Phylum', 'Unknown')
                else:
                    close = difflib.get_close_matches(csp, taxa_lookup.keys(), n=1, cutoff=0.7)
                    if close:
                        ph = taxa_lookup[close[0]].get('Phylum', 'Unknown')
                    else:
                        ph = 'Unknown'
                phylum_map[sp] = ph
            ph_series = pd.Series(phylum_map)
            # Per-sample phylum abundance
            phylum_abund = df_rel.groupby(ph_series).sum()
            ph_prevalence = (phylum_abund > 0).sum(axis=1) / df_rel.shape[1]
            ph_mean_abund = phylum_abund.mean(axis=1)

            plt.figure(figsize=(10, 8))
            sizes = (ph_mean_abund / ph_mean_abund.max()) * 2000
            sc = plt.scatter(ph_prevalence, ph_mean_abund, s=sizes, c=ph_mean_abund, cmap='plasma', alpha=0.9, edgecolor='k')
            plt.colorbar(sc, label='Mean Relative Abundance')
            for ph in ph_prevalence.index:
                plt.text(ph_prevalence[ph] + 0.01, ph_mean_abund[ph], ph, fontsize=10, weight='semibold')
            plt.xlabel('Prevalence (Fraction of Samples)')
            plt.ylabel('Mean Relative Abundance')
            plt.title('Phylum Core Analysis: Prevalence vs Abundance')
            plt.grid(True, linestyle='--', alpha=0.3)
            save_plot('plot_14_phylum_core_analysis.png')
    else:
        print('  -> Taxonomy file not found; skipping Phylum Core Analysis plot.')
except Exception as e:
    print(f'  -> Error generating phylum plot: {e}')

print(f"\nAll processing complete. All files (CSV, TXT, PNG) are in: {abundance_dir}")