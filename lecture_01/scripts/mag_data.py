import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import os
import re
import difflib

# ==========================================
# 0. CONFIGURATION & SETUP
# ==========================================
input_dir = '/Users/user/Documents/metagenomics/lecture_01/input/'
output_dir = '/Users/user/Documents/metagenomics/lecture_01/output/'

# Define specific file paths
abundance_path = os.path.join(input_dir, 'abundance_table.csv')
taxonomy_path = '/Users/user/Documents/metagenomics/lecture_01/input/mag_data_taxa.xlsx'
stats_file = os.path.join(output_dir, 'taxonomic_statistics.txt')

# Ensure output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Matching & mapping configuration
FUZZY_CUTOFF = 0.7  # Lower cutoff to be more permissive (0..1)
MAPPING_OVERRIDE_FILE = os.path.join(input_dir, 'species_mapping_overrides.csv')  # optional: abundance_species,mapped_species
MAPPING_REPORT_FILE = os.path.join(output_dir, 'species_mapping_report.csv')

# Initialize Statistics File
with open(stats_file, 'w') as f:
    f.write("=== TAXONOMIC ANALYSIS REPORT ===\n")
    f.write(f"Generated on: {pd.Timestamp.now()}\n\n")

def log_stat(title, content):
    """Helper to write to text file and print to console"""
    with open(stats_file, 'a') as f:
        f.write(f"--- {title} ---\n")
        f.write(str(content))
        f.write("\n\n")
    print(f"Processed: {title}")

# Set Visual Style
sns.set_style("whitegrid")
plt.rcParams.update({'font.size': 12, 'figure.figsize': (12, 8)})

# ==========================================
# 1. DATA LOADING
# ==========================================
print("Loading data...")

# Load Abundance (CSV) -> samples x species (rows: samples, cols: species)
# index_col=0 assumes the first column is the sample ID (e.g., 'MB-2860')
abund_df = pd.read_csv(abundance_path, index_col=0).fillna(0)

# Load Taxonomy (Excel)
# Ensure you have openpyxl installed: pip install openpyxl
taxa_df = pd.read_excel(taxonomy_path, index_col=0)

# Map taxonomy by species name (taxonomy file uses 'sp' column)
species_col = 'sp' if 'sp' in taxa_df.columns else ('Species' if 'Species' in taxa_df.columns else None)
if species_col is None:
    log_stat('ERROR', "Taxonomy file missing 'sp' or 'Species' column.")
    print("CRITICAL ERROR: Taxonomy file missing species column ('sp' or 'Species').")
    exit()

# Create taxa mapping with one row per species (taxonomy is consistent across MAGs)
taxa_by_sp = taxa_df.drop_duplicates(subset=species_col).set_index(species_col)
# Rename taxonomy columns to Title case expected downstream
taxa_by_sp = taxa_by_sp.rename(columns={
    'phylum': 'Phylum',
    'family': 'Family',
    'genus': 'Genus',
    'species': 'Species'
})

# Transpose abundance to species x samples (index: species)
species_abund = abund_df.T

# --- Name normalization, fuzzy matching, overrides & mapping report ---
from difflib import SequenceMatcher

def normalize_name(s):
    """Make names comparable: remove underscores/dashes, collapse whitespace, lowercase."""
    s = str(s)
    s = s.replace('_', ' ').replace('-', ' ')
    s = re.sub(r"\s+", ' ', s)
    return s.strip().lower()

# Normalized indexes
taxa_norm = taxa_by_sp.index.to_series().apply(normalize_name)
taxa_norm_to_orig = {}
for orig, norm in taxa_norm.items():
    taxa_norm_to_orig.setdefault(norm, []).append(orig)

species_norm = species_abund.index.to_series().apply(normalize_name)

# Build mapping report rows
rows = []
available_norms = list(taxa_norm_to_orig.keys())
for orig, norm in species_norm.items():
    row = {'abundance_species': orig, 'normalized_abundance': norm, 'matched_taxa': None, 'matched_norm': None, 'match_method': None, 'match_score': None}
    # Exact normalized match
    if norm in taxa_norm_to_orig:
        matched_orig = taxa_norm_to_orig[norm][0]
        row.update({'matched_taxa': matched_orig, 'matched_norm': norm, 'match_method': 'exact', 'match_score': 1.0})
    else:
        # Fuzzy match (top 3 candidates), choose best by SequenceMatcher ratio
        candidates = difflib.get_close_matches(norm, available_norms, n=3, cutoff=FUZZY_CUTOFF)
        best = None
        best_score = 0.0
        best_norm = None
        for cand in candidates:
            score = SequenceMatcher(None, norm, cand).ratio()
            if score > best_score:
                best_score = score
                best = taxa_norm_to_orig[cand][0]
                best_norm = cand
        if best is not None:
            row.update({'matched_taxa': best, 'matched_norm': best_norm, 'match_method': 'fuzzy', 'match_score': float(best_score)})
    rows.append(row)

mapping_df = pd.DataFrame(rows)
# If override file exists, apply overrides
if os.path.exists(MAPPING_OVERRIDE_FILE):
    try:
        overrides = pd.read_csv(MAPPING_OVERRIDE_FILE)
        # Expect columns: abundance_species,mapped_species
        overrides = overrides[['abundance_species', 'mapped_species']].dropna()
        overrides_map = dict(zip(overrides['abundance_species'], overrides['mapped_species']))
        applied = 0
        for idx, r in mapping_df.iterrows():
            a = r['abundance_species']
            if a in overrides_map:
                mapped = overrides_map[a]
                # Validate mapped exists in taxonomy
                if mapped in taxa_by_sp.index:
                    mapping_df.at[idx, 'matched_taxa'] = mapped
                    mapping_df.at[idx, 'matched_norm'] = normalize_name(mapped)
                    mapping_df.at[idx, 'match_method'] = 'override'
                    mapping_df.at[idx, 'match_score'] = 1.0
                    applied += 1
                else:
                    log_stat('Mapping Override Warning', f"Override for {a} maps to unknown taxa species '{mapped}'; ignored.")
        log_stat('Mapping Overrides', f"Overrides file applied. Overrides used: {applied}")
    except Exception as e:
        log_stat('Mapping Overrides Error', str(e))

# Save mapping report
mapping_df.to_csv(MAPPING_REPORT_FILE, index=False)
log_stat('Species Mapping Report', f"Saved mapping report to {MAPPING_REPORT_FILE}")

# Build final mapping dict (only rows with matched_taxa)
valid_map = mapping_df.dropna(subset=['matched_taxa']).set_index('abundance_species')['matched_taxa'].to_dict()
if len(valid_map) == 0:
    print("CRITICAL ERROR: No mapping found after normalization/fuzzy matching/overrides. Check species names or provide overrides.")
    exit()

# Build matched abundance dataframe, rename indices to taxa species names and aggregate duplicates
present_keys = [k for k in valid_map.keys() if k in species_abund.index]
species_abund_matched = species_abund.loc[present_keys].copy()
species_abund_matched.index = [valid_map[o] for o in species_abund_matched.index]
# If multiple original species map to the same taxa species, sum their counts
species_abund_matched = species_abund_matched.groupby(species_abund_matched.index).sum()

# Keep taxonomy rows for matched species (ensure order)
# Some mapped taxa may not be present in taxa_by_sp index (shouldn't happen but guard)
taxa_by_sp = taxa_by_sp.loc[species_abund_matched.index.intersection(taxa_by_sp.index)]

# Combine on species index
combined = species_abund_matched.join(taxa_by_sp, how='inner')
log_stat("Data Dimensions", f"Abundance (species x samples): {species_abund_matched.shape}\nTaxonomy: {taxa_by_sp.shape}\nMerged: {combined.shape}")

# ==========================================
# 2. DATA PROCESSING & STATISTICS
# ==========================================

# A. Phylum Level Statistics
# Group by Phylum and sum abundance across all samples
if 'Phylum' in combined.columns:
    phylum_counts = combined.groupby('Phylum')[species_abund.columns].sum()
    
    # Calculate Relative Abundance (Global)
    total_reads = phylum_counts.sum().sum()
    phylum_rel = phylum_counts.sum(axis=1) / total_reads
    
    # Sort and Log
    top_phyla = phylum_rel.sort_values(ascending=False)
    log_stat("Top Phyla (Global Abundance)", top_phyla.apply(lambda x: f"{x:.2%}"))
else:
    print("Error: 'Phylum' column missing.")
    exit()

# B. Family & Genus Level Statistics
# Identify the most dominant Genera
if 'Genus' in combined.columns:
    genus_counts = combined.groupby('Genus')[species_abund.columns].sum()
    genus_rel = genus_counts.div(genus_counts.sum(axis=0), axis=1) # Normalize per sample
    
    # Get Top 20 Genera by mean abundance
    top_20_genera = genus_rel.mean(axis=1).sort_values(ascending=False).head(20)
    log_stat("Top 20 Genera (Mean Relative Abundance)", top_20_genera.apply(lambda x: f"{x:.4f}"))

# ==========================================
# 3. ADVANCED VISUALIZATIONS
# ==========================================

# --- Plot 1: Stacked Bar Plot (Phylum Composition) ---
plt.figure(figsize=(12, 6))
# Calculate relative abundance per sample
phylum_sample_rel = phylum_counts.div(phylum_counts.sum(axis=0), axis=1)
# Sort to put largest phyla at bottom
sort_order = phylum_rel.sort_values().index
phylum_sample_rel.loc[sort_order].T.plot(kind='bar', stacked=True, colormap='tab20', width=0.9, figsize=(14, 7))

plt.title('Microbial Composition by Phylum', fontsize=16)
plt.ylabel('Relative Abundance', fontsize=14)
plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', title='Phylum')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'taxa_1_stacked_bar_phylum.png'))
plt.close()


# --- Plot 2: Sunburst Plot (Interactive Hierarchy) ---
# Requires Phylum -> Family -> Genus columns
if {'Family', 'Genus'}.issubset(combined.columns):
    # Melt/Aggregate for Plotly
    df_hier = combined.groupby(['Phylum', 'Family', 'Genus'])[species_abund.columns].sum().sum(axis=1).reset_index()
    df_hier.columns = ['Phylum', 'Family', 'Genus', 'Abundance']
    
    fig = px.sunburst(
        df_hier, 
        path=['Phylum', 'Family', 'Genus'], 
        values='Abundance', 
        color='Phylum',
        title='Taxonomic Hierarchy & Abundance',
        color_discrete_sequence=px.colors.qualitative.Pastel
    )
    fig.write_html(os.path.join(output_dir, 'taxa_2_sunburst.html'))
    log_stat("Sunburst Plot", "Saved interactive HTML to taxa_2_sunburst.html")


# --- Plot 3: Clustered Heatmap (Top 20 Genera) ---
# Use the genus_rel data calculated earlier
top_genera_names = top_20_genera.index
heatmap_data = genus_rel.loc[top_genera_names]

# z_score=0 standardizes rows (shows if abundance is higher/lower than average for that genus)
g = sns.clustermap(heatmap_data, cmap='viridis', z_score=0, 
                   method='ward', metric='euclidean',
                   figsize=(12, 10), dendrogram_ratio=(.1, .2))
g.fig.suptitle('Clustered Heatmap (Top 20 Genera)', fontsize=16, y=1.02)
g.savefig(os.path.join(output_dir, 'taxa_3_heatmap_genus.png'))
plt.close()


# --- Plot 4: Firmicutes/Bacteroidetes (F/B) Ratio ---
if 'Firmicutes' in phylum_counts.index and 'Bacteroidetes' in phylum_counts.index:
    firmicutes = phylum_counts.loc['Firmicutes']
    bacteroidetes = phylum_counts.loc['Bacteroidetes']
    
    # Calculate Ratio (add epsilon to avoid div/0)
    fb_ratio = firmicutes / (bacteroidetes + 1e-9)
    
    # Save Stats
    log_stat("F/B Ratio Statistics", fb_ratio.describe())
    
    # Plot Distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(fb_ratio, bins=20, kde=True, color='purple')
    plt.axvline(fb_ratio.mean(), color='red', linestyle='--', label=f'Mean: {fb_ratio.mean():.2f}')
    plt.title('Distribution of Firmicutes/Bacteroidetes Ratio')
    plt.xlabel('F/B Ratio')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'taxa_4_fb_ratio.png'))
    plt.close()


# --- Plot 5: Core Phylum Analysis (Prevalence vs Abundance) ---
# Prevalence = Fraction of samples where count > 0
phylum_prev = (phylum_counts > 0).sum(axis=1) / phylum_counts.shape[1]
phylum_mean = phylum_sample_rel.mean(axis=1)

plt.figure(figsize=(10, 8))
sns.scatterplot(x=phylum_prev, y=phylum_mean, hue=phylum_mean, 
                palette='magma', size=phylum_mean, sizes=(100, 2000), legend=False, alpha=0.7)

# Label points
for name in phylum_prev.index:
    # Only label if prevalence > 0.1 or abundance > 0.01 to avoid clutter
    if phylum_prev[name] > 0.1 or phylum_mean[name] > 0.01:
        plt.text(phylum_prev[name], phylum_mean[name], name, fontsize=10, fontweight='bold', ha='right')

plt.title('Phylum Core Analysis: Prevalence vs Abundance', fontsize=16)
plt.xlabel('Prevalence (Fraction of Samples)', fontsize=14)
plt.ylabel('Mean Relative Abundance', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.5)
plt.savefig(os.path.join(output_dir, 'taxa_5_core_scatter.png'))
plt.close()

print(f"\nAnalysis Complete. All Statistics saved to: {stats_file}")
print(f"All Plots saved to: {output_dir}")