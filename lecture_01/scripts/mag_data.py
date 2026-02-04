import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
import os
import difflib

# ==========================================
# 0. CONFIGURATION & SETUP
# ==========================================
input_dir = '/Users/user/Documents/metagenomics/lecture_01/input/'
base_output_dir = '/Users/user/Documents/metagenomics/lecture_01/output/'

# --- Create 'taxonomy' subfolder ---
output_dir = os.path.join(base_output_dir, 'taxonomy')
if not os.path.exists(output_dir):
    try:
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")
    except OSError as e:
        print(f"Error creating directory: {e}")

# Files
abundance_file = os.path.join(input_dir, 'abundance_table.csv')
taxonomy_file = os.path.join(input_dir, 'mag_data_taxa.xlsx')

# Visual Setup
sns.set_style("whitegrid")
plt.rcParams.update({'font.size': 12, 'figure.figsize': (10, 6)})

def log_to_file(message):
    logfile = os.path.join(output_dir, 'taxonomy_report.txt')
    with open(logfile, 'a') as f:
        f.write(message + "\n")
    print(message)

if os.path.exists(os.path.join(output_dir, 'taxonomy_report.txt')):
    os.remove(os.path.join(output_dir, 'taxonomy_report.txt'))

log_to_file("=== TAXONOMIC TABLE ANALYSIS REPORT ===\n")

# ==========================================
# 1. LOAD DATA
# ==========================================
print("Loading data...")

# Load Abundance
if os.path.exists(abundance_file):
    abundance_table = pd.read_csv(abundance_file, index_col=0).fillna(0)
else:
    print(f"Error: {abundance_file} not found.")
    exit()

# Load Taxonomy
if os.path.exists(taxonomy_file):
    taxa_df = pd.read_excel(taxonomy_file, index_col=0)
else:
    print(f"Error: {taxonomy_file} not found. Taxonomy analysis requires this file.")
    exit()

# Auto-Orient Abundance (make rows = Species, cols = Samples)
# If index looks like sample IDs (e.g., 'MB-') and columns look like species names
idx_is_samples = abundance_table.index.astype(str).str.match(r'^MB-').any()
cols_have_spaces = abundance_table.columns.astype(str).str.contains(' ').any()
if idx_is_samples and cols_have_spaces:
    abundance_table = abundance_table.T

# Calculate Relative Abundance (rows=Species, cols=Samples)
df_rel = abundance_table.div(abundance_table.sum(axis=0), axis=1)
log_to_file(f"Loaded {df_rel.shape[0]} species and {df_rel.shape[1]} samples.")

# ==========================================
# 2. ROBUST MAPPING (Species -> Taxonomy)
# ==========================================
log_to_file("\n--- Mapping Species to Taxonomy ---")

# Standardize column names in taxonomy file (Title Case)
taxa_df.columns = [c.capitalize() for c in taxa_df.columns]
required_levels = ['Phylum', 'Class', 'Order', 'Family', 'Genus']

# Create mapping dataframe
mapped_df = pd.DataFrame(index=df_rel.index)
for level in required_levels:
    mapped_df[level] = 'Unknown'

# Fuzzy Matching Function
def clean_name(name):
    return str(name).lower().replace('_', ' ').strip()

# Determine which column contains species names (accept 'sp' or 'species')
species_col = None
for c in taxa_df.columns:
    if c.lower() in ('sp', 'species'):
        species_col = c
        break

if species_col is None:
    log_to_file("ERROR: Taxonomy file does not contain a 'sp' or 'species' column. Cannot map species.")
    mapped_count = 0
else:
    # Build taxa lookup from the species column (deduplicate by species name)
    taxa_by_sp = taxa_df[~taxa_df[species_col].isna()].drop_duplicates(subset=species_col).set_index(species_col)
    taxa_lookup = {clean_name(idx): taxa_by_sp.loc[idx] for idx in taxa_by_sp.index}

    mapped_count = 0
    unmatched_examples = []
    for species in df_rel.index:
        clean_sp = clean_name(species)

        # Try Exact or Fuzzy Match
        match = None
        if clean_sp in taxa_lookup:
            match = taxa_lookup[clean_sp]
        else:
            # Fuzzy match
            close = difflib.get_close_matches(clean_sp, taxa_lookup.keys(), n=1, cutoff=0.7)
            if close:
                match = taxa_lookup[close[0]]

        if match is not None:
            mapped_count += 1
            for level in required_levels:
                if level in match:
                    mapped_df.loc[species, level] = match[level]
        else:
            if len(unmatched_examples) < 10:
                unmatched_examples.append(species)

    if len(unmatched_examples) > 0:
        log_to_file(f"Note: Examples of unmapped species: {unmatched_examples}")

log_to_file(f"Successfully mapped {mapped_count} out of {len(df_rel.index)} species.")

# Combine Abundance with Taxonomy
# Calculate Mean Relative Abundance across all samples for the 'Average Human' view
mean_abundance = df_rel.mean(axis=1)
full_data = mapped_df.copy()
full_data['Mean_Abundance'] = mean_abundance

# ==========================================
# 3. AGGREGATION & SUMMARY STATISTICS
# ==========================================
log_to_file("\n--- Phylum Composition Statistics ---")

# Group by Phylum
phylum_stats = full_data.groupby('Phylum')['Mean_Abundance'].sum().sort_values(ascending=False)
log_to_file("Top Phyla by Mean Relative Abundance:")
log_to_file(phylum_stats.to_string())

# Save aggregated tables
for level in required_levels:
    agg_df = full_data.groupby(level)['Mean_Abundance'].sum().sort_values(ascending=False)
    agg_df.to_csv(os.path.join(output_dir, f'stats_by_{level}.csv'))

# ==========================================
# 4. VISUALIZATIONS
# ==========================================
print("Generating visualizations...")

# --- A. STACKED BAR PLOT (Phylum Level) ---
# We need per-sample data for this, not mean
df_with_taxa = df_rel.join(mapped_df[['Phylum']])
phylum_sample_data = df_with_taxa.groupby('Phylum').sum()

# Filter small phyla for cleaner plot (Top 7 + Others)
top_7_phyla = phylum_sample_data.sum(axis=1).nlargest(7).index
phylum_sample_data_plot = phylum_sample_data.loc[top_7_phyla]
# Calculate 'Others'
others = phylum_sample_data.loc[~phylum_sample_data.index.isin(top_7_phyla)].sum(axis=0)
if others.sum() > 0:
    phylum_sample_data_plot.loc['Others'] = others

plt.figure(figsize=(12, 6))
phylum_sample_data_plot.T.plot(kind='bar', stacked=True, colormap='Paired', width=1.0, figsize=(14, 7))
plt.title('Taxonomic Composition by Phylum (Top 7)', fontsize=15)
plt.ylabel('Relative Abundance')
plt.xlabel('Samples (x-axis labels hidden)')
plt.xticks([]) # Hide dense sample labels
plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', title='Phylum')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, '1_stacked_bar_phylum.png'))
plt.close()

# --- B. SUNBURST PLOT (Hierarchical) ---
# Using the mean abundance data (The "Average" Cohort Profile)
# Reset index to make 'Species' a column
viz_data = full_data.reset_index().rename(columns={'index': 'Species'})
# Filter out zero abundance rows just in case
viz_data = viz_data[viz_data['Mean_Abundance'] > 0]

fig_sun = px.sunburst(
    viz_data,
    path=['Phylum', 'Family', 'Genus'], # Hierarchy path
    values='Mean_Abundance',
    color='Phylum', # Color segments by Phylum
    title='Hierarchical Microbiome Composition (Average Cohort)',
    width=800, height=800
)
fig_sun.write_html(os.path.join(output_dir, '2_sunburst_hierarchy.html'))
# Save static image too (requires kaleido, skipping if not installed, saving html is safer)
try:
    fig_sun.write_image(os.path.join(output_dir, '2_sunburst_hierarchy.png'))
except:
    print("  -> Could not save static Sunburst (install 'kaleido'). HTML saved.")

# --- C. SANKEY DIAGRAM (Phylum -> Family) ---
# Logic: Link Phylum Source to Family Target with width = Abundance
print("  -> Generating Sankey Diagram...")

# 1. Aggregate data: Phylum -> Family
sankey_data = viz_data.groupby(['Phylum', 'Family'])['Mean_Abundance'].sum().reset_index()

# 2. Create unique labels for all nodes
# Note: A family "Unknown" and a Phylum "Unknown" must be distinct nodes
sankey_data['Phylum_Label'] = sankey_data['Phylum'] + " (P)"
sankey_data['Family_Label'] = sankey_data['Family'] + " (F)"

all_nodes = list(pd.concat([sankey_data['Phylum_Label'], sankey_data['Family_Label']]).unique())
node_map = {node: i for i, node in enumerate(all_nodes)}

# 3. Create Source, Target, Value lists
sources = [node_map[p] for p in sankey_data['Phylum_Label']]
targets = [node_map[f] for f in sankey_data['Family_Label']]
values = sankey_data['Mean_Abundance'].values

# 4. Plot
fig_sankey = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=[node.replace(" (P)", "").replace(" (F)", "") for node in all_nodes], # Clean labels for display
        color="blue"
    ),
    link=dict(
        source=sources,
        target=targets,
        value=values
    )
)])

fig_sankey.update_layout(title_text="Taxonomic Flow: Phylum -> Family", font_size=12)
fig_sankey.write_html(os.path.join(output_dir, '3_sankey_phylum_family.html'))
try:
    fig_sankey.write_image(os.path.join(output_dir, '3_sankey_phylum_family.png'))
except:
    pass

print(f"\nProcessing complete. Check {output_dir}")