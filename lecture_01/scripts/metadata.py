import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import datetime

# --- Step 0: Configuration ---
output_dir = '/Users/user/Documents/metagenomics/lecture_01/output/'
input_dir = '/Users/user/Documents/metagenomics/lecture_01/input/'
metadata_file = os.path.join(input_dir, 'metadata_table.csv')
abundance_file = os.path.join(input_dir, 'abundance_table.csv')

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Set visual style
sns.set_style("whitegrid")
plt.rcParams.update({'font.size': 12, 'figure.figsize': (10, 6)})


def log_to_file(message, logfile=None):
    """Append a timestamped message to the analysis log and print to stdout."""
    if logfile is None:
        logfile = os.path.join(output_dir, 'analysis.log')
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    line = f"{timestamp} - {message}\n"
    with open(logfile, 'a') as f:
        f.write(line)
    print(message)

# --- Step 1: Data Loading & Generation ---

# 1. Load Abundance Table (Needed for the last part of your request)
# If file exists, load it. If not, lift error.
if os.path.exists(abundance_file):
    abundance_table = pd.read_csv(abundance_file, index_col=0)
else:
    raise FileNotFoundError(f"Abundance table file not found at {abundance_file}")

# 2. Load Metadata
# If file exists, load it. If not, lift error.
if os.path.exists(metadata_file):
    metadata = pd.read_csv(metadata_file, index_col=0)
else:
    raise FileNotFoundError(f"Metadata file not found at {metadata_file}")

# Ensure 'Sex' column exists (compatibility with older naming using 'Gender')
if 'Sex' not in metadata.columns and 'Gender' in metadata.columns:
    # If Gender stored as integers (0/1), map to 'M'/'F' (default mapping: 1 -> 'F', 0 -> 'M')
    if pd.api.types.is_integer_dtype(metadata['Gender']) or pd.api.types.is_float_dtype(metadata['Gender']):
        metadata['Sex'] = metadata['Gender'].map({1: 'F', 0: 'M'})
    else:
        metadata['Sex'] = metadata['Gender'].astype(str)

# Align Data
common_samples = abundance_table.columns.intersection(metadata.index)
# If no overlap, try transposing abundance table (samples may be rows instead of columns)
if len(common_samples) == 0:
    abundance_table = abundance_table.T
    common_samples = abundance_table.columns.intersection(metadata.index)

if len(common_samples) == 0:
    log_to_file("ERROR: No overlapping sample IDs between abundance table and metadata. Check sample naming/orientation.")
    raise ValueError("No overlapping sample IDs between abundance table and metadata.")

abundance_table = abundance_table[common_samples]
metadata = metadata.loc[common_samples]
df_rel_abund = abundance_table.div(abundance_table.sum(axis=0), axis=1)

log_to_file(f"Data Loaded. N Samples: {len(common_samples)}")
log_to_file("-" * 40)

# --- ANALYSIS 1: Diversity vs Sex ---
log_to_file("\n--- 1. Alpha Diversity by Sex ---")

# Calculate Shannon Diversity
shannon = - (df_rel_abund * np.log(df_rel_abund + 1e-9)).sum(axis=0)
metadata['Shannon'] = shannon

# Stats Test
f_div = metadata[metadata['Sex']=='F']['Shannon']
m_div = metadata[metadata['Sex']=='M']['Shannon']
stat, p_div = stats.mannwhitneyu(f_div, m_div)

log_to_file(f"Female Mean Shannon: {f_div.mean():.2f}")
log_to_file(f"Male Mean Shannon:   {m_div.mean():.2f}")
log_to_file(f"Mann-Whitney U Test: p-value = {p_div:.5f}")
if p_div < 0.05:
    log_to_file("RESULT: Significant difference in diversity between sexes.")
else:
    log_to_file("RESULT: No significant difference in diversity.")

# Plot
plt.figure(figsize=(6, 6))
sns.boxplot(x='Sex', y='Shannon', data=metadata, palette='Set2')
sns.stripplot(x='Sex', y='Shannon', data=metadata, color='black', alpha=0.3)
plt.title(f'Shannon Diversity by Sex\n(p={p_div:.4f})')
plt.savefig(os.path.join(output_dir, 'diversity_vs_sex.png'))
plt.close()

# --- ANALYSIS 2: Correlation (Age vs Species) ---
log_to_file("\n--- 2. Species correlation with Age ---")

# We only correlate the top 10 most abundant species to reduce noise
top_10 = df_rel_abund.mean(axis=1).sort_values(ascending=False).head(10).index
correlations = []

for species in top_10:
    # Spearman correlation (non-parametric, safer for biology)
    corr, p_val = stats.spearmanr(df_rel_abund.loc[species], metadata['Age'])
    correlations.append({'Species': species, 'Rho': corr, 'P_value': p_val})

corr_df = pd.DataFrame(correlations).sort_values('P_value')
log_to_file(corr_df.to_string(index=False))

# Plot Heatmap of these correlations
plt.figure(figsize=(8, 6))
# Create a matrix for heatmap (Species x 1 column for Age correlation)
heatmap_data = corr_df.set_index('Species')[['Rho']]
sns.heatmap(heatmap_data, annot=True, cmap='coolwarm', center=0, vmin=-1, vmax=1)
plt.title('Correlation: Top 10 Species vs Age')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'correlation_species_age.png'))
plt.close()

# --- ANALYSIS 3: Core Microbiome ---
log_to_file("\n--- 3. Core Microbiome Analysis ---")
# Definition: Species present in > 50% of samples
prevalence = (abundance_table > 0).sum(axis=1) / len(common_samples)
core_species = prevalence[prevalence > 0.5].index.tolist()

log_to_file(f"Total Species: {len(prevalence)}")
log_to_file(f"Core Species (>50% prevalence): {len(core_species)}")
log_to_file(f"List of Core Species: {', '.join(core_species)}")

# Venn Diagram-like logic (Overlap between Sexes)
core_F = (abundance_table[metadata[metadata['Sex']=='F'].index] > 0).sum(axis=1) / metadata[metadata['Sex']=='F'].shape[0]
core_M = (abundance_table[metadata[metadata['Sex']=='M'].index] > 0).sum(axis=1) / metadata[metadata['Sex']=='M'].shape[0]

# Strict core (>90% prevalence) for comparison
strict_F = set(core_F[core_F > 0.9].index)
strict_M = set(core_M[core_M > 0.9].index)
shared = strict_F.intersection(strict_M)
unique_F = strict_F - strict_M
unique_M = strict_M - strict_F

log_to_file(f"\nStrict Core (>90%) Comparison:")
log_to_file(f"Shared by both sexes: {len(shared)}")
log_to_file(f"Unique to Females: {len(unique_F)} {list(unique_F)}")
log_to_file(f"Unique to Males:   {len(unique_M)} {list(unique_M)}")

log_to_file("\nAnalysis Complete. Report saved.")
print(f"\nAll outputs saved to: {output_dir}")