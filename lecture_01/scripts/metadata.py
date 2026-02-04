import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import datetime

# ==========================================
# 0. CONFIGURATION & SETUP
# ==========================================
input_dir = '/Users/user/Documents/metagenomics/lecture_01/input/'
base_output_dir = '/Users/user/Documents/metagenomics/lecture_01/output/'

# --- Create 'metadata' subfolder ---
output_dir = os.path.join(base_output_dir, 'metadata')
if not os.path.exists(output_dir):
    try:
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")
    except OSError as e:
        print(f"Error creating directory: {e}")

# Files
metadata_file = os.path.join(input_dir, 'metadata_table.csv')
abundance_file = os.path.join(input_dir, 'abundance_table.csv')

# Visual Setup
sns.set_style("whitegrid")
plt.rcParams.update({'font.size': 12, 'figure.figsize': (10, 6)})

# Logging function
def log_to_file(message, logfile=None):
    if logfile is None:
        logfile = os.path.join(output_dir, 'metadata_summary_report.txt')
    # Use 'a' to append, create if not exists
    with open(logfile, 'a') as f:
        f.write(message + "\n")
    print(message)

# Clear old log file if needed
log_path = os.path.join(output_dir, 'metadata_summary_report.txt')
if os.path.exists(log_path):
    os.remove(log_path)

log_to_file("=== METADATA & COHORT ANALYSIS REPORT ===\n")

# ==========================================
# 1. DATA LOADING
# ==========================================
# Load Metadata
if os.path.exists(metadata_file):
    metadata = pd.read_csv(metadata_file, index_col=0)
else:
    print(f"Error: {metadata_file} not found.")
    exit()

# Load Abundance
if os.path.exists(abundance_file):
    abundance_table = pd.read_csv(abundance_file, index_col=0).fillna(0)
else:
    print(f"Error: {abundance_file} not found.")
    exit()

# Standardize Sex/Gender column
if 'Sex' not in metadata.columns and 'Gender' in metadata.columns:
    if pd.api.types.is_integer_dtype(metadata['Gender']) or pd.api.types.is_float_dtype(metadata['Gender']):
        metadata['Sex'] = metadata['Gender'].map({1: 'F', 0: 'M'})
    else:
        metadata['Sex'] = metadata['Gender'].astype(str)

# Align Data (Intersection of samples)
common_samples = abundance_table.columns.intersection(metadata.index)
# Check orientation if no overlap
if len(common_samples) == 0:
    abundance_table = abundance_table.T
    common_samples = abundance_table.columns.intersection(metadata.index)

if len(common_samples) == 0:
    log_to_file("ERROR: No overlapping sample IDs.")
    exit()

# Filter data
abundance_table = abundance_table[common_samples]
metadata = metadata.loc[common_samples]
# Relative Abundance
df_rel = abundance_table.div(abundance_table.sum(axis=0), axis=1)

log_to_file(f"Samples analyzed: {len(common_samples)}")

# ==========================================
# 2. DEMOGRAPHICS (Counts & Age Stats)
# ==========================================
log_to_file("\n--- 1. Demographics Summary ---")

# Sex Distribution
sex_counts = metadata['Sex'].value_counts()
log_to_file("Sex Distribution:")
log_to_file(sex_counts.to_string())

# Age Distribution Stats
age_stats = metadata['Age'].describe()
log_to_file("\nAge Statistics (Entire Cohort):")
log_to_file(f"Mean: {age_stats['mean']:.2f}, Median: {age_stats['50%']:.2f}, Std: {age_stats['std']:.2f}")
log_to_file(f"Min: {age_stats['min']}, Max: {age_stats['max']}")

# --- PLOT 1: Sex Pie Chart ---
plt.figure(figsize=(6, 6))
plt.pie(sex_counts, labels=sex_counts.index, autopct='%1.1f%%', colors=['#ff9999','#66b3ff'], startangle=90)
plt.title('Cohort Composition by Sex')
plt.savefig(os.path.join(output_dir, '1_sex_distribution_pie.png'))
plt.close()

# --- PLOT 2: Age Histogram ---
plt.figure(figsize=(8, 6))
sns.histplot(metadata['Age'], bins=15, kde=True, color='skyblue')
# Add vertical mean line and legend
mean_age = metadata['Age'].mean()
plt.axvline(mean_age, color='red', linestyle='--', linewidth=2, label=f"Mean: {mean_age:.1f}")
plt.title('Age Distribution of Cohort')
plt.xlabel('Age (Years)')
plt.legend()
plt.savefig(os.path.join(output_dir, '2_age_histogram.png'))
plt.close()

# ==========================================
# 3. AGE COMPARISON M vs F (Statistical Tests)
# ==========================================
log_to_file("\n--- 2. Age Comparison: Males vs Females ---")

males_age = metadata[metadata['Sex'] == 'M']['Age']
females_age = metadata[metadata['Sex'] == 'F']['Age']

log_to_file(f"Males (n={len(males_age)}): Mean Age = {males_age.mean():.2f}")
log_to_file(f"Females (n={len(females_age)}): Mean Age = {females_age.mean():.2f}")

# Statistical Test Selection for AGE
# 1. Check Normality (Shapiro-Wilk)
stat_m, p_m = stats.shapiro(males_age)
stat_f, p_f = stats.shapiro(females_age)

log_to_file(f"\nNormality Check for Age (Shapiro-Wilk):")
log_to_file(f"Males p={p_m:.4f}, Females p={p_f:.4f}")

if p_m > 0.05 and p_f > 0.05:
    log_to_file("  -> Data looks normal. Using Student's T-test.")
    stat_test, p_val = stats.ttest_ind(males_age, females_age)
    test_name = "T-test"
else:
    log_to_file("  -> Data deviates from normality. Using Mann-Whitney U Test.")
    stat_test, p_val = stats.mannwhitneyu(males_age, females_age)
    test_name = "Mann-Whitney U"

log_to_file(f"\nAge Difference Test ({test_name}): p-value = {p_val:.5f}")
if p_val < 0.05:
    log_to_file("RESULT: Significant difference in age between sexes.")
else:
    log_to_file("RESULT: No significant difference in age between sexes.")

# --- PLOT 3: Age Boxplot ---
plt.figure(figsize=(6, 6))
ax = sns.boxplot(x='Sex', y='Age', data=metadata, hue='Sex', palette={'F': '#ff9999', 'M': '#66b3ff'}, dodge=False)
# Add dots (jittered) on top of the boxplot to show individual samples
sns.stripplot(x='Sex', y='Age', data=metadata,
              color='0.2', # dark grey points
              size=6,
              jitter=0.25,
              alpha=0.7,
              edgecolor='white',
              linewidth=0.3,
              zorder=10)
# Remove legend since hue duplicates x
if ax.get_legend() is not None:
    ax.get_legend().remove()
plt.title(f'Age Comparison by Sex\n (Student\'s t-test (Parametric), p={p_val:.4f})')
plt.savefig(os.path.join(output_dir, '3_age_boxplot_sex.png'))
plt.close()

# --- PLOT 4: Population Pyramid ---
# Create bins
bins = np.arange(0, 100, 5) # 0, 5, 10 ...
m_hist, _ = np.histogram(males_age, bins=bins)
f_hist, _ = np.histogram(females_age, bins=bins)

plt.figure(figsize=(10, 6))
# Males on left (negative), Females on right (positive)
plt.barh(bins[:-1], -m_hist, height=4, label='Males', color='cornflowerblue', align='center')
plt.barh(bins[:-1], f_hist, height=4, label='Females', color='lightpink', align='center')
plt.axvline(0, color='black', linewidth=0.8)
plt.title('Population Pyramid (Age Structure)')
plt.xlabel('Count')
plt.ylabel('Age Group')
plt.legend()
# Fix x-axis labels to be positive on both sides (use set_xticks before set_xticklabels to avoid UserWarning)
ax = plt.gca()
xticks = ax.get_xticks()
ax.set_xticks(xticks)
ax.set_xticklabels([str(abs(int(x))) for x in xticks])
plt.savefig(os.path.join(output_dir, '4_population_pyramid.png'))
plt.close()

# ==========================================
# 4. MICROBIOME COMPARISON M vs F
# ==========================================
log_to_file("\n--- 3. Microbiome Comparison: M vs F ---")

# Split abundance tables
df_M = df_rel.loc[:, metadata['Sex'] == 'M']
df_F = df_rel.loc[:, metadata['Sex'] == 'F']

# Function to get top stats
def get_top_stats(df):
    # Prevalence: Fraction of samples where count > 0
    prev = (df > 0).sum(axis=1) / df.shape[1]
    # Abundance: Mean relative abundance
    abund = df.mean(axis=1)
    
    top_prev = prev.sort_values(ascending=False).head(5).index.tolist()
    top_abund = abund.sort_values(ascending=False).head(5).index.tolist()
    return top_prev, top_abund

top_prev_M, top_abund_M = get_top_stats(df_M)
top_prev_F, top_abund_F = get_top_stats(df_F)

log_to_file("\nMost Prevalent Species (Top 5):")
log_to_file(f"Males:   {top_prev_M}")
log_to_file(f"Females: {top_prev_F}")
intersection_prev = set(top_prev_M).intersection(set(top_prev_F))
log_to_file(f"Shared Prevalent: {len(intersection_prev)}/5")

log_to_file("\nMost Abundant Species (Top 5):")
log_to_file(f"Males:   {top_abund_M}")
log_to_file(f"Females: {top_abund_F}")
intersection_abund = set(top_abund_M).intersection(set(top_abund_F))
log_to_file(f"Shared Abundant: {len(intersection_abund)}/5")

# --- STRICT CORE ANALYSIS (>90%) ---
log_to_file("\n--- Strict Core Analysis (>90% Prevalence) ---")
# Calculate prevalence for each group
core_F = (df_F > 0).sum(axis=1) / df_F.shape[1]
core_M = (df_M > 0).sum(axis=1) / df_M.shape[1]

strict_F = set(core_F[core_F > 0.9].index)
strict_M = set(core_M[core_M > 0.9].index)
shared_core = strict_F.intersection(strict_M)
unique_F_core = strict_F - strict_M
unique_M_core = strict_M - strict_F

log_to_file(f"Strict Core Species (>90% prevalence):")
log_to_file(f"Shared by both: {len(shared_core)}")
log_to_file(f"Unique to Females: {len(unique_F_core)} {list(unique_F_core)}")
log_to_file(f"Unique to Males:   {len(unique_M_core)} {list(unique_M_core)}")

# ==========================================
# 5. DIVERSITY & CORRELATION ANALYSIS
# ==========================================
log_to_file("\n--- 4. Alpha Diversity Analysis ---")

# Calculate Shannon Diversity
shannon = - (df_rel * np.log(df_rel + 1e-9)).sum(axis=0)
metadata['Shannon'] = shannon

f_div = metadata[metadata['Sex']=='F']['Shannon']
m_div = metadata[metadata['Sex']=='M']['Shannon']

# Diversity Test (Mann-Whitney U is standard for indices)
stat, p_div = stats.mannwhitneyu(f_div, m_div)

log_to_file(f"Female Mean Shannon: {f_div.mean():.2f}")
log_to_file(f"Male Mean Shannon:   {m_div.mean():.2f}")
log_to_file(f"Diversity Test (Mann-Whitney): p={p_div:.5f}")

# --- PLOT 5: Diversity Boxplot ---
plt.figure(figsize=(6, 6))
ax = sns.boxplot(x='Sex', y='Shannon', data=metadata, hue='Sex', palette={'F': '#ff9999', 'M': '#66b3ff'}, dodge=False)
sns.stripplot(x='Sex', y='Shannon', data=metadata, color='black', alpha=0.3)
# Remove redundant legend
if ax.get_legend() is not None:
    ax.get_legend().remove()
plt.title(f'Shannon Diversity by Sex\n(p={p_div:.4f})')
plt.savefig(os.path.join(output_dir, '5_diversity_vs_sex.png'))
plt.close()

# --- ANALYSIS: Correlation Age vs Species ---
log_to_file("\n--- 5. Age vs Species Correlation ---")
top_10 = df_rel.mean(axis=1).sort_values(ascending=False).head(10).index
correlations = []

for species in top_10:
    corr, p_val = stats.spearmanr(df_rel.loc[species], metadata['Age'])
    correlations.append({'Species': species, 'Rho': corr, 'P_value': p_val})

corr_df = pd.DataFrame(correlations).sort_values('P_value')
log_to_file(corr_df.to_string(index=False))

# --- PLOT 6: Correlation Heatmap ---
plt.figure(figsize=(8, 6))
heatmap_data = corr_df.set_index('Species')[['Rho']]
sns.heatmap(heatmap_data, annot=True, cmap='coolwarm', center=0, vmin=-1, vmax=1)
plt.title('Correlation: Top 10 Species vs Age')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, '6_correlation_species_age.png'))
plt.close()

print(f"\nProcessing complete. All results saved to: {output_dir}")


print("\n--- Generating IMPROVED Visualizations ---")

# ==========================================
# 6. IMPROVED: SPECIES COMPARISON (Bar Chart)
# ==========================================
# 1. Identify the Union of Top Species (Male Top 5 + Female Top 5)
top_M_indices = df_rel.loc[:, metadata['Sex'] == 'M'].mean(axis=1).nlargest(5).index
top_F_indices = df_rel.loc[:, metadata['Sex'] == 'F'].mean(axis=1).nlargest(5).index
combined_top_species = list(set(top_M_indices) | set(top_F_indices))

# 2. Prepare Data for Plotting
# Filter original relative abundance table for these species only
df_viz = df_rel.loc[combined_top_species].T
df_viz['Sex'] = metadata['Sex'] # Add Sex column

# Melt (Unpivot) for Seaborn
df_melted = df_viz.melt(id_vars='Sex', var_name='Species', value_name='Relative Abundance')

# 3. Plot Grouped Bar Chart
plt.figure(figsize=(12, 8))
sns.barplot(
    data=df_melted, 
    y='Species', 
    x='Relative Abundance', 
    hue='Sex', 
    palette={'M': '#66b3ff', 'F': '#ff9999'}, # Blue for M, Pink for F
    orient='h',
    errorbar=('ci', 68), # Show Standard Error
    capsize=0.1
)

plt.title('Contrast in Top Bacterial Species: Males vs Females', fontsize=14, fontweight='bold')
plt.xlabel('Mean Relative Abundance')
plt.ylabel('')
plt.legend(title='Sex', loc='lower right')
plt.tight_layout()

# Save
save_path = os.path.join(output_dir, '7_improved_species_contrast.png')
plt.savefig(save_path, dpi=300)
print(f"Saved Improved Species Chart to: {save_path}")

# ==========================================
# 7. IMPROVED: RAINCLOUD STYLE PLOT (Diversity)
# ==========================================
plt.figure(figsize=(8, 6))

# 1. Violin (Density)
ax = sns.violinplot(x='Sex', y='Shannon', data=metadata, 
               hue='Sex', palette={'M': '#e6f2ff', 'F': '#ffe6e6'}, # Very light colors
               inner=None, linewidth=0, dodge=False)
# Remove legend created by hue (we already show Sex on x-axis)
if ax.get_legend() is not None:
    ax.get_legend().remove()

# 2. Strip (Raw Points) - Jittered (draw on same ax)
sns.stripplot(x='Sex', y='Shannon', data=metadata, 
              hue='Sex', palette={'M': '#0059b3', 'F': '#cc0000'}, # Darker colors
              size=4, alpha=0.6, jitter=0.2, dodge=False, ax=ax)
# Remove legend if added
if ax.get_legend() is not None:
    ax.get_legend().remove()

# 3. Boxplot (Summary) - Transparent, narrow
sns.boxplot(x='Sex', y='Shannon', data=metadata, 
            width=0.2, boxprops={'facecolor':'none', 'edgecolor':'black'},
            whiskerprops={'color':'black'}, medianprops={'color':'black'}, ax=ax)

plt.title('Detailed Diversity Distribution (Raincloud Style)')
save_path_rain = os.path.join(output_dir, '8_improved_diversity_raincloud.png')
plt.savefig(save_path_rain, dpi=300)
print(f"Saved Raincloud Plot to: {save_path_rain}")

# ==========================================
# 8. ADVANCED ANALYSIS (Enterotyping, Beta-Diversity, Rare Biosphere, F/B Ratio)
# ==========================================
print("\n--- Running Advanced Analysis ---")
log_to_file("\n=== ADVANCED ANALYSIS REPORT ===")

from sklearn.cluster import KMeans
from sklearn.manifold import MDS
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import MinMaxScaler

# ---------------------------------------------------------
# A. ENTEROTYPING (CLUSTERING)
# ---------------------------------------------------------
log_to_file("\n--- A. Enterotyping (K-Means Clustering) ---")

# We use relative abundance for clustering
X_cluster = df_rel.T # Samples as rows

# Try to find 2 natural clusters (Enterotypes)
kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
clusters = kmeans.fit_predict(X_cluster)
metadata['Enterotype'] = clusters

# Check if Enterotypes match Sex
crosstab = pd.crosstab(metadata['Sex'], metadata['Enterotype'])
log_to_file("Enterotype vs Sex Distribution:")
log_to_file(crosstab.to_string())

# ---------------------------------------------------------
# B. BETA-DIVERSITY (PCoA with Bray-Curtis)
# ---------------------------------------------------------
log_to_file("\n--- B. Beta-Diversity (PCoA / Bray-Curtis) ---")

# Calculate Bray-Curtis Dissimilarity Matrix
# (metric='braycurtis' requires scipy)
dist_matrix = pdist(df_rel.T, metric='braycurtis')
dist_square = squareform(dist_matrix)

# PCoA is essentially MDS on a distance matrix
# Use explicit parameters to avoid future sklearn warnings: use metric='precomputed', set n_init and init
mds = MDS(n_components=2, metric='precomputed', random_state=42, normalized_stress='auto', n_init=4, init='random')
pcoa_coords = mds.fit_transform(dist_square)

# Plot PCoA
plt.figure(figsize=(10, 8))
pcoa_df = pd.DataFrame(pcoa_coords, columns=['PCoA1', 'PCoA2'], index=metadata.index)
pcoa_df['Sex'] = metadata['Sex']
pcoa_df['Enterotype'] = metadata['Enterotype'].astype(str)

sns.scatterplot(data=pcoa_df, x='PCoA1', y='PCoA2', hue='Sex', style='Enterotype', 
                palette={'M': '#66b3ff', 'F': '#ff9999'}, s=80, alpha=0.8)
plt.title('Beta-Diversity: PCoA (Bray-Curtis)')
plt.xlabel('PCoA Axis 1')
plt.ylabel('PCoA Axis 2')
plt.legend(title='Sex / Enterotype', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, '9_beta_diversity_pcoa.png'))
plt.close()
print(f"Saved Beta-Diversity Plot to: {os.path.join(output_dir, '9_beta_diversity_pcoa.png')}")

# ---------------------------------------------------------
# C. RARE BIOSPHERE ANALYSIS
# ---------------------------------------------------------
log_to_file("\n--- C. Rare Biosphere Analysis ---")

# Definition: Species with max relative abundance < 1% across all samples
rare_threshold = 0.01 
max_abundance = df_rel.max(axis=1)
rare_species = max_abundance[max_abundance < rare_threshold].index.tolist()

log_to_file(f"Rare Species Count (<1% max abundance): {len(rare_species)}")
if len(rare_species) > 0:
    log_to_file(f"List: {rare_species}")
    
    # Check if rare species are unique to M or F
    rare_df = df_rel.loc[rare_species]
    rare_in_F = (rare_df.loc[:, metadata['Sex']=='F'] > 0).any(axis=1)
    rare_in_M = (rare_df.loc[:, metadata['Sex']=='M'] > 0).any(axis=1)
    
    unique_rare_F = rare_in_F[rare_in_F & ~rare_in_M].index.tolist()
    unique_rare_M = rare_in_M[rare_in_M & ~rare_in_F].index.tolist()
    
    log_to_file(f"Rare Unique to F: {unique_rare_F}")
    log_to_file(f"Rare Unique to M: {unique_rare_M}")
else:
    log_to_file("No species qualified as 'Rare' under this threshold.")

# ---------------------------------------------------------
# D. FIRMICUTES / BACTEROIDETES (F/B) RATIO
# ---------------------------------------------------------
log_to_file("\n--- D. F/B Ratio Analysis ---")

# Since we don't have a taxonomy file, we map based on common names manually.
# This is a heuristic approach for this specific dataset.
taxonomy_map = {
    'Firmicutes': [
        'Lactobacillus', 'Roseburia', 'Ruminococcus', 'Eubacterium', 
        'Clostridium', 'Enterococcus', 'Dialister', 'Streptococcus', 
        'Veillonella', 'Faecalibacterium', 'Blautia', 'Coprococcus'
    ],
    'Bacteroidetes': [ # Note: Phocaeicola is the new name for some Bacteroides
        'Bacteroides', 'Prevotella', 'Phocaeicola', 'Parabacteroides', 'Alistipes', 'Odoribacter'
    ],
    'Actinobacteria': ['Bifidobacterium', 'Collinsella'],
    'Proteobacteria': ['Escherichia', 'Klebsiella', 'Desulfovibrio']
}

def get_phylum(species_name):
    for phylum, keywords in taxonomy_map.items():
        for key in keywords:
            if key in species_name:
                return phylum
    return 'Other'

# Create Phylum Abundance Table
phylum_counts = pd.DataFrame(index=metadata.index)
for phylum in taxonomy_map.keys():
    # Sum relative abundance of all species belonging to this phylum
    species_in_phylum = [sp for sp in df_rel.index if get_phylum(sp) == phylum]
    if species_in_phylum:
        phylum_counts[phylum] = df_rel.loc[species_in_phylum].sum(axis=0)
    else:
        phylum_counts[phylum] = 0

# Calculate Ratio (add small epsilon to avoid division by zero)
phylum_counts['FB_Ratio'] = phylum_counts['Firmicutes'] / (phylum_counts['Bacteroidetes'] + 1e-6)
metadata['FB_Ratio'] = phylum_counts['FB_Ratio']

# Log stats
mean_fb = metadata.groupby('Sex')['FB_Ratio'].mean()
log_to_file(f"Mean F/B Ratio by Sex:\n{mean_fb.to_string()}")

# Test difference
fb_m = metadata[metadata['Sex']=='M']['FB_Ratio']
fb_f = metadata[metadata['Sex']=='F']['FB_Ratio']
stat, p_fb = stats.mannwhitneyu(fb_m, fb_f)
log_to_file(f"F/B Ratio Mann-Whitney p-value: {p_fb:.5f}")

# Plot F/B Ratio
plt.figure(figsize=(6, 6))
# Use log scale for Y axis because ratios can be skewed
ax = sns.boxplot(x='Sex', y='FB_Ratio', data=metadata, hue='Sex', palette={'F':'#ffcc99','M':'#99ccff'}, dodge=False)
# Remove redundant legend
if ax.get_legend() is not None:
    ax.get_legend().remove()
plt.yscale('log') 
plt.title(f'Firmicutes/Bacteroidetes Ratio (Log Scale)\n(p={p_fb:.4f})')
plt.ylabel('F/B Ratio (log)')
plt.savefig(os.path.join(output_dir, '10_fb_ratio_boxplot.png'))
plt.close()
print(f"Saved F/B Ratio Plot to: {os.path.join(output_dir, '10_fb_ratio_boxplot.png')}")

print("\nAdvanced analysis complete.")