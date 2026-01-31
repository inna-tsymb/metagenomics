import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --- Step 0: Configuration & Setup ---
# Define your Input and Output paths clearly here
input_file = '/Users/user/Documents/metagenomics/KSE_microbiome/lecture_01/input/abundance_table.csv'
output_dir = '/Users/user/Documents/metagenomics/KSE_microbiome/lecture_01/output/'

# Ensure output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Created output directory: {output_dir}")

# Set visual style
sns.set_style("whitegrid")
plt.rcParams.update({'font.size': 12, 'figure.figsize': (10, 6)})

# --- Step 1: Data Loading ---
print(f"Loading data from: {input_file}")
try:
    # index_col=0 assumes the first column contains Species names
    abundance_table = pd.read_csv(input_file, index_col=0)
    print("Data loaded successfully.")
except FileNotFoundError:
    print(f"ERROR: File not found at {input_file}. Please check the path.")
    exit()

# --- Step 1.1: Safety Check ---
# It is crucial to check if the data loaded correctly before processing
print("\n--- Data Overview ---")
print(f"Dimensions: {abundance_table.shape[0]} Species x {abundance_table.shape[1]} Samples")
print("First 5 rows:")
print(abundance_table.head())

# --- Step 2: Data Preprocessing ---
# Fill missing values with 0
abundance_table.fillna(0, inplace=True)

# Force numeric conversion. 
# errors='coerce' turns non-numeric text (like 'Bacteria;Firmicutes...') into NaN, which we then fill with 0
df_counts = abundance_table.apply(pd.to_numeric, errors='coerce').fillna(0)
num_samples = df_counts.shape[1]

# --- Step 3: Statistical Summaries ---

print("\nCalculating statistics...")

# --- A. SPECIES Statistics ---
# 1. Relative Abundance
df_rel_abund = df_counts.div(df_counts.sum(axis=0), axis=1)

# 2. Prevalence & Mean Abundance
species_stats = pd.DataFrame({
    'Prevalence': (df_counts > 0).sum(axis=1) / num_samples,
    'Mean_Rel_Abund': df_rel_abund.mean(axis=1)
})

# --- B. SAMPLE Statistics (Alpha Diversity) ---
# 1. Richness & Total Reads
sample_richness = (df_counts > 0).sum(axis=0)
sample_total_reads = df_counts.sum(axis=0)

# 2. Shannon Diversity
# Mask zeros to avoid log(0) errors
rel_nonzero = df_rel_abund.where(df_rel_abund > 0) 
shannon_diversity = - (rel_nonzero * np.log(rel_nonzero)).sum(axis=0).fillna(0)

# 3. Pielou's Evenness
with np.errstate(divide='ignore', invalid='ignore'):
    pielou_evenness = shannon_diversity / np.log(sample_richness)
# Fix undefined cases where richness <= 1
pielou_evenness[sample_richness <= 1] = np.nan

sample_stats = pd.DataFrame({
    'Richness': sample_richness,
    'Total_Reads': sample_total_reads,
    'Shannon': shannon_diversity,
    'Evenness': pielou_evenness
})

# --- Step 4: Visualizations ---

def save_plot(filename):
    """Helper to save and close plots cleanly"""
    path = os.path.join(output_dir, filename)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()
    print(f"Saved plot: {filename}")

def plot_distribution(data, column, title, color, filename, x_label):
    """Helper for histograms"""
    plt.figure()
    sns.histplot(data[column], bins=20, kde=True, color=color, edgecolor='black', alpha=0.7)
    mean_val = data[column].mean()
    plt.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.2f}')
    plt.title(title, fontsize=15)
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.legend()
    save_plot(filename)

print("\nGenerating visualizations...")

# ==========================================
# PART 1: SPECIES PLOTS
# ==========================================

# 1.1 Distribution Histograms
plot_distribution(species_stats, 'Prevalence', 'Distribution of Species Prevalence', 
                  'teal', 'dist_species_prevalence.png', 'Prevalence')
plot_distribution(species_stats, 'Mean_Rel_Abund', 'Distribution of Mean Relative Abundance', 
                  'steelblue', 'dist_species_mean_abundance.png', 'Mean Rel. Abundance')

# 1.2 Bar Plot: Top Prevalent Species
top_n = 10
top_prevalent = species_stats.sort_values(by='Prevalence', ascending=False).head(top_n)
plt.figure()
sns.barplot(x=top_prevalent.index, y=top_prevalent['Prevalence'], palette='viridis')
plt.title(f'Top {top_n} Most Prevalent Bacterial Species', fontsize=16, fontweight='bold')
plt.ylabel('Prevalence', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.ylim(0, 1.05)
plt.grid(axis='y', linestyle='--', alpha=0.6)
save_plot(f'bar_species_top_{top_n}_prevalent.png')

# 1.3 Scatter Plot: Prevalence vs. Abundance (Core Microbiome Analysis)
plt.figure()
sns.scatterplot(
    data=species_stats, x='Prevalence', y='Mean_Rel_Abund',
    hue='Mean_Rel_Abund', palette='viridis', size='Mean_Rel_Abund', sizes=(50, 300),
    legend=False, alpha=0.8, edgecolor='black'
)
# Label top species
if not species_stats.empty:
    top_sp = species_stats['Mean_Rel_Abund'].idxmax()
    plt.text(species_stats.loc[top_sp, 'Prevalence']+0.02, species_stats.loc[top_sp, 'Mean_Rel_Abund'], 
             top_sp, fontsize=11, fontweight='bold')
plt.title('Species: Prevalence vs. Mean Relative Abundance', fontsize=16, fontweight='bold')
plt.xlabel('Prevalence', fontsize=14)
plt.ylabel('Mean Relative Abundance', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.5)
save_plot('scatter_species_prevalence_vs_abundance.png')

# ==========================================
# PART 2: SAMPLE PLOTS
# ==========================================

# 2.1 Distribution Histograms
plot_distribution(sample_stats, 'Richness', 'Distribution of Sample Richness', 
                  'forestgreen', 'dist_sample_richness.png', 'Richness')
plot_distribution(sample_stats, 'Shannon', 'Distribution of Sample Shannon Diversity', 
                  'rebeccapurple', 'dist_sample_shannon.png', 'Shannon Index')
plot_distribution(sample_stats, 'Evenness', "Distribution of Sample Pielou's Evenness", 
                  'darkorange', 'dist_sample_pielou.png', "Pielou's Evenness")

# 2.2 Bar Plot: Top Samples by Diversity
top_n_samples = 10
top_diversity = sample_stats.sort_values(by='Shannon', ascending=False).head(top_n_samples)
plt.figure()
sns.barplot(x=top_diversity.index, y=top_diversity['Shannon'], palette='magma')
plt.title(f'Top {top_n_samples} Samples by Shannon Diversity', fontsize=16, fontweight='bold')
plt.ylabel('Shannon Diversity Index', fontsize=14)
plt.xlabel('Sample ID', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.grid(axis='y', linestyle='--', alpha=0.6)
save_plot(f'bar_sample_top_{top_n_samples}_diversity.png')

# 2.3 Scatter Plot: Richness vs. Sequencing Depth (QC Check)
plt.figure()
sns.scatterplot(
    data=sample_stats, x='Richness', y='Total_Reads',
    hue='Total_Reads', palette='magma', size='Total_Reads', sizes=(50, 300),
    legend=False, alpha=0.8, edgecolor='black'
)
# Label top sample
if not sample_stats.empty:
    top_sa = sample_stats['Total_Reads'].idxmax()
    plt.text(sample_stats.loc[top_sa, 'Richness']+0.1, sample_stats.loc[top_sa, 'Total_Reads'], 
             top_sa, fontsize=11, fontweight='bold')
plt.title('Samples: Richness vs. Total Read Counts', fontsize=16, fontweight='bold')
plt.xlabel('Richness (Number of Species)', fontsize=14)
plt.ylabel('Total Read Counts (Sequencing Depth)', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.5)
save_plot('scatter_sample_richness_vs_reads.png')

# --- Step 5: Save Statistics to CSV ---
species_stats.to_csv(os.path.join(output_dir, 'species_statistics.csv'))
sample_stats.to_csv(os.path.join(output_dir, 'sample_statistics.csv'))
print("\nAll processing complete. Statistics and plots saved to:", output_dir)


from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ==========================================
# PART 3: ADVANCED PLOTS
# ==========================================
print("\nGenerating advanced plots (PCA, Heatmap, Rank-Abundance)...")

# --- 3.1 PCA (Principal Component Analysis) ---
# We use Relative Abundance. 
# Log-transforming helps normalize the highly skewed microbiome data.
# We add a small constant (1e-6) to handle zeros before logging.
log_data = np.log(df_rel_abund + 1e-6)

# Transpose so rows = samples (sklearn expects samples as rows)
X = log_data.T 

# Run PCA
pca = PCA(n_components=2)
principal_components = pca.fit_transform(X)
pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=X.index)

# Calculate explained variance (How much info is in this plot?)
var_exp = pca.explained_variance_ratio_ * 100

plt.figure(figsize=(8, 6))
sns.scatterplot(x='PC1', y='PC2', data=pca_df, s=100, color='darkblue', alpha=0.7, edgecolor='black')

# Label points (optional, good if you have few samples)
if num_samples < 50:
    for sample_name, row in pca_df.iterrows():
        plt.text(row['PC1']+0.2, row['PC2']+0.2, sample_name, fontsize=9)

plt.title(f'PCA of Microbiome Data\n(PC1: {var_exp[0]:.1f}%, PC2: {var_exp[1]:.1f}%)', fontsize=14)
plt.xlabel(f'Principal Component 1 ({var_exp[0]:.1f}%)')
plt.ylabel(f'Principal Component 2 ({var_exp[1]:.1f}%)')
plt.grid(True, linestyle='--', alpha=0.5)
save_plot('pca_plot.png')


# --- 3.2 Clustered Heatmap (Top 20 Species) ---
# Plotting 1000+ species is messy. Let's filter for the Top 20 most abundant.
top_20_species = species_stats.sort_values('Mean_Rel_Abund', ascending=False).head(20).index
heatmap_data = df_rel_abund.loc[top_20_species]

# Use clustermap to group similar samples and species
# z_score=0 normalizes rows (shows if abundance is higher/lower than THAT species' average)
g = sns.clustermap(heatmap_data, cmap='viridis', standard_scale=None, z_score=0, 
                   figsize=(10, 10), method='ward', metric='euclidean')

g.fig.suptitle('Clustered Heatmap (Top 20 Species)', fontsize=16, y=1.02)
# Save directly using the clustergrid object
g.savefig(os.path.join(output_dir, 'clustered_heatmap.png'))
plt.close()
print("Saved plot: clustered_heatmap.png")


# --- 3.3 Rank-Abundance Curve ---
plt.figure(figsize=(10, 6))

# Sort species by abundance for every sample (aggregate mean)
sorted_abundance = species_stats['Mean_Rel_Abund'].sort_values(ascending=False).values

# Create the rank (0, 1, 2, ...)
ranks = np.arange(1, len(sorted_abundance) + 1)

plt.plot(ranks, sorted_abundance, marker='o', linestyle='-', color='black', markersize=4)

plt.yscale('log') # Log scale is standard for this plot
plt.title('Rank-Abundance Curve (Whittaker Plot)', fontsize=15)
plt.xlabel('Species Rank (Most to Least Abundant)', fontsize=12)
plt.ylabel('Mean Relative Abundance (Log Scale)', fontsize=12)
plt.grid(True, which="both", ls="--", alpha=0.4)

# Highlight the "Long Tail"
plt.text(len(ranks)*0.6, sorted_abundance[0], 'Dominant Species', fontsize=10, color='red')
plt.text(len(ranks)*0.6, sorted_abundance[-1], 'Rare Biosphere\n(Long Tail)', fontsize=10, color='blue', va='bottom')

save_plot('rank_abundance_curve.png')

print("\nAdvanced analysis complete.")