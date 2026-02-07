import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

# ==========================================
# 1. SETUP PATHS
# ==========================================
base_dir = '/Users/user/Documents/metagenomics/lecture_02'
input_dir = os.path.join(base_dir, 'input')
output_dir = os.path.join(base_dir, 'output', 'differential_abundance')
os.makedirs(output_dir, exist_ok=True)

# ==========================================
# 2. LOAD & PREPARE DATA (Robust Method)
# ==========================================
print("Loading data...")
try:
    meta_df = pd.read_csv(os.path.join(input_dir, 'environmental.tsv'), sep='\t')
    mp_path = os.path.join(input_dir, 'environmental_metaphlan4_2026-02-06.tsv')
    mp_df = pd.read_csv(mp_path, sep='\t')
    
    # Clean Column Names
    mp_df.columns = [c.replace('#', '').strip() for c in mp_df.columns]

    # Handle Long vs Wide Format
    if 'rel_abund' in mp_df.columns and 'sample_alias' in mp_df.columns:
        print("Detected LONG format. Pivoting...")
        if 'clade_name' in mp_df.columns: tax_col = 'clade_name'
        else: tax_col = mp_df.columns[1]
        
        # Filter Species (s__)
        mp_df = mp_df[mp_df[tax_col].str.contains('s__') & ~mp_df[tax_col].str.contains('t__')]
        abund_df = mp_df.pivot_table(index='sample_alias', columns=tax_col, values='rel_abund', aggfunc='mean').fillna(0)
    else:
        print("Detected WIDE format...")
        if 'clade_name' not in mp_df.columns: mp_df.rename(columns={mp_df.columns[0]: 'clade_name'}, inplace=True)
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

    print(f"Merged {len(full_df)} samples.")

except Exception as e:
    print(f"Error: {e}")
    exit()

# ==========================================
# 3. DEFINE GROUPS & BALANCE
# ==========================================
anthropogenic_studies = ['Rowe_2017_hospital_wastewater', 'Chu_2017_sludge']
municipal_study = 'Schulz_2017_wastewater'
environmental_studies = ['Chopyk_2020_pond', 'Lekunberri_2018_river_wastewater']

target_mask = full_df['study_code'].isin(anthropogenic_studies + environmental_studies + [municipal_study])
subset_df = full_df[target_mask].copy()

# Downsample Schulz
schulz_samples = subset_df[subset_df['study_code'] == municipal_study]
others = subset_df[subset_df['study_code'] != municipal_study]

if len(schulz_samples) > 20:
    schulz_subset = schulz_samples.sample(n=20, random_state=42)
else:
    schulz_subset = schulz_samples

final_df = pd.concat([others, schulz_subset])

def assign_group(study):
    if study in anthropogenic_studies or study == municipal_study: return 'Anthropogenic'
    elif study in environmental_studies: return 'Environmental'
    return 'Other'

final_df['Source_Type'] = final_df['study_code'].apply(assign_group)
print(f"Analysis Set: {len(final_df)} samples")

# ==========================================
# 4. DIFFERENTIAL ABUNDANCE CHECK
# ==========================================
print("Calculating Differential Abundance...")

# Identify Species Columns
metadata_cols = meta_df.columns.tolist() + ['Source_Type']
species_cols = [c for c in abund_df.columns if c in final_df.columns]

# Calculate Mean Abundance per Group
grouped_means = final_df.groupby('Source_Type')[species_cols].mean()

# Calculate Difference (Anthro - Env)
# Positive = Higher in Anthro
# Negative = Higher in Env
diff = grouped_means.loc['Anthropogenic'] - grouped_means.loc['Environmental']

# Extract Top 5 for each
top_anthro = diff.nlargest(5)
top_env = diff.nsmallest(5) # These are negative values

# Prepare Report Data
report_data = []

# Process Anthro Markers
for species, val in top_anthro.items():
    mean_anthro = grouped_means.loc['Anthropogenic', species]
    mean_env = grouped_means.loc['Environmental', species]
    report_data.append({
        'Species': species.split('|')[-1], # Short name
        'Group': 'Anthropogenic',
        'Mean_Abund_Anthro': mean_anthro,
        'Mean_Abund_Env': mean_env,
        'Difference': val
    })

# Process Environmental Markers
for species, val in top_env.items():
    mean_anthro = grouped_means.loc['Anthropogenic', species]
    mean_env = grouped_means.loc['Environmental', species]
    report_data.append({
        'Species': species.split('|')[-1],
        'Group': 'Environmental',
        'Mean_Abund_Anthro': mean_anthro,
        'Mean_Abund_Env': mean_env,
        'Difference': abs(val) # Magnitude of difference
    })

report_df = pd.DataFrame(report_data)

# Save Report
report_path = os.path.join(output_dir, 'top_species_report.txt')
with open(report_path, 'w') as f:
    f.write("TOP 5 MARKER SPECIES PER GROUP\n")
    f.write("==============================\n")
    f.write("(Based on difference in Mean Relative Abundance)\n\n")
    f.write(report_df.to_string(index=False))

print(f"Report saved to: {report_path}")

# ==========================================
# 5. VISUALIZATION
# ==========================================
# Prepare data for plotting (Long format of just the top 10 species)
top_10_species = top_anthro.index.tolist() + top_env.index.tolist()
plot_df = final_df.melt(id_vars=['Source_Type'], value_vars=top_10_species, 
                        var_name='Species', value_name='Abundance')

# Clean names for plot
plot_df['Species'] = plot_df['Species'].apply(lambda x: x.split('|')[-1])

plt.figure(figsize=(12, 8))
sns.barplot(data=plot_df, x='Abundance', y='Species', hue='Source_Type', palette={'Anthropogenic': '#d62728', 'Environmental': '#2ca02c'})
plt.title('Top 5 Indicator Species for Each Group')
plt.xlabel('Relative Abundance (%)')
plt.ylabel('')
plt.legend(title='Group')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'marker_species_barplot.png'))

print(f"Plot saved to: {output_dir}")