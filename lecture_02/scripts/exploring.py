import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# ==========================================
# 1. SETUP & DATA LOADING
# ==========================================
# Define paths (Adjust these if running elsewhere)
base_dir = '/Users/user/Documents/metagenomics/lecture_02'
input_dir = os.path.join(base_dir, 'input')
output_dir = os.path.join(base_dir, 'output', 'exploratory')
os.makedirs(output_dir, exist_ok=True)

print("--- Starting Comprehensive Metadata Analysis ---")

try:
    # Load files
    df_mapping = pd.read_csv(os.path.join(input_dir, 'environmental.tsv'), sep='\t')
    df_meta = pd.read_csv(os.path.join(input_dir, 'environmental_extended_wide.tsv'), sep='\t')

    # Merge Full Dataset (Inner Join)
    df = pd.merge(df_mapping, df_meta, 
                  left_on='ena_ers_sample_id', 
                  right_on='spire_sample_name', 
                  how='inner', suffixes=('', '_ext'))
    
    print(f"Data Loaded. Total Samples: {len(df)}")

except FileNotFoundError as e:
    print(f"Error: {e}")
    exit()

# Open a report file to write answers
report_path = os.path.join(output_dir, 'Full_Analysis_Report.txt')
with open(report_path, 'w') as f:
    
    # ==========================================
    # OBJECTIVE 1: Structure & Quality
    # ==========================================
    f.write("OBJECTIVE 1: DATASET STRUCTURE & QUALITY\n")
    f.write("========================================\n\n")
    
    # Q: Total Samples & Independent Studies
    total_samples = len(df)
    total_studies = df['study_code'].nunique()
    f.write(f"1. Total Samples: {total_samples}\n")
    f.write(f"2. Independent Studies: {total_studies}\n")
    f.write(f"   Studies: {', '.join(df['study_code'].unique())}\n\n")
    
    # Q: Min, Max, Mean samples per study
    study_counts = df['study_code'].value_counts()
    f.write("3. Samples per Study Statistics:\n")
    f.write(f"   - Minimum: {study_counts.min()}\n")
    f.write(f"   - Maximum: {study_counts.max()}\n")
    f.write(f"   - Mean:    {study_counts.mean():.2f}\n\n")
    
    # Q: Available Metadata Fields
    f.write(f"4. Available Metadata Fields ({len(df.columns)} total):\n")
    f.write(f"   {', '.join(df.columns[:10])} ... [truncated]\n\n")
    
    # Q: Missing Values Proportion
    f.write("5. Missing Data Proportion (Top 10 non-empty fields):\n")
    missing_prop = df.isnull().mean()
    # specific columns of interest
    cols_interest = ['collection_date', 'latitude', 'longitude', 'pH', 'temperature', 'depth_meters', 'instrument_model']
    for col in cols_interest:
        if col in df.columns:
            f.write(f"   - {col}: {missing_prop[col]*100:.1f}% missing\n")
        else:
            f.write(f"   - {col}: Not found in dataset\n")
    f.write("\n")

    # Q: Duplicates
    duplicates = df.duplicated(subset=['ena_ers_sample_id']).sum()
    f.write(f"9. Duplicate Samples (by ID): {duplicates}\n\n")

    # ==========================================
    # OBJECTIVE 2: Heterogeneity & Longitudinal
    # ==========================================
    f.write("OBJECTIVE 2: HETEROGENEITY & LONGITUDINAL SAMPLING\n")
    f.write("==================================================\n\n")

    # Q: Longitudinal Sampling (Repeated measures)
    # Logic: Look for same location/site + different dates
    f.write("7. Longitudinal Sampling Analysis:\n")
    if 'collection_date' in df.columns and 'location_name' in df.columns:
        # Group by study and location, count unique dates
        longitudinal = df.groupby(['study_code', 'location_name'])['collection_date'].nunique()
        repeated_sites = longitudinal[longitudinal > 1]
        
        if len(repeated_sites) > 0:
            f.write(f"   YES. Found {len(repeated_sites)} sites with repeated sampling over time.\n")
            f.write(f"   Examples:\n{repeated_sites.head().to_string()}\n")
        else:
            f.write("   NO. No sites found with multiple sampling timepoints.\n")
    else:
        f.write("   Cannot determine (Missing date or location metadata).\n")
    f.write("\n")

    # Q: Technical Heterogeneity
    # Logic: Check for 'instrument_model', 'library_strategy'
    f.write("8. Technical Heterogeneity:\n")
    tech_cols = [c for c in df.columns if 'instrument' in c or 'platform' in c or 'library' in c]
    if tech_cols:
        f.write(f"   Technical columns found: {tech_cols}\n")
        for c in tech_cols:
            if df[c].notnull().sum() > 0:
                f.write(f"   - {c} breakdown:\n{df[c].value_counts().to_string()}\n")
            else:
                f.write(f"   - {c} is present but EMPTY.\n")
    else:
        f.write("   No technical metadata (sequencing platform/library) available.\n")
        f.write("   (Cannot assess technical heterogeneity).\n")

print(f"Text Report generated: {report_path}")

# ==========================================
# 3. VISUALIZATIONS
# ==========================================
sns.set_theme(style="whitegrid")

# Visual 1: Bar Plot (Samples per Study)
plt.figure(figsize=(10, 6))
sns.barplot(x=study_counts.index, y=study_counts.values, palette='viridis')
plt.title('Number of Samples per Study')
plt.xlabel('Study Code')
plt.ylabel('Count')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, '1_Study_Distribution.png'))
plt.close()

# Visual 2: Missing Data Heatmap
# Select "Core" columns to make the map readable
core_cols = ['collection_date', 'latitude', 'longitude', 'pH', 'depth_meters', 'environment_biome', 'location_name']
valid_core_cols = [c for c in core_cols if c in df.columns]
if valid_core_cols:
    plt.figure(figsize=(12, 6))
    sns.heatmap(df[valid_core_cols].isnull(), cbar=False, cmap='YlGnBu', yticklabels=False)
    plt.title('Missing Data Patterns (Blue=Present, Yellow=Missing)')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '2_Missing_Data_Heatmap.png'))
    plt.close()

# Visual 3: Correlation Matrix (Numerical Variables)
# Select only numeric columns
numeric_df = df.select_dtypes(include=[np.number])
# Drop columns that are IDs or useless (like all NaNs)
numeric_df = numeric_df.dropna(axis=1, how='all')
# Filter out ID-like columns if they were detected as numbers (often TaxIDs)
cols_to_keep = [c for c in numeric_df.columns if 'id' not in c.lower() and 'tax' not in c.lower()]
numeric_df = numeric_df[cols_to_keep]

if numeric_df.shape[1] > 1:
    plt.figure(figsize=(10, 8))
    corr = numeric_df.corr()
    sns.heatmap(corr, annot=True, fmt=".2f", cmap='coolwarm', center=0)
    plt.title('Correlation Matrix of Environmental Variables')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '3_Correlation_Matrix.png'))
    plt.close()
else:
    print("Skipping Correlation Matrix: Not enough numerical data.")

# Visual 4: Timeline / Longitudinal
if 'collection_date' in df.columns:
    df['collection_date'] = pd.to_datetime(df['collection_date'], errors='coerce')
    df_sorted = df.sort_values('collection_date')
    
    plt.figure(figsize=(12, 6))
    sns.scatterplot(data=df_sorted, x='collection_date', y='study_code', hue='study_code', s=100, legend=False)
    plt.title('Temporal Coverage of Studies')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '4_Longitudinal_Timeline.png'))
    plt.close()

print(f"All visualizations saved to: {output_dir}")