import pandas as pd

df = pd.read_csv('/Users/user/Documents/metagenomics/lecture_02/input/environmental_extended_wide.tsv', sep='\t', usecols=['study_code'])

# Search for specific wastewater studies
target_studies = ['Schulz_2017_wastewater', 'Chu_2017_sludge', 'Rowe_2017_hospital_wastewater', 'Lekunberri_2018_river_wastewater', 'Chopyk_2020_pond']

print('Searching for target studies...')
for study in target_studies:
    mask = df['study_code'].str.contains(study, case=False, na=False)
    count = mask.sum()
    print(f'{study}: {count} samples')

# Also search for patterns in study names
print('\nSearching for wastewater/sludge/sewage patterns...')
patterns = ['wastewater', 'sludge', 'sewage']
for pattern in patterns:
    mask = df['study_code'].str.contains(pattern, case=False, na=False)
    count = mask.sum()
    if count > 0:
        print(f'\n{pattern}: {count} samples')
        studies = df[mask]['study_code'].unique()
        for s in studies[:10]:  # Show first 10
            s_count = (df['study_code'] == s).sum()
            print(f'  - {s}: {s_count}')
