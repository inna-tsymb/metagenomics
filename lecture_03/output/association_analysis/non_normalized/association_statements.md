# Species-Group Association Statements

## Update (2026-02-22): Dual-Mode Availability

- This report is aligned with two execution modes: **normalized** and **non_normalized**.
- Normalized mode applies lecture_02-style cohort balancing: **Schulz_2017_wastewater -> n=20**.
- Non-normalized mode preserves original cohort sizes for impact comparison.
- Full artifacts (CSV, PNG, TXT/MD reports) are available in mode-specific output folders.
- Re-run all 03/04 dual-mode analyses with: `./run_all_dual_mode.sh`

## Summary Statistics

- **Total species analyzed:** 567
- **Significant associations (FDR < 0.05):** 557
- **Statements generated for:** Top 30 species

## Interpretation Guide

- **Positive Association (✓):** Species is enriched (higher mean abundance) in this group
- **Negative Association (✗):** Species is depleted (lower mean abundance) in this group
- **β (Beta):** Effect size in CLR units; deviation from global mean abundance
  - β > 0 = enriched relative to average
  - β < 0 = depleted relative to average
- **Effect size range:** 0.95 to 15.25 CLR units
- **Mean R²:** 0.234 (average variance explained by group membership)

---

## Top Associated Species


### 1. **Limnohabitans sp Rim47**

**Statistical Significance:** p = 2.17e-102 | FDR = 1.23e-99 | Effect size = 15.247 CLR units | R² = 0.900

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +13.1976 |
| **Mean abundance (enriched)** | 12.9567 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.0497 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=2.17e-102). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +13.1976 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.0497 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 2. **Candidatus Planktophila sulfonica**

**Statistical Significance:** p = 3.45e-92 | FDR = 9.79e-90 | Effect size = 13.638 CLR units | R² = 0.875

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +11.6810 |
| **Mean abundance (enriched)** | 11.3478 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.9573 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=3.45e-92). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +11.6810 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.9573 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 3. **Limnohabitans sp 63ED37 2**

**Statistical Significance:** p = 2.63e-91 | FDR = 4.97e-89 | Effect size = 13.279 CLR units | R² = 0.872

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +11.3418 |
| **Mean abundance (enriched)** | 10.9880 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.9367 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=2.63e-91). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +11.3418 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.9367 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 4. **Limnohabitans sp G3 2**

**Statistical Significance:** p = 7.82e-89 | FDR = 1.11e-86 | Effect size = 12.847 CLR units | R² = 0.865

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.9351 |
| **Mean abundance (enriched)** | 10.5564 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.9119 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=7.82e-89). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.9351 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.9119 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 5. **GGB46525 SGB64386**

**Statistical Significance:** p = 5.57e-87 | FDR = 6.32e-85 | Effect size = 12.434 CLR units | R² = 0.859

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.5458 |
| **Mean abundance (enriched)** | 10.1434 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8882 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=5.57e-87). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.5458 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8882 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 6. **Limnohabitans sp 2KL 3**

**Statistical Significance:** p = 1.03e-86 | FDR = 9.74e-85 | Effect size = 12.491 CLR units | R² = 0.858

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.5993 |
| **Mean abundance (enriched)** | 10.2002 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8914 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.03e-86). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.5993 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8914 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 7. **Polynucleobacter sp es MAR 4**

**Statistical Significance:** p = 1.52e-86 | FDR = 1.23e-84 | Effect size = 11.982 CLR units | R² = 0.858

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.1202 |
| **Mean abundance (enriched)** | 9.6919 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8622 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.52e-86). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.1202 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8622 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 8. **Limnohabitans sp Rim11**

**Statistical Significance:** p = 5.54e-86 | FDR = 3.93e-84 | Effect size = 12.489 CLR units | R² = 0.856

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.5973 |
| **Mean abundance (enriched)** | 10.1980 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8913 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=5.54e-86). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.5973 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8913 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 9. **Candidatus Methylopumilus universalis**

**Statistical Significance:** p = 6.02e-85 | FDR = 3.56e-83 | Effect size = 12.291 CLR units | R² = 0.853

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.4106 |
| **Mean abundance (enriched)** | 10.0000 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8799 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=6.02e-85). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.4106 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8799 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 10. **Polynucleobacter sp MWH CaK5**

**Statistical Significance:** p = 6.27e-85 | FDR = 3.56e-83 | Effect size = 12.298 CLR units | R² = 0.853

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.4177 |
| **Mean abundance (enriched)** | 10.0075 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8804 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=6.27e-85). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.4177 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8804 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 11. **Polynucleobacter sp AP Ainpum 60 G11**

**Statistical Significance:** p = 8.10e-85 | FDR = 4.18e-83 | Effect size = 12.342 CLR units | R² = 0.852

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.4593 |
| **Mean abundance (enriched)** | 10.0516 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8829 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=8.10e-85). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.4593 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8829 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 12. **Polynucleobacter sp AM 26B4**

**Statistical Significance:** p = 5.20e-84 | FDR = 2.46e-82 | Effect size = 12.160 CLR units | R² = 0.850

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.2878 |
| **Mean abundance (enriched)** | 9.8698 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8725 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=5.20e-84). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.2878 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8725 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 13. **GGB24856 SGB81948**

**Statistical Significance:** p = 5.62e-83 | FDR = 2.45e-81 | Effect size = 12.100 CLR units | R² = 0.846

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.2308 |
| **Mean abundance (enriched)** | 9.8092 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8690 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=5.62e-83). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.2308 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8690 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 14. **Polynucleobacter cosmopolitanus**

**Statistical Significance:** p = 1.14e-82 | FDR = 4.63e-81 | Effect size = 11.900 CLR units | R² = 0.845

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.0425 |
| **Mean abundance (enriched)** | 9.6095 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8575 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.14e-82). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.0425 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8575 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 15. **Candidatus Fonsibacter ubiquis**

**Statistical Significance:** p = 4.07e-82 | FDR = 1.54e-80 | Effect size = 11.910 CLR units | R² = 0.843

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.0520 |
| **Mean abundance (enriched)** | 9.6196 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8581 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=4.07e-82). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.0520 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8581 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 16. **Candidatus Planktophila vernalis**

**Statistical Significance:** p = 1.26e-81 | FDR = 4.46e-80 | Effect size = 11.934 CLR units | R² = 0.841

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.0743 |
| **Mean abundance (enriched)** | 9.6432 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8595 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.26e-81). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.0743 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8595 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 17. **GGB46527 SGB64388**

**Statistical Significance:** p = 3.11e-81 | FDR = 1.04e-79 | Effect size = 11.707 CLR units | R² = 0.840

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.8601 |
| **Mean abundance (enriched)** | 9.4160 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8464 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=3.11e-81). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.8601 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8464 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 18. **Limnohabitans sp 2KL 27**

**Statistical Significance:** p = 1.69e-80 | FDR = 5.31e-79 | Effect size = 11.599 CLR units | R² = 0.837

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.7589 |
| **Mean abundance (enriched)** | 9.3086 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8402 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.69e-80). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.7589 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8402 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 19. **GGB34754 SGB82226**

**Statistical Significance:** p = 6.97e-79 | FDR = 2.08e-77 | Effect size = 11.705 CLR units | R² = 0.831

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.8583 |
| **Mean abundance (enriched)** | 9.4141 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8463 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=6.97e-79). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.8583 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8463 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 20. **Polynucleobacter sp MWH Jannik1A5**

**Statistical Significance:** p = 2.10e-78 | FDR = 5.95e-77 | Effect size = 11.530 CLR units | R² = 0.829

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.6940 |
| **Mean abundance (enriched)** | 9.2397 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8363 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=2.10e-78). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.6940 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8363 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 21. **Sphingorhabdus rigui**

**Statistical Significance:** p = 2.81e-78 | FDR = 7.58e-77 | Effect size = 13.642 CLR units | R² = 0.829

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +11.6349 |
| **Mean abundance (enriched)** | 11.3515 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.0071 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=2.81e-78). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +11.6349 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.0071 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 22. **GGB25723 SGB106723**

**Statistical Significance:** p = 4.29e-77 | FDR = 1.11e-75 | Effect size = 11.390 CLR units | R² = 0.824

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.5619 |
| **Mean abundance (enriched)** | 9.0996 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8282 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=4.29e-77). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.5619 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8282 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 23. **Polynucleobacter asymbioticus**

**Statistical Significance:** p = 3.35e-75 | FDR = 8.25e-74 | Effect size = 10.894 CLR units | R² = 0.817

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.0940 |
| **Mean abundance (enriched)** | 8.6032 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.7997 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=3.35e-75). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.0940 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.7997 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 24. **Limnohabitans sp Hippo3**

**Statistical Significance:** p = 7.65e-75 | FDR = 1.81e-73 | Effect size = 10.898 CLR units | R² = 0.815

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.0978 |
| **Mean abundance (enriched)** | 8.6072 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8000 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=7.65e-75). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.0978 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8000 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 25. **GGB34383 SGB60872**

**Statistical Significance:** p = 7.99e-74 | FDR = 1.81e-72 | Effect size = 11.128 CLR units | R² = 0.811

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.3148 |
| **Mean abundance (enriched)** | 8.8374 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8132 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=7.99e-74). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.3148 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8132 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 26. **Candidatus Planktophila limnetica**

**Statistical Significance:** p = 2.34e-73 | FDR = 5.10e-72 | Effect size = 10.862 CLR units | R² = 0.809

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.0636 |
| **Mean abundance (enriched)** | 8.5710 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.7979 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=2.34e-73). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.0636 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.7979 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 27. **GGB46797 SGB64687**

**Statistical Significance:** p = 5.55e-73 | FDR = 1.17e-71 | Effect size = 12.285 CLR units | R² = 0.807

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.3263 |
| **Mean abundance (enriched)** | 9.9949 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.9592 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=5.55e-73). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.3263 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.9592 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 28. **GGB32489 SGB48813**

**Statistical Significance:** p = 8.62e-73 | FDR = 1.75e-71 | Effect size = 10.824 CLR units | R² = 0.806

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.0286 |
| **Mean abundance (enriched)** | 8.5339 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.7958 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=8.62e-73). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.0286 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.7958 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 29. **GGB25723 SGB84803**

**Statistical Significance:** p = 1.14e-72 | FDR = 2.23e-71 | Effect size = 10.783 CLR units | R² = 0.806

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +8.9895 |
| **Mean abundance (enriched)** | 8.4923 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.7934 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.14e-72). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +8.9895 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.7934 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 30. **Candidatus Nanopelagicus SGB64093**

**Statistical Significance:** p = 1.73e-72 | FDR = 3.26e-71 | Effect size = 10.761 CLR units | R² = 0.805

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +8.9693 |
| **Mean abundance (enriched)** | 8.4709 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.7921 |
| **Mean abundance (depleted)** | -2.2905 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.73e-72). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +8.9693 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.7921 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.

