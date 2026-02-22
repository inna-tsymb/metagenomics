# Species-Group Association Statements

## Update (2026-02-22): Dual-Mode Availability

- This report is aligned with two execution modes: **normalized** and **non_normalized**.
- Normalized mode applies lecture_02-style cohort balancing: **Schulz_2017_wastewater -> n=20**.
- Non-normalized mode preserves original cohort sizes for impact comparison.
- Full artifacts (CSV, PNG, TXT/MD reports) are available in mode-specific output folders.
- Re-run all 03/04 dual-mode analyses with: `./run_all_dual_mode.sh`

## Summary Statistics

- **Total species analyzed:** 859
- **Significant associations (FDR < 0.05):** 758
- **Statements generated for:** Top 30 species

## Interpretation Guide

- **Positive Association (✓):** Species is enriched (higher mean abundance) in this group
- **Negative Association (✗):** Species is depleted (lower mean abundance) in this group
- **β (Beta):** Effect size in CLR units; deviation from global mean abundance
  - β > 0 = enriched relative to average
  - β < 0 = depleted relative to average
- **Effect size range:** 0.64 to 15.43 CLR units
- **Mean R²:** 0.237 (average variance explained by group membership)

---

## Top Associated Species


### 1. **Limnohabitans sp Rim47**

**Statistical Significance:** p = 7.30e-60 | FDR = 6.27e-57 | Effect size = 15.425 CLR units | R² = 0.942

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +12.9511 |
| **Mean abundance (enriched)** | 13.2259 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.4740 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=7.30e-60). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +12.9511 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.4740 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 2. **Candidatus Planktophila sulfonica**

**Statistical Significance:** p = 1.36e-54 | FDR = 4.86e-52 | Effect size = 13.816 CLR units | R² = 0.926

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +11.5333 |
| **Mean abundance (enriched)** | 11.6169 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.2829 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.36e-54). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +11.5333 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.2829 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 3. **Limnohabitans sp 63ED37 2**

**Statistical Significance:** p = 1.70e-54 | FDR = 4.86e-52 | Effect size = 13.456 CLR units | R² = 0.925

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +11.2162 |
| **Mean abundance (enriched)** | 11.2571 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.2401 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.70e-54). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +11.2162 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.2401 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 4. **Limnohabitans sp G3 2**

**Statistical Significance:** p = 2.58e-53 | FDR = 5.54e-51 | Effect size = 13.025 CLR units | R² = 0.921

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.8360 |
| **Mean abundance (enriched)** | 10.8256 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.1888 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=2.58e-53). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.8360 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.1888 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 5. **Polynucleobacter sp es MAR 4**

**Statistical Significance:** p = 4.78e-53 | FDR = 8.21e-51 | Effect size = 12.160 CLR units | R² = 0.920

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.0742 |
| **Mean abundance (enriched)** | 9.9611 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.0861 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=4.78e-53). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.0742 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.0861 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 6. **GGB46525 SGB64386**

**Statistical Significance:** p = 1.39e-52 | FDR = 1.99e-50 | Effect size = 12.612 CLR units | R² = 0.918

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.4720 |
| **Mean abundance (enriched)** | 10.4126 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.1398 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.39e-52). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.4720 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.1398 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 7. **Limnohabitans sp 2KL 3**

**Statistical Significance:** p = 2.82e-52 | FDR = 3.46e-50 | Effect size = 12.669 CLR units | R² = 0.917

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.5220 |
| **Mean abundance (enriched)** | 10.4693 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.1465 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=2.82e-52). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.5220 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.1465 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 8. **Limnohabitans sp Rim11**

**Statistical Significance:** p = 9.87e-52 | FDR = 1.06e-49 | Effect size = 12.666 CLR units | R² = 0.915

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.5202 |
| **Mean abundance (enriched)** | 10.4672 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.1463 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=9.87e-52). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.5202 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.1463 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 9. **Candidatus Methylopumilus universalis**

**Statistical Significance:** p = 2.87e-51 | FDR = 2.59e-49 | Effect size = 12.468 CLR units | R² = 0.913

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.3457 |
| **Mean abundance (enriched)** | 10.2692 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.1227 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=2.87e-51). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.3457 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.1227 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 10. **Polynucleobacter sp MWH CaK5**

**Statistical Significance:** p = 3.01e-51 | FDR = 2.59e-49 | Effect size = 12.476 CLR units | R² = 0.913

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.3523 |
| **Mean abundance (enriched)** | 10.2767 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.1236 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=3.01e-51). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.3523 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.1236 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 11. **Polynucleobacter sp AP Ainpum 60 G11**

**Statistical Significance:** p = 4.35e-51 | FDR = 3.40e-49 | Effect size = 12.520 CLR units | R² = 0.912

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.3912 |
| **Mean abundance (enriched)** | 10.3208 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.1289 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=4.35e-51). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.3912 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.1289 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 12. **Polynucleobacter sp AM 26B4**

**Statistical Significance:** p = 8.89e-51 | FDR = 6.36e-49 | Effect size = 12.338 CLR units | R² = 0.911

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.2309 |
| **Mean abundance (enriched)** | 10.1389 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.1073 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=8.89e-51). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.2309 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.1073 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 13. **Polynucleobacter cosmopolitanus**

**Statistical Significance:** p = 3.43e-50 | FDR = 2.27e-48 | Effect size = 12.078 CLR units | R² = 0.908

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.0016 |
| **Mean abundance (enriched)** | 9.8787 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.0763 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=3.43e-50). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.0016 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.0763 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 14. **GGB24856 SGB81948**

**Statistical Significance:** p = 4.35e-50 | FDR = 2.67e-48 | Effect size = 12.278 CLR units | R² = 0.908

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.1776 |
| **Mean abundance (enriched)** | 10.0784 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.1001 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=4.35e-50). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.1776 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.1001 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 15. **Candidatus Fonsibacter ubiquis**

**Statistical Significance:** p = 9.59e-50 | FDR = 5.49e-48 | Effect size = 12.088 CLR units | R² = 0.906

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.0104 |
| **Mean abundance (enriched)** | 9.8887 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.0775 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=9.59e-50). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.0104 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.0775 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 16. **GGB46527 SGB64388**

**Statistical Significance:** p = 2.05e-49 | FDR = 1.10e-47 | Effect size = 11.884 CLR units | R² = 0.905

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.8310 |
| **Mean abundance (enriched)** | 9.6852 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.0534 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=2.05e-49). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.8310 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.0534 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 17. **Candidatus Planktophila vernalis**

**Statistical Significance:** p = 2.46e-49 | FDR = 1.25e-47 | Effect size = 12.112 CLR units | R² = 0.905

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +10.0313 |
| **Mean abundance (enriched)** | 9.9124 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.0804 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=2.46e-49). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +10.0313 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.0804 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 18. **Limnohabitans sp 2KL 27**

**Statistical Significance:** p = 4.76e-49 | FDR = 2.27e-47 | Effect size = 11.777 CLR units | R² = 0.903

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.7364 |
| **Mean abundance (enriched)** | 9.5778 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.0406 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=4.76e-49). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.7364 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.0406 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 19. **Bradyrhizobium SGB11337**

**Statistical Significance:** p = 3.90e-48 | FDR = 1.76e-46 | Effect size = 10.291 CLR units | R² = 0.899

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +8.4267 |
| **Mean abundance (enriched)** | 8.0915 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.8640 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=3.90e-48). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +8.4267 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.8640 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 20. **GGB34754 SGB82226**

**Statistical Significance:** p = 1.27e-47 | FDR = 5.45e-46 | Effect size = 11.882 CLR units | R² = 0.896

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.8294 |
| **Mean abundance (enriched)** | 9.6833 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.0531 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.27e-47). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.8294 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.0531 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 21. **Polynucleobacter sp MWH Jannik1A5**

**Statistical Significance:** p = 1.49e-47 | FDR = 6.09e-46 | Effect size = 11.708 CLR units | R² = 0.896

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.6757 |
| **Mean abundance (enriched)** | 9.5089 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.0324 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.49e-47). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.6757 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.0324 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 22. **GGB25723 SGB106723**

**Statistical Significance:** p = 8.70e-47 | FDR = 3.40e-45 | Effect size = 11.568 CLR units | R² = 0.892

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.5522 |
| **Mean abundance (enriched)** | 9.3688 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -2.0158 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=8.70e-47). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.5522 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -2.0158 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 23. **Polynucleobacter asymbioticus**

**Statistical Significance:** p = 3.34e-46 | FDR = 1.25e-44 | Effect size = 11.072 CLR units | R² = 0.889

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.1148 |
| **Mean abundance (enriched)** | 8.8723 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.9568 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=3.34e-46). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.1148 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.9568 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 24. **Limnohabitans sp Hippo3**

**Statistical Significance:** p = 6.42e-46 | FDR = 2.30e-44 | Effect size = 11.076 CLR units | R² = 0.888

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.1183 |
| **Mean abundance (enriched)** | 8.8764 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.9573 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=6.42e-46). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.1183 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.9573 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 25. **GGB26951 SGB39152**

**Statistical Significance:** p = 7.05e-45 | FDR = 2.42e-43 | Effect size = 10.738 CLR units | R² = 0.882

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +8.8207 |
| **Mean abundance (enriched)** | 8.5386 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.9171 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=7.05e-45). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +8.8207 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.9171 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 26. **Candidatus Planktophila limnetica**

**Statistical Significance:** p = 7.81e-45 | FDR = 2.58e-43 | Effect size = 11.039 CLR units | R² = 0.882

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.0864 |
| **Mean abundance (enriched)** | 8.8401 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.9530 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=7.81e-45). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.0864 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.9530 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 27. **GGB34383 SGB60872**

**Statistical Significance:** p = 1.00e-44 | FDR = 3.19e-43 | Effect size = 11.306 CLR units | R² = 0.881

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.3212 |
| **Mean abundance (enriched)** | 9.1066 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.9846 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.00e-44). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.3212 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.9846 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 28. **Aquirufa lenticrescens**

**Statistical Significance:** p = 1.15e-44 | FDR = 3.54e-43 | Effect size = 10.614 CLR units | R² = 0.881

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +8.7115 |
| **Mean abundance (enriched)** | 8.4147 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.9024 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.15e-44). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +8.7115 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.9024 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 29. **GGB32489 SGB48813**

**Statistical Significance:** p = 1.85e-44 | FDR = 5.47e-43 | Effect size = 11.002 CLR units | R² = 0.880

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.0537 |
| **Mean abundance (enriched)** | 8.8030 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.9485 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.85e-44). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.0537 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.9485 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.


### 30. **GGB25723 SGB84803**

**Statistical Significance:** p = 1.96e-44 | FDR = 5.61e-43 | Effect size = 10.961 CLR units | R² = 0.880

| Metric | Value |
|--------|-------|
| **Enriched in** | Lekunberri_2018_river_wastewater |
| **β (enriched)** | +9.0171 |
| **Mean abundance (enriched)** | 8.7615 |
| **Depleted in** | Rowe_2017_hospital_wastewater |
| **β (depleted)** | -1.9436 |
| **Mean abundance (depleted)** | -2.1992 |

**Interpretation:** This species is strongly associated with wastewater source (p=1.96e-44). 
It is substantially enriched in **Lekunberri_2018_river_wastewater** 
(effect size: +9.0171 CLR units) and depleted in **Rowe_2017_hospital_wastewater** 
(effect size: -1.9436 CLR units). This species may serve as a biomarker 
for Lekunberri_2018_river_wastewater wastewater.

