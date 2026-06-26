# Glia_3740732_Data
Repository of series of data sets used in the publications of the paper "Dorso-Ventral and Night-Day Regulation of Extracellular K+ Dynamics in Mouse Hippocampal Astrocytes" in Glia

# Electrophysiological and Structural Dynamics of Astrocytic Currents across Circadian Rhythms (ZT) and Brain Regions

##  General Information

* Dataset DOI: [10.5281/zenodo.20933398] (https://doi.org/10.5281/zenodo.20933397)
* Associated Publication: Manuscript 3740732, titled "Dorso-Ventral and Night-Day Regulation of Extracellular K+ Dynamics in Mouse Hippocampal Astrocytes"
* Principal Author: Nariman KIANI, Aix Marseille Univ, INSERM, INS, Inst Neurosci Syst, Marseille, France
* Date of Data Collection: Between September, 2022 and May, 2026
* Funding Sources: 
European Union’s Horizon 2020 Research and Innovation Program under Grant Agreement Number 956325, the Fondation pour la Recherche Médicale (FRM) FDT202404018111, the Croatian Science Foundation under project number IP-2022-10-8493, and Petroleum Technology Development Fund (PTDF/ED/OSS/PHD/KA/2015/22).

---

##  Dataset Overview & File Structure

This repository contains 5 master datasets investigating astrocytic properties, potassium ($K^+$) dynamics, and immunohistochemical markers across different Zeitgeber Times (ZT3, ZT8, ZT15) and hippocampal regions (Dorsal Hippocampus vs. Ventral Hippocampus).

### File Inventory
1.	`K_cntr_10Hz30s_norm` - Normalized K+ dynamic data, including LFP_integral, Peak amplitude, and Tau decay
2.	`Rise_time_Norm.CSV` - Normalized 25%-75% rise time assessment.
3.	`Control_undershoot.CSV` - Peak undershoot amplitude .
4.	`Coupling.CSV` - Astroglial dye-coupling cell count data.
5.	`Control_MFA_ephys.CSV` - Patch-clamp electrophysiology data involving MFA and Barium-chloride ($BaCl_2$) blockade ratios.
6.	`Microscopy.CSV` - Immunohistochemistry integrated density data for Kir4.1 and GFAP.
---

##  Data Dictionary & Variable Definitions

### 1. Factor/Design Variables (Common across datasets)
* `Region` / `Place`: Brain region analyzed. 
  * `DH`: Dorsal Hippocampus
  * `VH`: Ventral Hippocampus
* `ZT`: Zeitgeber Time point of tissue collection / recording (`3`, `8`, `15`).
* `Date`: Date of experiment execution (Formatted as YYYY-MM-DD via `lubridate`).
* `Animal`: Unique identifier code for the subject animal.

### 2. File-Specific Variables

### `K_cntr_10Hz30s_norm`
* `LFP_Integral`: Local field potential (LFP) integral, measure of synaptic strength
* `Peak`: Peak K+ amplitude in response to 10 Hz 30 second stimulation
* `Tau`: Monoexponentially fitted decay of the K+ transient after reaching peak
Values are baseline-normalized to the mean of `DH` at `ZT3`.
R_script: K_dynamics_10Hz30s_R

#### `Rise_time_Norm.CSV`
* `Slope`: Normalized 25%–75% rise time slope value. Values are baseline-normalized to the mean of `DH` at `ZT3`.
R_script: Rise_time_R

#### `Control_undershoot.CSV`
* `U_shoot`: Raw membrane potential undershoot value (mM).
* `U_norm`: Raw undershoot values normalized to the mean of `DH` at `ZT3`.
* `T_U`: Time (in seconds) required to reach the peak undershoot.
R_script: Undershoot_analysis_R

#### `Coupling.CSV`
* `Count`: Total number of biocytin/dye-coupled astrocytes surrounding the recorded cell.
R_script: Coupling_R

#### `Control_MFA_ephys.CSV`
* `Em1`: Initial resting membrane potential (RMP, measured in mV).
* `Diff`: Barium-dependent depolarization value Em2 – Em1 (in mV).
* `PercgBa`: Barium-sensitive conductance at -130mV holding potential representing Kir4.1 current magnitude (expressed as a decimal/percentage).
* `Blockade_ratio`: Absolute ratio of current inhibition following pharmacology.
* `slope_ratio`: Absolute value representing the rectification profile of astrocytic currents (Values > 1.0 indicate rectification).
* `Rectification`: Binary logical indicator (`TRUE` = Rectifying cell, `FALSE` = Non-rectifying cell).
Additional variables:
* `t_incub`: incubation time in Meclofenamic acid (MFA) in minutes
*  `Rp`:  pipette resistance in MΩ
* `Rs1`: Initial series resistance in MΩ
* `Rs2.Rs1`: Ratio of series resistance after barium application to baseline
* `Rt1`: initial input resistance calculated from inverting the slope of holding potential to current in MΩ
R_script: Astrocyte_MFA_Barium_R

#### `Microscopy.CSV`
* `Int_norm`: Normalized integrated density of Kir4.1/GFAP fluorescence, batch-normalized against dorsal control samples to minimize inter-staining variance.
* `LogGFAP`: $Log_{10}$ transformed values of raw GFAP Integrated Density.
* `Area_F`: Surface area percentage of Kir4.1-positive pixels overlapping with GFAP domains.
R_script: Immunohistochemistry_R
---

##  Methodological & Analytical Notes

* **Data Transformations:** Continuous variables like `slope_ratio` and `Blockade_ratio` were converted to their absolute values (`abs()`) to accurately reflect rectification indexes.
* **Distribution Handling:** Distributions were systematically verified via AIC modeling using the `fitdistrplus` package. Lognormal, Gamma, and Student-t distributions were preferentially utilized within Bayesian Generalized Linear Models (`brms`) where data significantly deviated from Gaussian normality.
* **Statistical Dependencies:** To reproduce the statistical reports, ensure the companion R scripts are run with the minimized package parameters (`tidyverse`, `car`, `permuco`, `brms`, `emmeans`, `fitdistrplus`, `dabestr`, `lubridate`, and `MASS`).

---

##  Raw Electrophysiology Data

Astrocytes were recorded in the whole-cell patch-clamp configuration. 

### Folder Hierarchy & Conditions
* **Baseline Condition (`MFA` folders):** Represents the baseline astrocytic current baseline control.
* **Experimental Manipulation (`MFA_Ba` folders):** Represents the active pharmacological application of Barium Chloride ($BaCl_2$) to block Kir4.1 channels.
* **Sub-folder Structure:** Within each region folder (`DH` and `VH`), files are grouped into separate `MFA` and `MFA_Ba` directories, further categorized by Zeitgeber Time point (**ZT3, ZT8, ZT15**).

### File Format & Naming Convention
* **Format:** Raw files are saved in standard **ATF (Axon Text File)** format, fully compatible with **pCLAMP software (v10.4)**.
* **Naming Metadata:** Individual file names explicitly denote the **experimental date (DDMMYYYY)** and the **cell order** recorded on that day, allowing seamless cross-referencing with the processed master data frames.

## 📜 Licensing & Usage Terms

This data is made available under the **Creative Commons Attribution 4.0 International License (CC-BY 4.0)**. You are free to share, copy, and adapt this data for any purpose, provided appropriate credit is given via citation of the paper or the dataset DOI listed above.
