# GALAH_DR3

This repository accompanies GALAH DR3 and its release paper.

If you are only looking for the catalog, please go here: [GALAH DR2 (not yet DR3)](https://datacentral.org.au/teamdata/GALAH/public/)

## Content of this repository:  
  
### 1) input  
This directory includes the information that WG4 has received from WG3 (Reduction) as well as additional information on performed crossmatches and external information (such as asteroseismic values). Here we also plot a first overview of the data:

![alt text](https://github.com/svenbuder/GALAH_DR3/tree/master/release_paper/figuresplot_parallax_quality_and_cmds.png "Astro- and photometric overview of the observed stars")

**Do not forget the t-SNE input for flagging!**
**figures not fully migrated yet!**

### 2) analysis

#### a) stellar parameters
- The IDL script used to run the stellar parameter estimation, [galah_sp.pro](https://github.com/svenbuder/GALAH_DR3/tree/master/analysis/stellar_parameters/galah_sp.pro)
- The final selection of lines used for the [DR3_Sp](https://github.com/svenbuder/GALAH_DR3/tree/master/analysis/stellar_parameters/DR3_Sp.dat) estimation, as well as the continuum segments, [DR3_Segm](https://github.com/svenbuder/GALAH_DR3/tree/master/analysis/stellar_parameters/DR3_Segm.dat).

#### b) abundance_estimation
- The IDL script used to run the element abundance estimation, [galah_ab.pro](https://github.com/svenbuder/GALAH_DR3/tree/master/analysis/abundances/galah_ab.pro)
- The final selection of masks (line center, line mask, continuum mask): [mode_DR3](https://github.com/svenbuder/GALAH_DR3/tree/master/analysis/abundances/mode_DR3)
- **convert_mode_dr3_to_tex.ipynb**
- The IDL script used to run the upper limit estimation for element abundance, [galah_limits.pro](https://github.com/svenbuder/GALAH_DR3/tree/master/analysis/abundances/galah_limits.pro)

#### c) collection routine
- The IDL script used to collect SME output files into FITS files, [galah_collect.pro](https://github.com/svenbuder/GALAH_DR3/tree/master/analysis/galah_collect.pro)

#### c) tests used to improve the setup
- Which lines should be used and which continuum regions should be selected. A test was performed for 10000 randomly chosen spectra, see [select_masks](https://github.com/svenbuder/GALAH_DR3/tree/master/analysis/abundances/tests/select_masks), see **plot_abundance_masks_for_10k_spectra.ipynb**


### 3) Processing

#### galah_output_structure.ipynb
- This notebook sets up the keywords for the *processed_files*, saved as *galah_dr3_output_structure.fits*

#### process_sme_results.ipynb
- This notebook takes the *sme_result_files* files and processes them, i.e. estimates final uncertainties, combines abundances measurements, and performs flagging. results are saved under *processed_files*.
- formerly known as create_data_products.ipynb

#### create_main_catalog.ipynb
- This notebook takes the *processed_files*, combines them, and post-processes them to create the final GALAH_DR3_main.fits catalog.
- formerly known as combine_10k_subsets.ipynb

**data_products not fully migrated yet!**
**figures not fully migrated yet!**

### 4) validation

**CHECK WHICH FILES FROM GALAH_iDR3 MUST BE COPIED**

#### a) repeat_observations
- estimate_repeat_uncertainty.ipynb

#### b) stellar_parameters
- GBS
- gaia_radial_velocity_differences.ipynb
- asteroseismology: see *input/asteroseismic_information/*

#### c) abundances
- **abundance_zeropoints.ipynb**: Notebook that includes the final decisions on abundance zeropoints
**Reformat for final format**
- abundance zeropoints: we have estimated abundances from a skyflat as well as arcturus
- We have performed tests around the best fit for both sun and arcturus (*SME_sun_arcturus_abundance_variations/* including SME setup to *test_all_abundance_variations.pro*)
- We have created an overview of the abundance zeropoints (*abundance_zeropoints_overviews/*)
- **Check again: What does Plot_line_by_line_abundances.ipynb do?**
- **Check again: What does abundance_differences_line_of_same_element.ipynb do? Plots histograms of differences between lines of same element**

#### d) flagging
- *flagging/binarity/birko_binarity.ipynb*: see Zwitter email from 200205 re. Birko et al. 2019 (AJ, 158, 155)

#### e) comparisons with other data
- Clusters: **comparison_clusters/** include HERMES_OC_Gaia_distances.fits from Spina (email 191205)
- GALAH DR2
- SDSS-IV / APOGEE DR16 **redo with DR16 abundances** see /validation/APOGEE_comparison/Comparison_GALAH_iDR3_APOGEE_DR16beta.ipynb
- LAMOST VAC by Xiang et al.
- Solar twins: **already copied spina/bedell. do plots**

### 4) catalogs
This directory includes either links to VAC descriptions or the direct scripts that were used to compute them.

### 5) spectra
This directory includes tools to plot spectra, identified via **sobject_id**. This includes both *observed* as well as *synthetic* spectra (if available).