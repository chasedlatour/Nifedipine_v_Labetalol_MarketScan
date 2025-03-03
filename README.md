This repository contains the SAS code and reference files used for the analysis by Latour et al. 2025 (forthcoming) to compare initiation of nifedipine versus labetalol prior to 23 weeks' gestation. 

The base set of pregnancies identified for this algorithm were identified using the code in the **prospective_pregnancy_marketscan** repository. 

All of the code contains substantial documentation, but I have provided some additional detail here to get started.

To start, run the following files to create the SAS datasets needed throughout the rest:
- `00_create_reference_files.sas` -- This creates SAS datasets based upon the codelists saved in the Excel file `Codelist.xlsx`.
- `01a_import_preg_data.sas` -- All of the relevant pregnancy files were saved in a different location on UNC's server, so I imported the necessary temporary files from that location into the current project folder via this code.
- `01b_clean_preg_cohort.sas` -- The initial pregnancy cleaning code was developed for this project, not the **prospective_pregnancy_marketscan** repository. As such, I implemented my own cleaning based upon the pregnancy cohort output after the `Step11___PregPlusLMP.sas` program, prior to implementing `Step12___clean_pregnancies.sas`.

I then derived the cohorts and necessary covariates using the following files:
- `02_derive_cohort.sas` -- This derived the cohorts for the primary analyses and sensitivity analyses that assumed different gestational ages at the indexing prenatal encounter.
- `02_derive_cohort_3wk.sas` -- This derived the cohorts for th esensitivity analyses that assumed that estimated LMPs for non-live births were +/- 3 weeks of those initially estimated by the algorithm.
- `03_derive_covaraites.sas` -- This derived the the versions of covariates that were implemented in analyses.

Analyses were conducted the the subsequent files:
- `04_primary_analyses.sas` -- Primary analyses.
- `04_primary_analyses_noweights.sas` -- Primary analyses, applying no weights and only applying standardized mortality ratio weights.
- `05_full_bounds_loss.sas` and `_prb.sas` -- Sensitivity analyses constructing the full bounds on the pregnancy loss and preterm birth outcomes if outcomes for all censored pregnancies were missing not at random.
- `06_enrl_bounds_loss.sas`, `_ptb.sas` and `_ptb_loss.sas` -- Sensitivity analyses leveraging enrollment information to identify pregnancies with outcomes missing not at random. The final program was used to construct those bounds if pregnancy loss was incorporated into the bounds for preterm birth.
- `07_complete_case.sas` -- Complete case analyses (i.e., only pregnancies with observed outcomes).
- `08_sens_12gwAtPNC.sas` and `_6gwAtPNC.sas` - Sensitivity analyses where we assumed that censored pregnancies with no gestational age information were 6 weeks' and 12 weeks' gestation at their indexing prenatal encounter.
- `09_sens_after2019.sas` - Sensitivity analyses restricted to pregnancies with indexing fills on or after January 1, 2019.
- `10_sens_dt_lmp.sas` - Sensitivity analyses where the lookback date for covariate information was indexing on the LMP for the pregnancies.
- `12_sens_chtn.sas` - Sensitivity analyses restricting to pregnancies with diagnosis codes for chronic hypertension.
- `13_sens_lmp_mins32.sas` and `_plus3w.sas` - Sensitivity analyses substracting and adding, respectively 21 days (3 weeks) from the estimated LMP for non-live births.
- `14_sens_uab.sas` - Sensitivity analyses assuming that unspecified abortions were induced abortions.
- `15_sens_censor_disenroll.sas` - Sensitivity analyses that censored pregnancies at their last prenatal encounter prior to their first disenrollment after the index date.
- `16_sens_covariates_pre.sas` - Sensitivity analyses where covariate information was identified only using claims prior to the index date.
- `17_sens_udl.sas` - Sensitivity analyses that assumed uncategorized deliveries were stillbirths.

Some other relevant files include:
- `FormatStatements_CDWH.sas` - SAS file with multiple format statements that may need to be run to view the initial pregnancy dataset.
- `bound_analyses.sas` -- Macro to conduct the full bound analyses.
- `combine_boot_estimates.sas` - Macro to combine estimates from the bootstrapped samples.
- `combine_point_estimates.sas` - Macro to combine the point estimates from the original sample.
- `competing2risk_weights.sas` - Macro to implement Aalen-Johansen estimator with 2 competing events. Includes application of standardized mortality ratio and inverse probability of censoring weights.
- `competing3risk_weights.sas` - Macro to implement Aalen-Johansen estimator with 3 competing events. Includes application of standardized mortality ratio and inverse probability of censoring weights.
- `count_missing_zero.sas` - Macro that counts the number of bootstrapped sample with estimated risks of 0 or missing values for estimates.
- `full_fup_weights.sas` - Macro that conducts weighted proportions to calculate effect estimates in a sample with full follow-up (e.g., complete case analyses).
- `getga.sas` - Macro used only for the file `01b_clean_preg_cohort.sas` to re-estimate LMPs under different assumptions about the gestational age at the claim date.
- `logparse.sas` - Macro used in setup.
- `overall_esetimates_w_ci.sas` - Macro that outputs the estimates in the overall sample with the 95% confidence interval based upon bootstrapped standard errors.
- `passinfo.sas` - Macro used in setup.
- `repeat_in_subset.sas` - Macro that repeats the analysis in a subset of the cohort.
- `setup.sas` - Macro run at the beginning of the program that sets necessary libraries and calls in macro files.
- `stddiff.sas` - Macro used to calculate standardized mean differences within levels of a covariate. Estimates across levels of the same covariate are incorrect.
- `stat_estimates_w_ci.sas` - Macro that returns the stratified risk, risk difference, and risk ratio estimates with 95% confidence intervals based upon bootstrapped standard errors.
- `table1.sas` - Macro that returns a weighted and unweighted set of descriptive characteristics fo the study population.

Data are accessible after payment to Merative with an appropriate data use agreement. All analyses were approved by UNC's Institutional Review Board. No data are uploaded to this repository.
