/********************************************************************************************************************************************
PROGRAM: 04_primary_analyses.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to clean the pregnancy cohort that we derived from the MarketScan claims data
and prepare it for analyses.
	
Goal: 
Output data: 

Date: 12.19.2024
********************************************************************************************************************************************/






/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 01 - CONDUCT PREGNANCY LOSS ANALYSIS
	- 04 - CONDUCT PRETERM BIRTH ANALYSIS
	- 05 - AD-HOC DESCRIPTIVES

********************************************************************************************************************************************/











/********************************************************************************************************************************************

															00 - SET UP LIBRARIES

********************************************************************************************************************************************/

/*run local:*/
/*SIGNOFF;*/
/*%LET server=n2.schsr.unc.edu 1234; */
/*options comamid=tcp remote=server; */
/*signon username=_prompt_;*/


*Run setup macro and define libnames;

options sasautos=(SASAUTOS "/local/projects/marketscan_preg/Latour_23_2322/programs/macros");
/*options mprint;*/

/*change "saveLog=" to "Y" when program is closer to complete*/
%setup(sample=full, programname=25_primary_more_cov_loss, savelog=Y);

******************************************************************************************************************************************;
/*Create local mirros of the libraries from the set up macro - Run locally*/
/*libname lout slibref=out server=server; */
/*libname lana slibref=ana server=server;*/
/*libname lwork slibref=work server=server;*/
/*libname ltemp slibref=temp server=server;*/
/*libname lexpref slibref=expref server=server;*/
/*libname lcovref slibref=covref server=server;*/


*Create formats;
proc format;
	VALUE treated
		2 = "Diagnosis, medication fill"
		1 = "Daignosis, no fill"
		0 = "No diagnosis";
	VALUE $pregoutc
		"IAB" = "Induced abortion"
		"LBM" = "Live birth, Multiple"
		"LBS" = "Live birth, Singleton"
		"SAB" = "Spontaneous abortion"
		"SB" = "Stillbirth"
		"UAB" = "Unspecified abortion"
		"UDL" = "Uncategorized delivery"
		"UNK" = "Unknown outcome";
	VALUE $diabetes
		"EOND" = "Early-onset GDM or newly diagnosed"
		"GDM" = "Gestational diabetes"
		"NA" = "No diabetes"
		"NOS" = "Pregestational diabetes, not otherwise specified"
		"T1DM" = "Type-1 diabetes"
		"T2DM" = "Type-2 diabetes"
		"UNSP" = "Diabetes, type unknown";
	VALUE $simpdiab
		"PRE" = "Pregestational diabetes"
		"OTH" = "Other diabetes"
		"NA" = "No diabetes";
	VALUE unknown
		. = "Unknown"
		0 = "Metropolitan area"
		1 = "Rural area";
	/*CDL: ADDED region 1.5.2026*/
	value $region
		"1" = "Northeast"
		"2" = "North Central"
		"3" = "South"
		"4" = "West"
		"5" = "Unknown";
run;


%*13 people with region unknown. Lump these into the largest category;
data primary;
set ana.primary_cohort;
	where region ne '5';
run;




	
/********************************************************************************************************************************************

														01 - CONDUCT PREGNANCY LOSS ANALYSIS

********************************************************************************************************************************************/
	
		
*Get the point estimate;
%competing2risk_weights(
	boot=0, 
	inds=primary, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_clean,
	eventDT=dt_GApreg, 
	event = 'SAB' 'UAB' 'SB' 'MLS',
	cr1='IAB', 
	cr2='LBM' 'LBS' 'UDL',
	psvars=ga_quartile age_at_index age_at_index_2 year_index2017 region
			year_index2017*region t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post num_outptpnc num_outptpnc_2
			fibroids_post ldarx_pre nsaidrx_pre,
	dovars=ga_quartile age_at_index year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre depressi_post anxiety_post antideprx_post 
			benzorx_post teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	dovarsmodel = ga_quartile age_at_index year_index4
			t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post
			num_outptpnc num_outptpnc_2 rural2,
	psclassvars=ga_quartile year_index2017 region t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post teratrx_pre
			fibroids_post ldarx_pre nsaidrx_pre, 
	doclassvars=ga_quartile year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre rural2,
	trtvar=trt,
	numiterations=1,
	initialseed=23244, 
	outds=ana.moreadj,
	outds_dist=ga_dist_moreadj,
	outds_ps = NA,
	outds_surv = NA
	);


*Look at the stratified estimates

Datasets:
- ana.primary_point_noipcw_1
- ana.primary_point_noipcw_2;

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_moreadj,
	numgastrat = 2,
	est = ana.moreadj_,
	outds = ana.moreadj_overall
	);

	

*********Finally, conduct the bootstrap;
	
*Now incorporate the IPCW;
%competing2risk_weights(
	boot=1, 
	inds=primary, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_clean,
	eventDT=dt_GApreg, 
	event = 'SAB' 'UAB' 'SB' 'MLS',
	cr1='IAB', 
	cr2='LBM' 'LBS' 'UDL',
	psvars=ga_quartile age_at_index age_at_index_2 year_index2017 region
			year_index2017*region t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post num_outptpnc num_outptpnc_2
			fibroids_post ldarx_pre nsaidrx_pre,
	dovars=ga_quartile age_at_index year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre depressi_post anxiety_post antideprx_post 
			benzorx_post teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	dovarsmodel = ga_quartile age_at_index year_index4
			t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post
			num_outptpnc num_outptpnc_2 rural2,
	psclassvars=ga_quartile year_index2017 region t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post teratrx_pre
			fibroids_post ldarx_pre nsaidrx_pre, 
	doclassvars=ga_quartile year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre rural2,
	trtvar=trt, 
	numiterations=1000,
	initialseed=23244, 
	outds=ana.moreadj_boot,
	outds_dist=ga_dist_moreadj_boot,
	outds_ps = NA,
	outds_surv = NA
	);
	
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.moreadj_boot, 
			inputDist= ana.ga_dist_moreadj_boot,
			numStrata= 2, 
			output_stratified= ana.moreadj_boot_strat,
			output_overall= ana.moreadj_boot_overall);


*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.moreadj_boot_1, pointdsn=ana.moreadj_1, output=ana.moreadj_boot_1_ci);
%strat_estimates_w_CI(bootdsn=ana.moreadj_boot_2, pointdsn=ana.moreadj_2, output=ana.moreadj_boot_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.moreadj_boot_overall, pointdsn=ana.moreadj_overall, output=ana.moreadj_overall_boot_ci);
	
	
%count_missing_zero(inds1=ana.moreadj_boot_1, inds2=ana.moreadj_boot_2);








