/********************************************************************************************************************************************
PROGRAM: 07_complete_case.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to repeat the primary analysis among those pregnancies with observed outcomes.
	
Goal: 
Output data: 

Date: 
********************************************************************************************************************************************/






/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - PREPARE DATASET
	- 02 - CONDUCT LOSS ANALYSIS
	- 03 - CONDUCT PRETERM BIRTH ANALYSIS

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
%setup(sample=full, programname=07_complete_case, savelog=N);

******************************************************************************************************************************************;
/*Create local mirros of the libraries from the set up macro - Run locally*/
/*libname lout slibref=out server=server; */
/*libname lana slibref=ana server=server;*/
/*libname lwork slibref=work server=server;*/
/*libname ltemp slibref=temp server=server;*/


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



/********************************************************************************************************************************************

														01 - PREPARE DATASET

********************************************************************************************************************************************/

*Subset to complete cases where pregnancy outcome is unknown;
data cca;
set ana.primary_cohort;
	where preg_outcome_clean ne 'UNK';

	*Create a dichotomous variable for pregnancy loss;
	if preg_outcome_clean in ('SAB' 'UAB' 'IAB') then loss = 1;
		else loss = 0;

	*Create a dichotmous variable for preterm birth prior to 37 weeks gestation;
	if preg_outcome_ptb = "PTB" then preterm = 1;
		else preterm = 0;

run;

*Check for any unexpected missingness;
proc freq data=cca;
	table loss  preterm / missing;
run;

proc freq data=cca;
	table trt / missing;
run;





		

/********************************************************************************************************************************************

													02 - CONDUCT LOSS ANALYSIS

********************************************************************************************************************************************/

*First, get the point estimate;
%full_fup_weights(
	boot = 0,
	inds = cca,
	gacatvar = ga_index_cat,
	outcomevar = loss,
	psvars=ga_quartile age_at_index age_at_index_2
					year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
					year_index2017*t2dmrx_post diabetes_simp*t2dmrx_post
					nausea_pre recurlos_pre obesity_post chronichypertension_pre
					depressi_post anxiety_post antideprx_post benzorx_post 
					teratrx_pre num_outptpnc num_outptpnc_2,
	psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
					nausea_pre recurlos_pre obesity_post chronichypertension_pre 
					depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre,
	trtvar = trt,
	numiterations=1,
	initialseed=23244,
	outds=ana.cca_loss_point,
	outds_dist=ga_cca_loss_point
);

*Second, run the bootstrap to get standard error;
%full_fup_weights(
	boot = 1,
	inds = cca,
	gacatvar = ga_index_cat,
	outcomevar = loss,
	psvars=ga_quartile age_at_index age_at_index_2
					year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
					year_index2017*t2dmrx_post diabetes_simp*t2dmrx_post
					nausea_pre recurlos_pre obesity_post chronichypertension_pre
					depressi_post anxiety_post antideprx_post benzorx_post 
					teratrx_pre num_outptpnc num_outptpnc_2,
	psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
					nausea_pre recurlos_pre obesity_post chronichypertension_pre 
					depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre,
	trtvar = trt,
	numiterations=1000,
	initialseed=23244,
	outds=ana.cca_loss_boot,
	outds_dist=ga_cca_loss_boot
);

options mlogic mprint symbolgen notes;

*Combine the point estimates;
%combine_point_estimates(
	gadist=ana.ga_cca_loss_point,
	numgastrat=2,
	est=ana.cca_loss_point_,
	outds=cca_loss_point_overall
);

*Combine all the bootstrapped estimates;
%combine_boot_estimates(
	inputEst=ana.cca_loss_boot,
	inputDist=ana.ga_cca_loss_boot,
	numStrata=2,
	output_stratified=cca_loss_boot_strat,
	output_overall=cca_loss_boot_overall
);

*Stratified estimates with confidence intervals;
%strat_estimates_w_ci(
	bootdsn=ana.cca_loss_boot_1,
	pointdsn=ana.cca_loss_point_1,
	output=ana.cca_loss_1_ci
);
%strat_estimates_w_ci(
	bootdsn=ana.cca_loss_boot_2,
	pointdsn=ana.cca_loss_point_2,
	output=ana.cca_loss_2_ci
);

*overall estiamte with confidence interval;
%overall_estimates_w_ci(
	stderrdsn=cca_loss_boot_overall,
	pointdsn=cca_loss_point_overall,
	output=ana.cca_loss_all_ci
);



%count_missing_zero(inds1=ana.cca_loss_boot_1, inds2=ana.cca_loss_boot_2);









/********************************************************************************************************************************************

													03 - CONDUCT PRETERM BIRTH ANALYSIS

********************************************************************************************************************************************/


*First, get the point estimate;
%full_fup_weights(
	boot = 0,
	inds = cca,
	gacatvar = ga_index_cat,
	outcomevar = preterm,
	psvars=ga_quartile age_at_index age_at_index_2
					year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
					nausea_pre recurlos_pre obesity_post chronichypertension_pre
					depressi_post anxiety_post antideprx_post benzorx_post 
					teratrx_pre num_outptpnc num_outptpnc_2,
	psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
					nausea_pre recurlos_pre obesity_post chronichypertension_pre 
					depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre,
	trtvar = trt,
	numiterations=1,
	initialseed=23244,
	outds=ana.cca_ptb_point,
	outds_dist=ga_cca_ptb_point
);

*Second, run the bootstrap to get standard error;
%full_fup_weights(
	boot = 1,
	inds = cca,
	gacatvar = ga_index_cat,
	outcomevar = preterm,
	psvars=ga_quartile age_at_index age_at_index_2
					year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
					nausea_pre recurlos_pre obesity_post chronichypertension_pre
					depressi_post anxiety_post antideprx_post benzorx_post 
					teratrx_pre num_outptpnc num_outptpnc_2,
	psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
					nausea_pre recurlos_pre obesity_post chronichypertension_pre 
					depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre,
	trtvar = trt,
	numiterations=1000,
	initialseed=23244,
	outds=ana.cca_ptb_boot,
	outds_dist=ga_cca_ptb_boot
);

options mlogic mprint symbolgen notes;

*Combine the point estimates;
%combine_point_estimates(
	gadist=ana.ga_cca_ptb_point,
	numgastrat=2,
	est=ana.cca_ptb_point_,
	outds=cca_ptb_point_overall
);

*Combine all the bootstrapped estimates;
%combine_boot_estimates(
	inputEst=ana.cca_ptb_boot,
	inputDist=ana.ga_cca_ptb_boot,
	numStrata=2,
	output_stratified=cca_ptb_boot_strat,
	output_overall=cca_ptb_boot_overall
);

*Stratified estimates with confidence intervals;
%strat_estimates_w_ci(
	bootdsn=ana.cca_ptb_boot_1,
	pointdsn=ana.cca_ptb_point_1,
	output=ana.cca_ptb_1_ci
);
%strat_estimates_w_ci(
	bootdsn=ana.cca_ptb_boot_2,
	pointdsn=ana.cca_ptb_point_2,
	output=ana.cca_ptb_2_ci
);

*overall estiamte with confidence interval;
%overall_estimates_w_ci(
	stderrdsn=cca_ptb_boot_overall,
	pointdsn=cca_ptb_point_overall,
	output=ana.cca_ptb_all_ci
);


%count_missing_zero(inds1=ana.cca_ptb_boot_1, inds2=ana.cca_ptb_boot_2);
