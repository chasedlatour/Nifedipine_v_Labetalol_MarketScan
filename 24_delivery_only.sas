/********************************************************************************************************************************************
PROGRAM: 24_delivery_only.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to repeat the primary analysis among those pregnancies that ended in delivery only.
	
Goal: To see if our results move closer to the null once we emulate the requirements of Leonard et al. 2025
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
%setup(sample=full, programname=24_delivery_only, savelog=Y);

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

%*Create the cohort like we did the primary cohort -- implement relevant inclusion and exclusion criteria.;

%implement_exclusions(inds=ana.preg_covar_63_dt_lmp_180_del, outds=delivery);

%*Require continuous enrollment (rx) from LMP-180 through the index date and then any disenrollment +28 days after delivery;
data delivery;
set delivery;

	%*Require continuous enrollment between lmp-180 and dt_index -- start and end centered around dt_index;
	if dt_lmp - 180 >= cont_enrl_start_ffsrx then pre_enrl = 1;
		else pre_enrl = 0;

	%*Require continuous enrollment between dt_index and dt_gapreg + 28 -- start and end centered around dt_index;
	if dt_gapreg + 28 <= cont_enrl_end_any then end_enrl = 1;
		else end_enrl = 0;

	%*Retain those that meet the enrollment requirements;
	if pre_enrl = 1 and end_enrl = 1 then output;

run;

%*Create a preterm birth and preeclampsia variable;

data preterm;
set delivery;

	*Preterm birth;
	if dt_gapreg - dt_lmp < 259 then preterm = 1;
		else preterm = 0;

	*First get the minimum date;
	dt_preeclampsia = min(dt_preec_outcInpt, dt_preec_outcOutpt);

	*Preeclampsia;
	if . < dt_preeclampsia < dt_gapreg + 14 then preeclampsia = 1;
		else preeclampsia = 0;

run;


%*For preeclampsia analysis, remove those with preeclampsia dxs prior to study entry.;
data preeclampsia;
set preterm;
	where preeclampsia_pre = 0;
run;

/*proc freq data=preeclampsia;*/
/*	table preeclampsia/missing;*/
/*run;*/
		

/********************************************************************************************************************************************

													02 - CONDUCT PRETERM BIRTH ANALYSIS

********************************************************************************************************************************************/

*Get counts -- overall;
proc freq data=preterm;
	table trt ga_index_cat*trt;
run;

*First, get the point estimate;
%full_fup_weights(
	boot = 0,
	inds = preterm,
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
	outds=ana.del_ptb_point,
	outds_dist=ga_del_ptb_point
);

/**Second, run the bootstrap to get standard error;*/
%full_fup_weights(
	boot = 1,
	inds = preterm,
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
	outds=ana.del_ptb_boot,
	outds_dist=ga_del_ptb_boot
);


options mlogic mprint symbolgen notes;

*Combine the point estimates;
%combine_point_estimates(
	gadist=ana.ga_del_ptb_point,
	numgastrat=2,
	est=ana.del_ptb_point_,
	outds=ana.del_ptb_point_overall
);

*Combine all the bootstrapped estimates;
%combine_boot_estimates(
	inputEst=ana.del_ptb_boot,
	inputDist=ana.ga_del_ptb_boot,
	numStrata=2,
	output_stratified=del_ptb_boot_strat,
	output_overall=del_ptb_boot_overall
);

*Stratified estimates with confidence intervals;
%strat_estimates_w_ci(
	bootdsn=ana.del_ptb_boot_1,
	pointdsn=ana.del_ptb_point_1,
	output=ana.del_ptb_1_ci
);
%strat_estimates_w_ci(
	bootdsn=ana.del_ptb_boot_2,
	pointdsn=ana.del_ptb_point_2,
	output=ana.del_ptb_2_ci
);

*overall estiamte with confidence interval;
%overall_estimates_w_ci(
	stderrdsn=del_ptb_boot_overall,
	pointdsn=ana.del_ptb_point_overall,
	output=ana.del_ptb_all_ci
);



%count_missing_zero(inds1=ana.del_ptb_boot_1, inds2=ana.del_ptb_boot_2);









/********************************************************************************************************************************************

													03 - CONDUCT PREECLAMPSIA ANALYSIS

********************************************************************************************************************************************/


*Get counts -- overall;
proc freq data=preeclampsia;
	table trt ga_index_cat*trt;
run;



*First, get the point estimate;
%full_fup_weights(
	boot = 0,
	inds = preeclampsia,
	gacatvar = ga_index_cat,
	outcomevar = preeclampsia,
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
	outds=ana.del_preec_point,
	outds_dist=ga_del_preec_point /*cca_del_preec_point*/
);


/**Second, run the bootstrap to get standard error;*/
%full_fup_weights(
	boot = 1,
	inds = preeclampsia,
	gacatvar = ga_index_cat,
	outcomevar = preeclampsia,
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
	outds=ana.del_preec_boot,
	outds_dist=ga_del_preec_boot
);

options mlogic mprint symbolgen notes;

*Combine the point estimates;
%combine_point_estimates(
	gadist=ana.ga_del_preec_point,
	numgastrat=2,
	est=ana.del_preec_point_,
	outds=ana.del_preec_point_overall
);

*Combine all the bootstrapped estimates;
%combine_boot_estimates(
	inputEst=ana.del_preec_boot,
	inputDist=ana.ga_del_preec_boot,
	numStrata=2,
	output_stratified=del_preec_boot_strat,
	output_overall=del_preec_boot_overall
);

*Stratified estimates with confidence intervals;
%strat_estimates_w_ci(
	bootdsn=ana.del_preec_boot_1,
	pointdsn=ana.del_preec_point_1,
	output=ana.del_preec_1_ci
);
%strat_estimates_w_ci(
	bootdsn=ana.del_preec_boot_2,
	pointdsn=ana.del_preec_point_2,
	output=ana.del_preec_2_ci
);

*overall estiamte with confidence interval;
%overall_estimates_w_ci(
	stderrdsn=del_preec_boot_overall,
	pointdsn=ana.del_preec_point_overall,
	output=ana.del_preec_all_ci
);



%count_missing_zero(inds1=ana.del_preec_boot_1, inds2=ana.del_preec_boot_2);

