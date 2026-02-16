/********************************************************************************************************************************************
PROGRAM: 04_preec_outcome.sas
PROGRAMMER: Chase Latour
PURPOSE: Conduct post-hoc analysis where preeclampsia (up to 2 weeks post-pregnancy outcome) is the study outcome
	
Goal: 
Output data: 

Date: 3.29.25
********************************************************************************************************************************************/






/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - IMPLEMENT EXCLUSION CRITERIA AND GET RELEVANT COUNTS
	- 02 - INVESTIGATE COVARIATE DISTRIBUTIONS
	- 03 - CONDUCT PREGNANCY LOSS ANALYSIS
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
%setup(sample=full, programname=04_preec_outcome, savelog=N);

******************************************************************************************************************************************;
/*Create local mirros of the libraries from the set up macro - Run locally*/
/*libname lout slibref=out server=server; */
/*libname lana slibref=ana server=server;*/
/*libname lwork slibref=work server=server;*/
/*libname ltemp slibref=temp server=server;*/
/*libname lexpref slibref=expref server=server;*/


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

										01 - IMPLEMENT EXCLUSION CRITERIA AND GET RELEVANT COUNTS

********************************************************************************************************************************************/

*First: Exclude individuals with prior preeclampsia diagnosis codes upon cohort entry;
data no_preeclampsia;
set ana.primary_cohort;
	where preeclampsia_pre = 0;
run;

/*proc contents data=no_preeclampsia; run;*/

*Second: Create outcome indicators.

In this case, people will be censored at their first disenrollment from a MarketScan plan;
data ana.preeclampsia_outc;
set no_preeclampsia;

	***Create inpatient preeclampsia outcome;

	*First, create the dates ignoring censoring;
	*if preeclampsia occurs prior to their last prenatal encounter, then they have preeclampsia;
	if preg_outcome_clean = "UNK" and . < dt_preec_outcInpt < dt_gapreg then do;
		preec_inpt_nocensor = "PRE";
		preec_inpt_dt_nocensor = dt_preec_outcInpt;
	end;
		*Otherwise, if preeclampsia does not occur prior to the last prenatal encounter, then censoredf;
		else if preg_outcome_clean = "UNK" then do;
			preec_inpt_nocensor = "UNK";
			preec_inpt_dt_nocensor = dt_gapreg;
		end;
		*otherwise, if they have preeclampsia prior to 14 days aftre their pregnancy outcome, then they have preeclampsia.;
		else if . < dt_preec_outcInpt < dt_gapreg + 14 then do;
			preec_inpt_nocensor = "PRE";
			preec_inpt_dt_nocensor = dt_preec_outcInpt;
		end;
		*otherwise, they do not have an outcome.;
		else do;
			preec_inpt_nocensor = preg_outcome_clean;
			preec_inpt_dt_nocensor = dt_gapreg + 14;
		end;

	*Now implement the censoring due to disenrollment.;
	*If censored before their preeclampsia outcome, then censored. otherwise it is the same;
	if . < cont_enrl_end_any < preec_inpt_dt_nocensor then do;
		preec_inpt_outc = "UNK";
		preec_inpt_dt = cont_enrl_end_any;
	end;
		else do;
			preec_inpt_outc = preec_inpt_nocensor;
			preec_inpt_dt = preec_inpt_dt_nocensor;
		end;


	***Create any preeclampsia outcome;

	*First get the minimum date;
	dt_preeclampsia = min(dt_preec_outcInpt, dt_preec_outcOutpt);

	*First, create the dates ignoring censoring;
	*if preeclampsia occurs prior to their last prenatal encounter, then they have preeclampsia;
	if preg_outcome_clean = "UNK" and . < dt_preeclampsia < dt_gapreg then do;
		preec_any_nocensor = "PRE";
		preec_any_dt_nocensor = dt_preeclampsia;
	end;
		*Otherwise, if preeclampsia does not occur prior to the last prenatal encounter, then censoredf;
		else if preg_outcome_clean = "UNK" then do;
			preec_any_nocensor = "UNK";
			preec_any_dt_nocensor = dt_gapreg;
		end;
		*otherwise, if they have preeclampsia prior to 14 days aftre their pregnancy outcome, then they have preeclampsia.;
		else if . < dt_preeclampsia < dt_gapreg + 14 then do;
			preec_any_nocensor = "PRE";
			preec_any_dt_nocensor = dt_preeclampsia;
		end;
		*otherwise, they do not have an outcome.;
		else do;
			preec_any_nocensor = preg_outcome_clean;
			preec_any_dt_nocensor = dt_gapreg + 14;
		end;

	*Now implement the censoring due to disenrollment.;
	*If censored before their preeclampsia outcome, then censored. otherwise it is the same;
	if . < cont_enrl_end_any < preec_any_dt_nocensor then do;
		preec_any_outc = "UNK";
		preec_any_dt = cont_enrl_end_any;
	end;
		else do;
			preec_any_outc = preec_any_nocensor;
			preec_any_dt = preec_any_dt_nocensor;
		end;

run;


/*proc freq data=ana.preeclampsia_outc;*/
/*	table preec_inpt_outc preec_any_outc / missing;*/
/*run;*/
/**/
/*proc means data=ana.preeclampsia_outc nmiss;*/
/*	var preec_inpt_dt preec_any_dt;*/
/*run;*/

	



	
/********************************************************************************************************************************************

													03 - ANY PREECLAMPSIA ANALYSIS

Challenges to defining preeclampsia:
- Outpatient encounters could represent rule-out diagnoses.
- INpatient encounters may represent more severe cases but could just be deliveries with mild preeclampsia.

Going to do both and see if particularly concerned about it. 
Date of capture is important in time-to-event analyses.

********************************************************************************************************************************************/
	

	
*Now incorporate the IPCW;
%competing3risk_weights(
	boot=0, 
	inds=ana.preeclampsia_outc, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preec_any_outc,
	eventDT=preec_any_dt, 
	event = 'PRE',
	cr1='IAB', 
	cr2='SAB' 'UAB' 'SB' 'MLS',
	cr3='LBS' 'UDL' 'LBM',
	psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars=ga_quartile age_at_index year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre depressi_post anxiety_post antideprx_post 
			benzorx_post teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	dovarsmodel = ga_quartile age_at_index year_index4
			t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post
			teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre, 
	doclassvars=ga_quartile year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre rural2,
	trtvar=trt,
	numiterations=1,
	initialseed=23244, 
	outds=ana.preec_any_point,
	outds_dist=ga_dist_preec_any_point,
	outds_ps = ps_preec_any_point,
	outds_surv = ana.surv_preec_any
	);


*Put out the trimmed cohorts for at-risk counts in cumulative incidence figure;
data ana.preeclampsia_any_trim_1;
format preec_any_dt MMDDYY10. preec_inpt_dt MMDDYY10.;
set _anacohort_trim1;
run;
data ana.preeclampsia_any_trim_2;
format preec_any_dt MMDDYY10. preec_inpt_dt MMDDYY10.;
set _anacohort_trim2;
run;

*Get counts for disenrollment;
data enroll_preec;
set ana.preeclampsia_outc;
	
	days_disenroll = dt_disenroll_post_any - preec_any_dt;
	daysg31 = days_disenroll > 31;

	ga_at_end = (preec_any_dt - dt_lmp) / 7;
	ga_at_index = (dt_index - dt_lmp) / 7;

run;

proc sort data=enroll_preec;
	by trt;
run;


*Overall cohort;
proc freq data=enroll_preec;
	by trt;
	table daysg31*preec_any_outc;
run;

*GA at index: <14 w;
proc freq data=enroll_preec (where = (ga_index_cat = 1));
	by trt;
	table daysg31*preec_any_outc;
run;

*GA at index: >=14 w;
proc freq data=enroll_preec (where = (ga_index_cat = 2));
	by trt;
	table daysg31*preec_any_outc;
run;

****GA stats;

*Overall cohort;
proc means data=enroll_preec (where = (preec_any_outc = 'UNK')) median p25 p75;
	class trt daysg31;
	var ga_at_index ga_at_end;
run;

*Stratified by gestational age at index;

*Overall cohort;
proc means data=enroll_preec (where = (preec_any_outc = 'UNK')) median p25 p75;
	class trt ga_index_cat daysg31;
	var ga_at_index ga_at_end;
run;




*Look at the stratified estimates;

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_preec_any_point,
	numgastrat = 2,
	est = ana.preec_any_point_,
	outds = ana.preec_any_point_overall
	);

	


*********Finally, conduct the bootstrap;
	
*Now incorporate the IPCW;
%competing3risk_weights(
	boot=1, 
	inds=ana.preeclampsia_outc, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preec_any_outc,
	eventDT=preec_any_dt, 
	event = 'PRE',
	cr1='IAB', 
	cr2='SAB' 'UAB' 'SB' 'MLS',
	cr3='LBS' 'UDL' 'LBM',
	psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars=ga_quartile age_at_index year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre depressi_post anxiety_post antideprx_post 
			benzorx_post teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	dovarsmodel = ga_quartile age_at_index year_index4
			t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post
			teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre, 
	doclassvars=ga_quartile year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre rural2,
	trtvar=trt,
	numiterations=1000,
	initialseed=23244, 
	outds=ana.preec_any_boot,
	outds_dist=ga_dist_preec_any_boot,
	outds_ps = NA,
	outds_surv = NA
	);
	
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.preec_any_boot,
			inputDist= ana.ga_dist_preec_any_boot,
			numStrata= 2, 
			output_stratified= ana.preec_any_boot_strat,
			output_overall= ana.preec_any_boot_overall);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.preec_any_boot_1, pointdsn=ana.preec_any_point_1, output=ana.preec_any_boot_1_ci);
%strat_estimates_w_CI(bootdsn=ana.preec_any_boot_2, pointdsn=ana.preec_any_point_2, output=ana.preec_any_boot_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.preec_any_boot_overall, pointdsn=ana.preec_any_point_overall, output=ana.preec_any_boot_ci);

	
%count_missing_zero(inds1=ana.preec_any_boot_1, inds2=ana.preec_any_boot_2);











/********************************************************************************************************************************************

												03 - INPATIENT PREECLAMPSIA ANALYSIS

Challenges to defining preeclampsia:
- Outpatient encounters could represent rule-out diagnoses.
- INpatient encounters may represent more severe cases but could just be deliveries with mild preeclampsia.

This section is only going to focus on the inpatient records for preeclampsia.

********************************************************************************************************************************************/
	

	
*Now incorporate the IPCW;
%competing3risk_weights(
	boot=0, 
	inds=ana.preeclampsia_outc, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preec_inpt_outc,
	eventDT=preec_inpt_dt, 
	event = 'PRE',
	cr1='IAB', 
	cr2='SAB' 'UAB' 'SB' 'MLS',
	cr3='LBS' 'UDL' 'LBM',
	psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars=ga_quartile age_at_index year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre depressi_post anxiety_post antideprx_post 
			benzorx_post teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	dovarsmodel = ga_quartile age_at_index year_index4
			t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post
			teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre, 
	doclassvars=ga_quartile year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre rural2,
	trtvar=trt,
	numiterations=1,
	initialseed=23244, 
	outds=ana.preec_inpt_point,
	outds_dist=ga_dist_preec_inpt_point,
	outds_ps = ps_preec_inpt_point,
	outds_surv = ana.surv_preec_inpt
	);


*Look at the stratified estimates;

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_preec_inpt_point,
	numgastrat = 2,
	est = ana.preec_inpt_point_,
	outds = ana.preec_inpt_point_overall
	);

	


*********Finally, conduct the bootstrap;
	
*Now incorporate the IPCW;
%competing3risk_weights(
	boot=1, 
	inds=ana.preeclampsia_outc, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preec_inpt_outc,
	eventDT=preec_inpt_dt, 
	event = 'PRE',
	cr1='IAB', 
	cr2='SAB' 'UAB' 'SB' 'MLS',
	cr3='LBS' 'UDL' 'LBM',
	psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars=ga_quartile age_at_index year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre depressi_post anxiety_post antideprx_post 
			benzorx_post teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	dovarsmodel = ga_quartile age_at_index year_index4
			t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post
			teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre, 
	doclassvars=ga_quartile year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre rural2,
	trtvar=trt,
	numiterations=1000,
	initialseed=23244, 
	outds=ana.preec_inpt_boot,
	outds_dist=ga_dist_preec_inpt_boot,
	outds_ps = NA,
	outds_surv = NA
	);
	
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.preec_inpt_boot,
			inputDist= ana.ga_dist_preec_inpt_boot,
			numStrata= 2, 
			output_stratified= ana.preec_inpt_boot_strat,
			output_overall= ana.preec_inpt_boot_overall);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.preec_inpt_boot_1, pointdsn=ana.preec_inpt_point_1, output=ana.preec_inpt_boot_1_ci);
%strat_estimates_w_CI(bootdsn=ana.preec_inpt_boot_2, pointdsn=ana.preec_inpt_point_2, output=ana.preec_inpt_boot_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.preec_inpt_boot_overall, pointdsn=ana.preec_inpt_point_overall, output=ana.preec_inpt_boot_ci);

	
%count_missing_zero(inds1=ana.preec_inpt_boot_1, inds2=ana.preec_inpt_boot_2);
