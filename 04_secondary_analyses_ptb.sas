/********************************************************************************************************************************************
PROGRAM: 04_secondary_analyses.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to conduct the secondary analyses for the study.

(1) Preterm birth at <28, <32, and <34 weeks of gestation
(2) Pregnancy loss at <12, 12-20, and >20 weeks of gestation
	
Goal: 
Output data: 

Date: 12.19.2024
********************************************************************************************************************************************/






/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - CONDUCT PRETERM BIRTH ANALYSIS

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
%setup(sample=full, programname=04_secondary_analyses_ptb, savelog=Y);

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

												01 - CONDUCT PRETERM BIRTH ANALYSIS

********************************************************************************************************************************************/


*Create preterm birth variables based on the specific GAs;
data preterm;
set ana.primary_cohort;

	******Preterm prior to 28 weeks;

	*Event indicator;
	if preg_outcome_clean in ('LBM' 'LBS' 'UDL') and dt_gapreg - dt_lmp < 196 then preg_outcome_ptb28 = 'PTB';
		else if dt_gapreg - dt_lmp >= 196 then preg_outcome_ptb28 = 'TB';
		else preg_outcome_ptb28 = preg_outcome_clean;

	*Outcome time;
	if preg_outcome_ptb28 = 'TB' then dt_gapreg_ptb28 = dt_lmp + 195;
		else dt_gapreg_ptb28 = dt_gapreg;

	*******Preterm prior to 32 weeks;

	*Event indicator;
	if preg_outcome_clean in ('LBM' 'LBS' 'UDL') and dt_gapreg - dt_lmp < 224 then preg_outcome_ptb32 = 'PTB';
		else if dt_gapreg - dt_lmp >= 224 then preg_outcome_ptb32 = 'TB';
		else preg_outcome_ptb32 = preg_outcome_clean;

	*Outcome time;
	if preg_outcome_ptb32 = 'TB' then dt_gapreg_ptb32 = dt_lmp + 223;
		else dt_gapreg_ptb32 = dt_gapreg;

	*******Preterm prior to 34 weeks;

	*Event indicator;
	if preg_outcome_clean in ('LBM' 'LBS' 'UDL') and dt_gapreg - dt_lmp < 238 then preg_outcome_ptb34 = 'PTB';
		else if dt_gapreg - dt_lmp >= 238 then preg_outcome_ptb34 = 'TB';
		else preg_outcome_ptb34 = preg_outcome_clean;

	*Outcome time;
	if preg_outcome_ptb34 = 'TB' then dt_gapreg_ptb34 = dt_lmp + 237;
		else dt_gapreg_ptb34 = dt_gapreg;

run;


*********PRETERM BIRTH PRIOR TO 28 WEEKS GESTATION;

*Point estimates;	
%competing3risk_weights(
	boot=0, 
	inds=preterm, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_ptb28,
	eventDT=dt_GApreg_ptb28, 
	event = 'PTB',
	cr1='IAB', 
	cr2='SAB' 'UAB' 'SB' 'MLS',
	cr3='TB',
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
	outds=ana.second_point_ptb28,
	outds_dist=ga_dist_second_point_ptb28,
	outds_ps = ps_second_point_ptb28,
	outds_surv = NA
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_second_point_ptb28,
	numgastrat = 2,
	est = ana.second_point_ptb28_,
	outds = ana.second_point_ptb28_overall
	);
	
*Conduct the bootstrap to get the variance;
%competing3risk_weights(
	boot=1, 
	inds=preterm, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_ptb28,
	eventDT=dt_GApreg_ptb28, 
	event = 'PTB',
	cr1='IAB', 
	cr2='SAB' 'UAB' 'SB' 'MLS',
	cr3='TB',
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
	outds=ana.second_boot_ptb28,
	outds_dist=ga_dist_second_boot_ptb28,
	outds_ps = ps_second_boot_ptb28,
	outds_surv = NA
	);
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.second_boot_ptb28, 
			inputDist= ana.ga_dist_second_boot_ptb28,
			numStrata= 2, 
			output_stratified= ana.second_ptb28_boot_strat,
			output_overall= ana.second_ptb28_boot_overall);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.second_boot_ptb28_1, pointdsn=ana.second_point_ptb28_1, output=ana.second_boot_ptb28_1_ci);
%strat_estimates_w_CI(bootdsn=ana.second_boot_ptb28_2, pointdsn=ana.second_point_ptb28_2, output=ana.second_boot_ptb28_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.second_ptb28_boot_overall, pointdsn=ana.second_point_ptb28_overall, output=ana.second_overall_ptb28_boot_ci);

	
*Count the number with zero risks, or missing effect esitmates;
%count_missing_zero(inds1=ana.second_boot_ptb28_1, inds2=ana.second_boot_ptb28_2);









*********PRETERM BIRTH PRIOR TO 32 WEEKS GESTATION;

*Point estimates;	
%competing3risk_weights(
	boot=0, 
	inds=preterm, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_ptb32,
	eventDT=dt_GApreg_ptb32, 
	event = 'PTB',
	cr1='IAB', 
	cr2='SAB' 'UAB' 'SB' 'MLS',
	cr3='TB',
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
	outds=ana.second_point_ptb32,
	outds_dist=ga_dist_second_point_ptb32,
	outds_ps = ps_second_point_ptb32,
	outds_surv = NA
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_second_point_ptb32,
	numgastrat = 2,
	est = ana.second_point_ptb32_,
	outds = ana.second_point_ptb32_overall
	);
	
*Conduct the bootstrap to get the variance;
%competing3risk_weights(
	boot=1, 
	inds=preterm, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_ptb32,
	eventDT=dt_GApreg_ptb32, 
	event = 'PTB',
	cr1='IAB', 
	cr2='SAB' 'UAB' 'SB' 'MLS',
	cr3='TB',
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
	outds=ana.second_boot_ptb32,
	outds_dist=ga_dist_second_boot_ptb32,
	outds_ps = ps_second_boot_ptb32,
	outds_surv = NA
	);
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.second_boot_ptb32, 
			inputDist= ana.ga_dist_second_boot_ptb32,
			numStrata= 2, 
			output_stratified= ana.second_ptb32_boot_strat,
			output_overall= ana.second_ptb32_boot_overall);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.second_boot_ptb32_1, pointdsn=ana.second_point_ptb32_1, output=ana.second_boot_ptb32_1_ci);
%strat_estimates_w_CI(bootdsn=ana.second_boot_ptb32_2, pointdsn=ana.second_point_ptb32_2, output=ana.second_boot_ptb32_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.second_ptb32_boot_overall, pointdsn=ana.second_point_ptb32_overall, output=ana.second_overall_ptb32_boot_ci);

%count_missing_zero(inds1=ana.second_boot_ptb32_1, inds2=ana.second_boot_ptb32_2);






*********PRETERM BIRTH PRIOR TO 34 WEEKS GESTATION;

*Point estimates;	
%competing3risk_weights(
	boot=0, 
	inds=preterm, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_ptb34,
	eventDT=dt_GApreg_ptb34, 
	event = 'PTB',
	cr1='IAB', 
	cr2='SAB' 'UAB' 'SB' 'MLS',
	cr3='TB',
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
	outds=ana.second_point_ptb34,
	outds_dist=ga_dist_second_point_ptb34,
	outds_ps = ps_second_point_ptb34,
	outds_surv = NA
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_second_point_ptb34,
	numgastrat = 2,
	est = ana.second_point_ptb34_,
	outds = ana.second_point_ptb34_overall
	);
	
*Conduct the bootstrap to get the variance;
%competing3risk_weights(
	boot=1, 
	inds=preterm, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_ptb34,
	eventDT=dt_GApreg_ptb34, 
	event = 'PTB',
	cr1='IAB', 
	cr2='SAB' 'UAB' 'SB' 'MLS',
	cr3='TB',
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
	outds=ana.second_boot_ptb34,
	outds_dist=ga_dist_second_boot_ptb34,
	outds_ps = ps_second_boot_ptb34,
	outds_surv = NA
	);
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.second_boot_ptb34, 
			inputDist= ana.ga_dist_second_boot_ptb34,
			numStrata= 2, 
			output_stratified= ana.second_ptb34_boot_strat,
			output_overall= ana.second_ptb34_boot_overall);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.second_boot_ptb34_1, pointdsn=ana.second_point_ptb34_1, output=ana.second_boot_ptb34_1_ci);
%strat_estimates_w_CI(bootdsn=ana.second_boot_ptb34_2, pointdsn=ana.second_point_ptb34_2, output=ana.second_boot_ptb34_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.second_ptb34_boot_overall, pointdsn=ana.second_point_ptb34_overall, output=ana.second_overall_ptb34_boot_ci);


%count_missing_zero(inds1=ana.second_boot_ptb34_1, inds2=ana.second_boot_ptb34_2);









