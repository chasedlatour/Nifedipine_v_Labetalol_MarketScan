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
	- 02 - DERIVE PREGNANCY LOSS DATA AND DESCRIPTIVES
	- 03 - CONDUCT PREGNANCY LOSS ANALYSIS

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
%setup(sample=full, programname=04_secondary_analyses, savelog=Y);

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











/********************************************************************************************************************************************

											02 - DERIVE PREGNANCY LOSS DATA AND DESCRIPTIVES

********************************************************************************************************************************************/

*Subset to the cohort of pregnancies to be included in the analysis -- those who enrolled prior to 12 weeks of gestation;
data loss;
set ana.primary_cohort;
	where dt_index - dt_lmp < 84; 

	*******LOSS PRIOR TO 12 WEEKS;

	*Event indicator;
	if preg_outcome_clean in ('SAB' 'UAB' 'SB' 'MLS') and dt_gapreg - dt_lmp < 84 then preg_outcome_loss12 = 'OUT';
		else if dt_gapreg - dt_lmp >= 84 then preg_outcome_loss12 = 'NO';
		else preg_outcome_loss12 = preg_outcome_clean;

	*Event time variable;
	if preg_outcome_loss12 = 'NO' then dt_gapreg_loss12 = dt_lmp + 83;
		else dt_gapreg_loss12 = dt_gapreg;

	*******LOSS 12 TO 20 WEEKS;

	*Event indicator;
	if preg_outcome_clean in ('SAB' 'UAB' 'SB' 'MLS') and dt_gapreg - dt_lmp < 82 then preg_outcome_loss1220 = 'EAR'; /*early losses - will be competing event*/
		else if preg_outcome_clean in ('SAB' 'UAB' 'SB' 'MLS') and dt_gapreg - dt_lmp < 140 then preg_outcome_loss1220 = 'OUT';
		else if dt_gapreg - dt_lmp >= 140 then preg_outcome_loss1220 = 'NO';
		else preg_outcome_loss1220 = preg_outcome_clean;

	*Event time variable;
	if preg_outcome_loss1220 = 'NO' then dt_gapreg_loss1220 = dt_lmp + 139;
		else dt_gapreg_loss1220 = dt_gapreg;

	*******LOSS AFTER 20 WEEKS;

	*Event indicator;
	if preg_outcome_clean in ('SAB' 'UAB' 'SB' 'MLS') and dt_gapreg - dt_lmp < 140 then preg_outcome_loss20 = 'MSC';
		else if preg_outcome_clean in ('SAB' 'UAB' 'SB' 'MLS') then preg_outcome_loss20 = 'OUT';
		else if preg_outcome_clean in ('LBS' 'UDL' 'LBM') then preg_outcome_loss20 = 'LB';
		else preg_outcome_loss20 = preg_outcome_clean;

	*Event time variable is just dt_gapreg;

run;


/**Check;*/
/*proc freq data=loss;*/
/*	table preg_outcome_loss12 preg_outcome_loss1220 preg_outcome_loss20 / missing;*/
/*run;*/

*First, describe the cohort with all of the variables that we had proposed including;
%table1(inds = loss, colVar = exposure,
	rowVars = ga_index_days ga_index_lt14 age_at_index year_index year_le2019 
		preg_outcome_clean chronichypertension_pre preeclampsia_pre
		substance_use_pre nausea_pre any_diabetes_post pregestation_diab_post diabetes_type_post 
		t2dmrx_post t1t2dmrx_post metforrx_post
		obesity_post glp1wgtrx_post otherwgtrx_post bariatric_post
		migraine_pre recurlos_pre ckd_post 
		thyroid_disorder_post thyroidrx_post thyroid_disorder_trt_post
		depressi_post anxiety_post antideprx_post adhd_post adhdrx_post
		bipolar_post moodstabrx_post bipolar_trt_post
		ptsd_post schizo_post antipsyrx_post schizo_trt_post  
		rural hyperlip_post teratrx_pre benzorx_post anticonvulrx_post
		num_OutptPNC numOutptPNC_lt3 num_InptAdm  numInptADM_gt1,
	outfile = Table 1: Secondary loss cohort overall without weights,
	title = Table 1: Secondary loss cohort);
	

*Fit PS model within each strata. ignore variables for inverse probability of censoring weights since this is just
	to derive the PSs and SMR weights;
%competing2risk_weights(
	boot=0, 
	inds=loss, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_loss12, /*does not matter for exposure model*/
	eventDT=dt_GApreg_loss12, 
	event = 'OUT',
	cr1='IAB', 
	cr2='NO',
	psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars= ,
	dovarsmodel =  ,
	psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre  
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre, 
	doclassvars= ,
	trtvar=trt, 
	numiterations=1,
	initialseed=23244, 
	outds=test,
	outds_dist= test,
	outds_ps = ps_secondary_point,
	outds_surv = NA
	);
	
*********Look at the PS distribution by GA strata;

*Less than 14 weeks;
proc sort data=ana.ps_secondary_point_num_1; by trt; run;
proc sgplot data=ana.ps_secondary_point_num_1;
	title "Propensity scores by treatment: Less than 14 weeks at index";
	histogram d / group=trt transparency=0.2;
/* 	density d / group = trt; */
run;

proc means data=ana.ps_secondary_point_num_1 min mean max;
	title "Gestational age <14 weeks";
	by trt;
	var d;
run;

*Attach these weights to peoples observations for a weighted Table 1 by strata;
proc sql;
	create table weighted as 
	select a.*, 1 as ga_strat, b.expwgt as smrw
	from ana.primary_cohort as a
	left join ana.ps_secondary_point_num_1 as b
	on a.idxpren=b.idxpren
	/*Subset to those that remained in the trimmed sample*/
	having idxpren in (select distinct idxpren 
						from ana.ps_secondary_point_num_1 )
	;
	quit;
	
proc sort data=weighted;
	by ga_strat trt;
run;
proc means data=weighted min p25 median mean p75 max nmiss;
	by ga_strat trt;
	var smrw;
run;


*Output updated Table 1;
%table1(inds = weighted, colVar = exposure,
	rowVars = ga_index_days ga_index_lt14 age_at_index year_index year_le2019 
		preg_outcome_clean chronichypertension_pre preeclampsia_pre
		substance_use_pre nausea_pre any_diabetes_post pregestation_diab_post diabetes_type_post 
		t2dmrx_post t1t2dmrx_post metforrx_post
		obesity_post glp1wgtrx_post otherwgtrx_post bariatric_post
		migraine_pre recurlos_pre ckd_post 
		thyroid_disorder_post thyroidrx_post thyroid_disorder_trt_post
		depressi_post anxiety_post antideprx_post adhd_post adhdrx_post
		bipolar_post moodstabrx_post bipolar_trt_post
		ptsd_post schizo_post antipsyrx_post schizo_trt_post  
		rural hyperlip_post teratrx_pre benzorx_post anticonvulrx_post
		num_OutptPNC numOutptPNC_lt3 num_InptAdm  numInptADM_gt1,
	outfile = Table 1: Secondary loss cohort overall with weights,
	title = Table 1: Secondary loss cohort with weights);

	
	
/********************************************************************************************************************************************

														03 - CONDUCT PREGNANCY LOSS ANALYSIS

********************************************************************************************************************************************/

*********PREGNANCY LOSS PRIOR TO 12 WEEKS GESTATION;

*Conduct analysis for point estimate;
%competing2risk_weights(
	boot=0, 
	inds=loss, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_loss12,
	eventDT=dt_GApreg_loss12, 
	event = 'OUT',
	cr1='IAB', 
	cr2='NO',
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
	outds=ana.secondary_point_12,
	outds_dist=ga_dist_secondary_point_12,
	outds_ps = ps_secondary_point_12,
	outds_surv = NA
	);

*Get the sample sizes across the two treatment groups;
proc freq data=ana.ps_secondary_point_12_num_1;
	table trt / missing;
run;

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_secondary_point_12,
	numgastrat = 1,
	est = ana.secondary_point_12_,
	outds = ana.secondary_point_12_overall
	);
	
*Conduct bootstrapping to get the variance;
%competing2risk_weights(
	boot=1, 
	inds=loss, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_loss12,
	eventDT=dt_GApreg_loss12, 
	event = 'OUT',
	cr1='IAB', 
	cr2='NO',
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
	outds=ana.secondary_boot_12,
	outds_dist=ga_dist_secondary_boot_12,
	outds_ps = ps_secondary_boot_12,
	outds_surv = NA
	);
	
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.secondary_boot_12, 
			inputDist= ana.ga_dist_secondary_boot_12,
			numStrata= 1, 
			output_stratified= ana.secondary_boot_strat_12,
			output_overall= ana.secondary_boot_overall_12);

*Overall estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.secondary_boot_12_1, pointdsn=ana.secondary_point_12_1, output=ana.secondary_boot_12_1_ci);


%count_missing_zero(inds1=ana.secondary_boot_12_1, inds2=NA);
	

	
	


*********PREGNANCY LOSS PRIOR >= 12 WEEKS OR <20 WEEKS OF GESTATION;

*Conduct analysis for point estimate;
%competing3risk_weights(
	boot=0, 
	inds=loss, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_loss1220,
	eventDT=dt_GApreg_loss1220, 
	event = 'OUT',
	cr1='IAB', 
	cr2='EAR',
	cr3='NO',
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
	outds=ana.secondary_point_120,
	outds_dist=ga_dist_secondary_point_120,
	outds_ps = ps_secondary_point_120,
	outds_surv = NA
	);

*Get the sample sizes across the two treatment groups;
proc freq data=ana.ps_secondary_point_120_num_1;
	table trt / missing;
run;

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_secondary_point_120,
	numgastrat = 1,
	est = ana.secondary_point_120_,
	outds = ana.secondary_point_120_overall
	);
	
*Conduct bootstrapping to get the variance;
%competing3risk_weights(
	boot=1, 
	inds=loss, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_loss1220,
	eventDT=dt_GApreg_loss1220, 
	event = 'OUT',
	cr1='IAB', 
	cr2='EAR',
	cr3='NO',
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
	outds=ana.secondary_boot_120,
	outds_dist=ga_dist_secondary_boot_120,
	outds_ps = ps_secondary_boot_120,
	outds_surv = NA
	);
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.secondary_boot_120, 
			inputDist= ana.ga_dist_secondary_boot_120,
			numStrata= 1, 
			output_stratified= ana.secondary_boot_strat_120,
			output_overall= ana.secondary_boot_overall_120);

*Overall estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.secondary_boot_120_1, pointdsn=ana.secondary_point_120_1, output=ana.secondary_boot_120_1_ci);


%count_missing_zero(inds1=ana.secondary_boot_120_1, inds2=NA);







*********PREGNANCY LOSS >= 20 WEEKS OF GESTATION;

*Conduct analysis for point estimate;
%competing3risk_weights(
	boot=0, 
	inds=loss, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_loss20,
	eventDT=dt_GApreg, 
	event = 'OUT',
	cr1='IAB', 
	cr2='MSC',
	cr3='LB',
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
	outds=ana.secondary_point_20,
	outds_dist=ga_dist_secondary_point_20,
	outds_ps = ps_secondary_point_20,
	outds_surv = NA
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_secondary_point_20,
	numgastrat = 1,
	est = ana.secondary_point_20_,
	outds = ana.secondary_point_20_overall
	);
	
*Conduct bootstrapping to get the variance;
%competing3risk_weights(
	boot=1, 
	inds=loss, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_loss20,
	eventDT=dt_GApreg, 
	event = 'OUT',
	cr1='IAB', 
	cr2='MSC',
	cr3='LBM' 'LBS' 'UDL',
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
	outds=ana.secondary_boot_20,
	outds_dist=ga_dist_secondary_boot_20,
	outds_ps = ps_secondary_boot_20,
	outds_surv = NA
	);
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.secondary_boot_20, 
			inputDist= ana.ga_dist_secondary_boot_20,
			numStrata= 1, 
			output_stratified= ana.secondary_boot_strat_20,
			output_overall= ana.secondary_boot_overall_20);

*Overall estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.secondary_boot_20_1, pointdsn=ana.secondary_point_20_1, output=ana.secondary_boot_20_1_ci);

%count_missing_zero(inds1=ana.secondary_boot_20_1, inds2=NA);
