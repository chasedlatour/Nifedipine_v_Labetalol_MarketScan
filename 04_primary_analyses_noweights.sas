/********************************************************************************************************************************************
PROGRAM: 04_primary_analyses_noweights.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to derive the primary study estimates without SMR weights and with only SMR weights.
	
Goal: 
Output data: 

Date: 12.19.2024
********************************************************************************************************************************************/






/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - CONDUCT PREGNANCY LOSS ANALYSIS WITHOUT WEIGHTS
	- 02 - CONDUCT PREGNANCY LOSS ANALYSIS WITH SMR WEIGHTS ONLY
	- 03 - CONDUCT PRETERM BIRTH ANALYSIS WITHOUT WEIGHTS
	- 04 - CONDUCT PRETERM BIRTH ANALYSIS WITH SMR WEIGHTS ONLY

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
%setup(sample=full, programname=04_primary_analyses_noweights, savelog=Y);

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

													01 - CONDUCT PREGNANCY LOSS ANALYSIS WITHOUT WEIGHTS

********************************************************************************************************************************************/

*Get the point estimate;
%competing2risk_weights(
	boot=0, 
	inds=ana.primary_cohort, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_clean,
	eventDT=dt_GApreg, 
	event = 'SAB' 'UAB' 'SB' 'MLS',
	cr1='IAB', 
	cr2='LBM' 'LBS' 'UDL',
	psvars=,
	dovars= ,
	dovarsmodel= ,
	psclassvars=, 
	doclassvars= ,
	trtvar=trt, 
	numiterations=1,
	initialseed=23244, 
	outds=ana.primary_point_nowgt,
	outds_dist=ga_dist_primary_point_nowgt,
	outds_ps = ps_primary_point_nowgt,
	outds_surv = NA,
	smr=0,
	ipcw=0
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_primary_point_nowgt,
	numgastrat = 2,
	est = ana.primary_point_nowgt_,
	outds = ana.primary_point_overall_nowgt
	);

*Conduct the boostrap so that we can get standard errors;
%competing2risk_weights(
	boot=1, 
	inds=ana.primary_cohort, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_clean,
	eventDT=dt_GApreg, 
	event = 'SAB' 'UAB' 'SB' 'MLS',
	cr1='IAB', 
	cr2='LBM' 'LBS' 'UDL',
	psvars=,
	dovars=,
	dovarsmodel = ,
	psclassvars=, 
	doclassvars=,
	trtvar=trt, 
	numiterations=1000,
	initialseed=23244, 
	outds=ana.primary_boot_nowgt,
	outds_dist=gadist_primary_boot_nowgt,
	outds_ps = NA,
	outds_surv = NA,
	smr=0,
	ipcw=0
	);
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.primary_boot_nowgt, 
			inputDist= ana.gadist_primary_boot_nowgt,
			numStrata= 2, 
			output_stratified= ana.primary_boot_strat_nowgt,
			output_overall= ana.primary_boot_overall_nowgt);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.primary_boot_nowgt_1, pointdsn=ana.primary_point_nowgt_1, output=ana.primary_boot_nowgt_1_ci);
%strat_estimates_w_CI(bootdsn=ana.primary_boot_nowgt_2, pointdsn=ana.primary_point_nowgt_2, output=ana.primary_boot_nowgt_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.primary_boot_overall_nowgt, pointdsn=ana.primary_point_overall_nowgt, output=ana.primary_overall_boot_nowgt_ci);
	
%count_missing_zero(inds1=ana.primary_boot_nowgt_1, inds2=ana.primary_boot_nowgt_2);










/********************************************************************************************************************************************

											02 - CONDUCT PREGNANCY LOSS ANALYSIS WITH SMR WEIGHTS ONLY

********************************************************************************************************************************************/

*Get the point estimate;
%competing2risk_weights(
	boot=0, 
	inds=ana.primary_cohort, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_clean,
	eventDT=dt_GApreg, 
	event = 'SAB' 'UAB' 'SB' 'MLS',
	cr1='IAB', 
	cr2='LBM' 'LBS' 'UDL',
	psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars= ,
	dovarsmodel= ,
	psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre, 
	doclassvars= ,
	trtvar=trt, 
	numiterations=1,
	initialseed=23244, 
	outds=ana.primary_point_noipcw,
	outds_dist=ga_dist_primary_point_noipcw,
	outds_ps = ps_primary_point_noipcw,
	outds_surv = NA,
	smr=1,
	ipcw=0
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_primary_point_noipcw,
	numgastrat = 2,
	est = ana.primary_point_noipcw_,
	outds = ana.primary_point_overall_noipcw
	);

*Conduct the boostrap so that we can get standard errors;
%competing2risk_weights(
	boot=1, 
	inds=ana.primary_cohort, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_clean,
	eventDT=dt_GApreg, 
	event = 'SAB' 'UAB' 'SB' 'MLS',
	cr1='IAB', 
	cr2='LBM' 'LBS' 'UDL',
	psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars= ,
	dovarsmodel= ,
	psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre,  
	doclassvars=,
	trtvar=trt, 
	numiterations=1000,
	initialseed=23244, 
	outds=ana.primary_boot_noipcw,
	outds_dist=gadist_primary_boot_noipcw,
	outds_ps = NA,
	outds_surv = NA,
	smr=1,
	ipcw=0
	);
	
options mlogic mprint symbolgen notes; 
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.primary_boot_noipcw, 
			inputDist= ana.gadist_primary_boot_noipcw,
			numStrata= 2, 
			output_stratified= ana.primary_boot_strat_noipcw,
			output_overall= ana.primary_boot_overall_noipcw);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.primary_boot_noipcw_1, pointdsn=ana.primary_point_noipcw_1, output=ana.primary_boot_noipcw_1_ci);
%strat_estimates_w_CI(bootdsn=ana.primary_boot_noipcw_2, pointdsn=ana.primary_point_noipcw_2, output=ana.primary_boot_noipcw_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.primary_boot_overall_noipcw, pointdsn=ana.primary_point_overall_noipcw, output=ana.primary_overall_boot_noipcw_ci);



%count_missing_zero(inds1=ana.primary_boot_noipcw_1, inds2=ana.primary_boot_noipcw_2);









/********************************************************************************************************************************************

													04 - CONDUCT PRETERM BIRTH ANALYSIS WITHOUT WEIGHTS

********************************************************************************************************************************************/

*Get the point estimate;
%competing3risk_weights(
	boot=0, 
	inds=ana.primary_cohort, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_ptb,
	eventDT=dt_GApreg_ptb, 
	event = 'PTB',
	cr1='IAB', 
	cr2='SAB' 'UAB' 'SB' 'MLS',
	cr3='TB',
	psvars=,
	dovars= ,
	dovarsmodel= ,
	psclassvars=, 
	doclassvars= ,
	trtvar=trt, 
	numiterations=1,
	initialseed=23244, 
	outds=ana.primary_point_ptb_nowgt,
	outds_dist=ga_dist_primary_point_ptb_nowgt,
	outds_ps = NA,
	outds_surv = NA,
	smr=0,
	ipcw=0
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_primary_point_ptb_nowgt,
	numgastrat = 2,
	est = ana.primary_point_ptb_nowgt_,
	outds = ana.primary_point_ptb_overall_nowgt
	);

*Conduct the boostrap so that we can get standard errors;
%competing3risk_weights(
	boot=1, 
	inds=ana.primary_cohort, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_ptb,
	eventDT=dt_GApreg_ptb, 
	event = 'PTB',
	cr1='IAB', 
	cr2='SAB' 'UAB' 'SB' 'MLS',
	cr3='TB',
	psvars=,
	dovars=,
	dovarsmodel = ,
	psclassvars=, 
	doclassvars=,
	trtvar=trt, 
	numiterations=1000,
	initialseed=23244, 
	outds=ana.primary_boot_ptb_nowgt,
	outds_dist=gadist_primary_boot_ptb_nowgt,
	outds_ps = NA,
	outds_surv = NA,
	smr=0,
	ipcw=0
	);
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.primary_boot_ptb_nowgt, 
			inputDist= ana.gadist_primary_boot_ptb_nowgt,
			numStrata= 2, 
			output_stratified= ana.primary_boot_ptb_strat_nowgt,
			output_overall= ana.primary_boot_ptb_overall_nowgt);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.primary_boot_ptb_nowgt_1, pointdsn=ana.primary_point_ptb_nowgt_1, output=ana.primary_boot_ptb_nowgt_1_ci);
%strat_estimates_w_CI(bootdsn=ana.primary_boot_ptb_nowgt_2, pointdsn=ana.primary_point_ptb_nowgt_2, output=ana.primary_boot_ptb_nowgt_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.primary_boot_ptb_overall_nowgt, pointdsn=ana.primary_point_ptb_overall_nowgt, output=ana.primary_all_boot_ptb_nowgt_ci);
	
%count_missing_zero(inds1=ana.primary_boot_ptb_nowgt_1, inds2=ana.primary_boot_ptb_nowgt_2);










/********************************************************************************************************************************************

											02 - CONDUCT PREGNANCY LOSS ANALYSIS WITH SMR WEIGHTS ONLY

********************************************************************************************************************************************/

*Get the point estimate;
%competing3risk_weights(
	boot=0, 
	inds=ana.primary_cohort, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_ptb,
	eventDT=dt_GApreg_ptb, 
	event = 'PTB',
	cr1='IAB', 
	cr2='SAB' 'UAB' 'SB' 'MLS',
	cr3='TB',
	psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars= ,
	dovarsmodel= ,
	psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre, 
	doclassvars= ,
	trtvar=trt, 
	numiterations=1,
	initialseed=23244, 
	outds=ana.primary_point_ptb_noipcw,
	outds_dist=gadist_primary_point_ptb_noipcw,
	outds_ps = NA,
	outds_surv = NA,
	smr=1,
	ipcw=0
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group; 
%combine_point_estimates(
	gadist = ana.gadist_primary_point_ptb_noipcw,
	numgastrat = 2,
	est = ana.primary_point_ptb_noipcw_,
	outds = ana.primary_point_overall_ptb_noipcw
	);

*Conduct the boostrap so that we can get standard errors; *No errors;
%competing3risk_weights(
	boot=1, 
	inds=ana.primary_cohort, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_ptb,
	eventDT=dt_GApreg_ptb, 
	event = 'PTB',
	cr1='IAB', 
	cr2='SAB' 'UAB' 'SB' 'MLS',
	cr3='TB',
	psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars= ,
	dovarsmodel= ,
	psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre,  
	doclassvars=,
	trtvar=trt, 
	numiterations=1000,
	initialseed=23244, 
	outds=ana.primary_boot_ptb_noipcw,
	outds_dist=gadist_primary_boot_ptb_noipcw,
	outds_ps = NA,
	outds_surv = NA,
	smr=1,
	ipcw=0
	);
	
options mlogic mprint symbolgen notes; 
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.primary_boot_ptb_noipcw, 
			inputDist= ana.gadist_primary_boot_ptb_noipcw,
			numStrata= 2, 
			output_stratified= ana.primary_boot_ptb_strat_noipcw,
			output_overall= ana.primary_boot_overall_ptb_noipcw);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.primary_boot_ptb_noipcw_1, pointdsn=ana.primary_point_ptb_noipcw_1, output=ana.primary_boot_ptb_noipcw_1_ci);
%strat_estimates_w_CI(bootdsn=ana.primary_boot_ptb_noipcw_2, pointdsn=ana.primary_point_ptb_noipcw_2, output=ana.primary_boot_ptb_noipcw_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.primary_boot_overall_ptb_noipcw, pointdsn=ana.primary_point_overall_ptb_noipcw, output=ana.primary_all_boot_ptb_noipcw_ci);

%count_missing_zero(inds1=ana.primary_boot_ptb_noipcw_1, inds2=ana.primary_boot_ptb_noipcw_2);

	


	
	
