/********************************************************************************************************************************************
PROGRAM: 04_posthoc_analyses_preec_noweights.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to derive the primary study estimates without SMR weights and with only SMR weights.
	
Goal: 
Output data: 

Date: 12.19.2024
********************************************************************************************************************************************/






/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - ANY PREECLAMPSIA ANALYSIS WITHOUT WEIGHTS
	- 02 - ANY PREECLAMPSIA ANALYSIS WITH SMR WEIGHTS ONLY
	- 03 - INPATIENT PREECLAMPSIA ANALYSIS WITHOUT WEIGHTS
	- 04 - INPATIENT PREECLAMPSIA ANALYSIS WITH SMR WEIGHTS ONLY

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
%setup(sample=full, programname=04_posthoc_analyses_preec_noweights, savelog=Y);

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

													01 - ANY PREECLAMPSIA WITHOUT WEIGHTS

********************************************************************************************************************************************/

*Get the point estimate;
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
	cr3='LBS' 'LBM' 'UDL',
	psvars=,
	dovars= ,
	dovarsmodel= ,
	psclassvars=, 
	doclassvars= ,
	trtvar=trt, 
	numiterations=1,
	initialseed=23244, 
	outds=ana.prim_point_pany_nowgt,
	outds_dist=ga_prim_point_pany_nowgt,
	outds_ps = NA,
	outds_surv = NA,
	smr=0,
	ipcw=0
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_prim_point_pany_nowgt,
	numgastrat = 2,
	est = ana.prim_point_pany_nowgt_,
	outds = ana.prim_point_pany_overall_nowgt
	);

*Conduct the boostrap so that we can get standard errors;
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
	cr3='LBS' 'LBM' 'UDL',
	psvars=,
	dovars=,
	dovarsmodel = ,
	psclassvars=, 
	doclassvars=,
	trtvar=trt, 
	numiterations=1000,
	initialseed=23244, 
	outds=ana.prim_boot_pany_nowgt,
	outds_dist=ga_prim_boot_pany_nowgt,
	outds_ps = NA,
	outds_surv = NA,
	smr=0,
	ipcw=0
	);
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.prim_boot_pany_nowgt, 
			inputDist= ana.ga_prim_boot_pany_nowgt,
			numStrata= 2, 
			output_stratified= ana.prim_boot_pany_strat_nowgt,
			output_overall= ana.prim_boot_pany_overall_nowgt);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.prim_boot_pany_nowgt_1, pointdsn=ana.prim_point_pany_nowgt_1, output=ana.prim_boot_pany_nowgt_1_ci);
%strat_estimates_w_CI(bootdsn=ana.prim_boot_pany_nowgt_2, pointdsn=ana.prim_point_pany_nowgt_2, output=ana.prim_boot_pany_nowgt_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.prim_boot_pany_overall_nowgt, pointdsn=ana.prim_point_pany_overall_nowgt, output=ana.prim_all_boot_pany_nowgt_ci);
	
%count_missing_zero(inds1=ana.prim_boot_pany_nowgt_1, inds2=ana.prim_boot_pany_nowgt_2);










/********************************************************************************************************************************************

											02 - ANY PREECLAMPSIA WITH SMR WEIGHTS ONLY

********************************************************************************************************************************************/

*Get the point estimate;
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
	cr3='LBS' 'LBM' 'UDL',
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
	outds=ana.prim_point_pany_noipcw,
	outds_dist=ga_prim_point_pany_noipcw,
	outds_ps = NA,
	outds_surv = NA,
	smr=1,
	ipcw=0
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group; 
%combine_point_estimates(
	gadist = ana.ga_prim_point_pany_noipcw,
	numgastrat = 2,
	est = ana.prim_point_pany_noipcw_,
	outds = ana.prim_point_overall_pany_noipcw
	);

*Conduct the boostrap so that we can get standard errors; *No errors;
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
	cr3='LBS' 'LBM' 'UDL',
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
	outds=ana.prim_boot_pany_noipcw,
	outds_dist=ga_prim_boot_pany_noipcw,
	outds_ps = NA,
	outds_surv = NA,
	smr=1,
	ipcw=0
	);
	
options mlogic mprint symbolgen notes; 
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.prim_boot_pany_noipcw, 
			inputDist= ana.ga_prim_boot_pany_noipcw,
			numStrata= 2, 
			output_stratified= ana.prim_boot_pany_strat_noipcw,
			output_overall= ana.prim_boot_overall_pany_noipcw);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.prim_boot_pany_noipcw_1, pointdsn=ana.prim_point_pany_noipcw_1, output=ana.prim_boot_pany_noipcw_1_ci);
%strat_estimates_w_CI(bootdsn=ana.prim_boot_pany_noipcw_2, pointdsn=ana.prim_point_pany_noipcw_2, output=ana.prim_boot_pany_noipcw_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.prim_boot_overall_pany_noipcw, pointdsn=ana.prim_point_overall_pany_noipcw, output=ana.prim_all_boot_pany_noipcw_ci);

%count_missing_zero(inds1=ana.prim_boot_pany_noipcw_1, inds2=ana.prim_boot_pany_noipcw_2);

	



















/********************************************************************************************************************************************

												03 - INPATIENT PREECLAMPSIA WITHOUT WEIGHTS

********************************************************************************************************************************************/

*Get the point estimate;
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
	cr3='LBS' 'LBM' 'UDL',
	psvars=,
	dovars= ,
	dovarsmodel= ,
	psclassvars=, 
	doclassvars= ,
	trtvar=trt, 
	numiterations=1,
	initialseed=23244, 
	outds=ana.prim_point_pinpt_nowgt,
	outds_dist=ga_prim_point_pinpt_nowgt,
	outds_ps = NA,
	outds_surv = NA,
	smr=0,
	ipcw=0
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_prim_point_pinpt_nowgt,
	numgastrat = 2,
	est = ana.prim_point_pinpt_nowgt_,
	outds = ana.prim_point_pinpt_overall_nowgt
	);

*Conduct the boostrap so that we can get standard errors;
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
	cr3='LBS' 'LBM' 'UDL',
	psvars=,
	dovars=,
	dovarsmodel = ,
	psclassvars=, 
	doclassvars=,
	trtvar=trt, 
	numiterations=1000,
	initialseed=23244, 
	outds=ana.prim_boot_pinpt_nowgt,
	outds_dist=ga_prim_boot_pinpt_nowgt,
	outds_ps = NA,
	outds_surv = NA,
	smr=0,
	ipcw=0
	);
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.prim_boot_pinpt_nowgt, 
			inputDist= ana.ga_prim_boot_pinpt_nowgt,
			numStrata= 2, 
			output_stratified= ana.prim_boot_pinpt_strat_nowgt,
			output_overall= ana.prim_boot_pinpt_overall_nowgt);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.prim_boot_pinpt_nowgt_1, pointdsn=ana.prim_point_pinpt_nowgt_1, output=ana.prim_boot_pinpt_nowgt_1_ci);
%strat_estimates_w_CI(bootdsn=ana.prim_boot_pinpt_nowgt_2, pointdsn=ana.prim_point_pinpt_nowgt_2, output=ana.prim_boot_pinpt_nowgt_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.prim_boot_pinpt_overall_nowgt, pointdsn=ana.prim_point_pinpt_overall_nowgt, output=ana.prim_all_boot_pinpt_nowgt_ci);
	
%count_missing_zero(inds1=ana.prim_boot_pinpt_nowgt_1, inds2=ana.prim_boot_pinpt_nowgt_2);










/********************************************************************************************************************************************

											04 - INPATIENT PREECLAMPSIA WITH SMR WEIGHTS ONLY

********************************************************************************************************************************************/

*Get the point estimate;
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
	cr3='LBS' 'LBM' 'UDL',
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
	outds=ana.prim_point_pinpt_noipcw,
	outds_dist=ga_prim_point_pinpt_noipcw,
	outds_ps = NA,
	outds_surv = NA,
	smr=1,
	ipcw=0
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group; 
%combine_point_estimates(
	gadist = ana.ga_prim_point_pinpt_noipcw,
	numgastrat = 2,
	est = ana.prim_point_pinpt_noipcw_,
	outds = ana.prim_point_overall_pinpt_noipcw
	);

*Conduct the boostrap so that we can get standard errors; *No errors;
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
	cr3='LBS' 'LBM' 'UDL',
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
	outds=ana.prim_boot_pinpt_noipcw,
	outds_dist=ga_prim_boot_pinpt_noipcw,
	outds_ps = NA,
	outds_surv = NA,
	smr=1,
	ipcw=0
	);
	
options mlogic mprint symbolgen notes; 
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.prim_boot_pinpt_noipcw, 
			inputDist= ana.ga_prim_boot_pinpt_noipcw,
			numStrata= 2, 
			output_stratified= ana.prim_boot_pinpt_strat_noipcw,
			output_overall= ana.prim_boot_overall_pinpt_noipcw);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.prim_boot_pinpt_noipcw_1, pointdsn=ana.prim_point_pinpt_noipcw_1, output=ana.prim_boot_pinpt_noipcw_1_ci);
%strat_estimates_w_CI(bootdsn=ana.prim_boot_pinpt_noipcw_2, pointdsn=ana.prim_point_pinpt_noipcw_2, output=ana.prim_boot_pinpt_noipcw_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.prim_boot_overall_pinpt_noipcw, pointdsn=ana.prim_point_overall_pinpt_noipcw, output=ana.prim_all_boot_pinpt_noipcw_ci);

%count_missing_zero(inds1=ana.prim_boot_pinpt_noipcw_1, inds2=ana.prim_boot_pinpt_noipcw_2);

	

	
	
