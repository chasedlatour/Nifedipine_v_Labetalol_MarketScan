/********************************************************************************************************************************************
PROGRAM: 16_sens_covariates_pre.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to conduct sensitivity analyses where we only control for covariate information that was 
captured prior to the index date.

NOTE: Using only pre-index information resulted in much fewer nifedipine initiators being flagged for T2DM rx fills. As a result, 
we could not model it as flexibly. We removed the interaction between the number of fills and year of index and diabetes. We otherwise
remained faithful to the primary analysis. We confirmed that the SMDs after SMR weighting were still reasonable.
	
Goal: 
Output data: 

Date: 12.19.2024
********************************************************************************************************************************************/






/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - PRINT UNWEIGHTED AND WEIGHTED TABLE 1
	- 02 - CONDUCT PREGNANCY LOSS ANALYSIS
	- 03 - CONDUCT PRETERM BIRTH ANALYSIS
	- 04 - AD-HOC DESCRIPTIVES

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
%setup(sample=full, programname=16_sens_covariates_pre, savelog=Y);

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

													01 - PRINT UNWEIGHTED AND WEIGHTED TABLE 1

********************************************************************************************************************************************/


%*Create pre version of diabetes_simp variable;
data ana.primary_cohort;
set ana.primary_cohort;

	/*Need a pre-version of this variable*/
	if diabetes_type_pre in ("T1DM" "T2DM") then diabetes_simp_pre = "PRE";
		else if diabetes_type_pre = "NA" then diabetes_simp_pre = "NA";
		else diabetes_simp_pre = "OTH";

run;




%*Stratify by gestational age at index date;
data primary_lt14 primary_ge14;
set ana.primary_cohort;
	if ga_index_lt14 = 1 then output primary_lt14;
		else output primary_ge14;
run;

*Now output a table for each strata;

%*Now do some table 1 descriptives of the cohort - ana.primary_cohort;
%table1(inds = primary_lt14, colVar = exposure,
	rowVars = ga_index_days age_at_index age_at_index_cat year_index year_le2019 chronichypertension_pre preeclampsia_pre
		substance_use_pre nausea_pre pregestation_diab_pre t2dmrx_pre t1t2dmrx_pre metforrx_pre
		obesity_pre glp1wgtrx_pre otherwgtrx_pre bariatric_pre migraine_pre recurlos_pre ckd_pre 
		thyroid_disorder_pre thyroidrx_pre
		depressi_pre anxiety_pre antideprx_pre adhd_pre adhdrx_pre bipolar_pre moodstabrx_pre
		ptsd_pre schizo_pre antipsyrx_pre rural hyperlip_pre teratrx_pre benzorx_pre anticonvulrx_pre
		num_OutptPNC numOutptPNC_lt3 num_InptAdm  numInptADM_gt1,
	outfile = Table 1: Unweighted primary cohort GA lt 14wk pre index covariates,
	title = Table 1: Primary cohort where GA less than 14w at index);
	
%*Now do some table 1 descriptives of the cohort - ana.primary_cohort;
%table1(inds = primary_ge14, colVar = exposure,
	rowVars = ga_index_days age_at_index age_at_index_cat year_index year_le2019 chronichypertension_pre preeclampsia_pre
		substance_use_pre nausea_pre pregestation_diab_pre t2dmrx_pre t1t2dmrx_pre metforrx_pre
		obesity_pre glp1wgtrx_pre otherwgtrx_pre bariatric_pre migraine_pre recurlos_pre ckd_pre 
		thyroid_disorder_pre thyroidrx_pre
		depressi_pre anxiety_pre antideprx_pre adhd_pre adhdrx_pre bipolar_pre moodstabrx_pre
		ptsd_pre schizo_pre antipsyrx_pre rural hyperlip_pre teratrx_pre benzorx_pre anticonvulrx_pre
		num_OutptPNC numOutptPNC_lt3 num_InptAdm  numInptADM_gt1,
	outfile = Table 1: Unweighted primary cohort GA ge 14wk pre index covariates,
	title = Table 1: Primary cohort where GA less than 14w at index);
	
	
	
*Fit PS model within each strata. Include variables for inverse probability of censoring weights;

%competing2risk_weights(
	boot=0, 
	inds=ana.primary_cohort, 
	gacatvar = ga_index_cat, /*ga_index_cat_pre,*/
	startDT=dt_index, 
	outcomevar = preg_outcome_clean,
	eventDT=dt_GApreg, 
	event = 'SAB' 'UAB' 'SB' 'MLS',
	cr1='IAB', 
	cr2='LBM' 'LBS' 'UDL',
	psvars=/*ga_quartile_pre*/ ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
			/*year_index2017*t2dmrx_pre diabetes_simp_pre*t2dmrx_pre */
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre
			depressi_pre anxiety_pre antideprx_pre benzorx_pre
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars= ,
	dovarsmodel =  ,
	psclassvars=/*ga_quartile_pre*/ ga_quartile year_index2017 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre  
			depressi_pre anxiety_pre antideprx_pre benzorx_pre teratrx_pre, 
	doclassvars= ,
	trtvar=trt, 
	numiterations=1,
	initialseed=23244, 
	outds=test, /*ana.primary_point,*/
	outds_dist=NA, /*ga_dist_primary_point,*/
	outds_ps = ps_sens_precov,
	outds_surv = NA /*ana.surv_primary_point*/
	);


*Attach these weights to peoples observations for a weighted Table 1 by strata;
proc sql;
	create table weighted as 
	select a.*, case when missing(c.expwgt) then 1 else 2 end as ga_strat,
			coalesce(b.expwgt, c.expwgt) as smrw
	from ana.primary_cohort as a
/*	left join ana.ps_sens_precov_num_1 as b*/
/*	on a.idxpren=b.idxpren*/
	left join ana.ps_sens_precov_num_1 as b
	on a.idxpren=b.idxpren
	left join ana.ps_sens_precov_num_2 as c
	on a.idxpren=c.idxpren
	/*Subset to those that remained in the trimmed sample*/
/*	having idxpren in (select distinct idxpren*/
/*						from ana.ps_sens_precov_num_1)*/
	having idxpren in (select distinct idxpren 
						from ana.ps_sens_precov_num_1 
						union
						select distinct idxpren
						from ana.ps_sens_precov_num_2)
	;
	quit;
	
proc sort data=weighted;
	by ga_strat trt;
run;
proc means data=weighted min p25 median mean p75 max nmiss;
	class ga_strat trt;
	var smrw;
run;


%*Stratify by gestational age at index date;
data primary_lt14 primary_ge14;
set weighted;
	if ga_strat = 1 then output primary_lt14;
		else output primary_ge14;
run;

*Now output a table for each strata;

%*Now do some table 1 descriptives of the cohort - ana.primary_cohort;
%table1(inds = primary_lt14, colVar = exposure, wgtvar=smrw,
	rowVars = ga_index_days age_at_index age_at_index_cat year_index year_le2019 chronichypertension_pre preeclampsia_pre
		substance_use_pre nausea_pre pregestation_diab_pre t2dmrx_pre t1t2dmrx_pre metforrx_pre
		obesity_pre glp1wgtrx_pre otherwgtrx_pre bariatric_pre migraine_pre recurlos_pre ckd_pre 
		thyroid_disorder_pre thyroidrx_pre
		depressi_pre anxiety_pre antideprx_pre adhd_pre adhdrx_pre bipolar_pre moodstabrx_pre
		ptsd_pre schizo_pre antipsyrx_pre rural hyperlip_pre teratrx_pre benzorx_pre anticonvulrx_pre
		num_OutptPNC numOutptPNC_lt3 num_InptAdm  numInptADM_gt1,
	outfile = Table 1: Weighted primary cohort lt 14w preindex covariates,
	title = Table 1: Primary cohort where GA less than 14w at index);
	
%*Now do some table 1 descriptives of the cohort - ana.primary_cohort;
%table1(inds = primary_ge14, colVar = exposure, wgtvar=smrw,
	rowVars = ga_index_days age_at_index age_at_index_cat year_index year_le2019 chronichypertension_pre preeclampsia_pre
		substance_use_pre nausea_pre pregestation_diab_pre t2dmrx_pre t1t2dmrx_pre metforrx_pre
		obesity_pre glp1wgtrx_pre otherwgtrx_pre bariatric_pre migraine_pre recurlos_pre ckd_pre 
		thyroid_disorder_pre thyroidrx_pre
		depressi_pre anxiety_pre antideprx_pre adhd_pre adhdrx_pre bipolar_pre moodstabrx_pre
		ptsd_pre schizo_pre antipsyrx_pre rural hyperlip_pre teratrx_pre benzorx_pre anticonvulrx_pre
		num_OutptPNC numOutptPNC_lt3 num_InptAdm  numInptADM_gt1,
	outfile = Table 1: Weighted primary cohort ge 14w preindex covariates,
	title = Table 1: Primary cohort where GA at least 14w at index);
	
	
	
	
/********************************************************************************************************************************************

														03 - CONDUCT PREGNANCY LOSS ANALYSIS

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
			year_index2017 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
/*			year_index2017*t2dmrx_post diabetes_simp*t2dmrx_post*/
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre
			depressi_pre anxiety_pre antideprx_pre benzorx_pre
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars=ga_quartile age_at_index year_index4 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre depressi_pre anxiety_pre antideprx_pre 
			benzorx_pre teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	dovarsmodel = ga_quartile age_at_index year_index4
			t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre /*year_index4*t2dmrx_pre*/
			/*diabetes_simp_pre*t2dmrx_pre*/
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre 
			depressi_pre anxiety_pre antideprx_pre benzorx_pre
			teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	psclassvars=ga_quartile year_index2017 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre 
			depressi_pre anxiety_pre antideprx_pre benzorx_pre teratrx_pre, 
	doclassvars=ga_quartile year_index4 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre
			depressi_pre anxiety_pre antideprx_pre benzorx_pre teratrx_pre rural2,
	trtvar=trt,
	numiterations=1,
	initialseed=23244, 
	outds=ana.senspre_point,
	outds_dist=ga_dist_senspre_point,
	outds_ps = ps_senspre_point,
	outds_surv = NA
	);


*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_senspre_point,
	numgastrat = 2,
	est = ana.senspre_point_,
	outds = ana.senspre_point_overall
	);

*********Finally, conduct the bootstrap;
	
*Now incorporate the IPCW;
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
			year_index2017 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
/*			year_index2017*t2dmrx_post diabetes_simp*t2dmrx_post*/
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre
			depressi_pre anxiety_pre antideprx_pre benzorx_pre
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars=ga_quartile age_at_index year_index4 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre depressi_pre anxiety_pre antideprx_pre 
			benzorx_pre teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	dovarsmodel = ga_quartile age_at_index year_index4
			t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre /*year_index4*t2dmrx_pre*/
			/*diabetes_simp_pre*t2dmrx_pre*/
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre 
			depressi_pre anxiety_pre antideprx_pre benzorx_pre
			teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	psclassvars=ga_quartile year_index2017 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre 
			depressi_pre anxiety_pre antideprx_pre benzorx_pre teratrx_pre, 
	doclassvars=ga_quartile year_index4 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre
			depressi_pre anxiety_pre antideprx_pre benzorx_pre teratrx_pre rural2,
	trtvar=trt, 
	numiterations=1000,
	initialseed=23244, 
	outds=ana.senspre_boot,
	outds_dist=ga_dist_senspre_boot,
	outds_ps = ps_senspre_boot,
	outds_surv = NA
	);
	
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.senspre_boot, 
			inputDist= ana.ga_dist_senspre_boot,
			numStrata= 2, 
			output_stratified= ana.senspre_boot_strat,
			output_overall= ana.senspre_boot_overall);


*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.senspre_boot_1, pointdsn=ana.senspre_point_1, output=ana.senspre_boot_1_ci);
%strat_estimates_w_CI(bootdsn=ana.senspre_boot_2, pointdsn=ana.senspre_point_2, output=ana.senspre_boot_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.senspre_boot_overall, pointdsn=ana.senspre_point_overall, output=ana.senspre_overall_boot_ci);
	
	
%count_missing_zero(inds1=ana.senspre_boot_1, inds2=ana.senspre_boot_2);
















	
/********************************************************************************************************************************************

														04 - CONDUCT PRETERM BIRTH ANALYSIS

********************************************************************************************************************************************/


*********Conduct analyses with full adjustment;
	
*Now incorporate the IPCW;
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
			year_index2017 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
/*			year_index2017*t2dmrx_post diabetes_simp*t2dmrx_post*/
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre
			depressi_pre anxiety_pre antideprx_pre benzorx_pre
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars=ga_quartile age_at_index year_index4 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre depressi_pre anxiety_pre antideprx_pre 
			benzorx_pre teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	dovarsmodel = ga_quartile age_at_index year_index4
			t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre /*year_index4*t2dmrx_pre*/
			/*diabetes_simp_pre*t2dmrx_pre*/
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre 
			depressi_pre anxiety_pre antideprx_pre benzorx_pre
			teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	psclassvars=ga_quartile year_index2017 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre 
			depressi_pre anxiety_pre antideprx_pre benzorx_pre teratrx_pre, 
	doclassvars=ga_quartile year_index4 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre
			depressi_pre anxiety_pre antideprx_pre benzorx_pre teratrx_pre rural2,
	trtvar=trt,
	numiterations=1,
	initialseed=23244, 
	outds=ana.senspre_point_ptb,
	outds_dist=ga_dist_senspre_point_ptb,
	outds_ps = ps_senspre_point_ptb,
	outds_surv = NA
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_senspre_point_ptb,
	numgastrat = 2,
	est = ana.senspre_point_ptb_,
	outds = ana.senspre_point_ptb_overall
	);

*********Finally, conduct the bootstrap;
	
*Now incorporate the IPCW;
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
			year_index2017 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
/*			year_index2017*t2dmrx_post diabetes_simp*t2dmrx_post*/
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre
			depressi_pre anxiety_pre antideprx_pre benzorx_pre
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars=ga_quartile age_at_index year_index4 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre depressi_pre anxiety_pre antideprx_pre 
			benzorx_pre teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	dovarsmodel = ga_quartile age_at_index year_index4
			t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre /*year_index4*t2dmrx_pre*/
			/*diabetes_simp_pre*t2dmrx_pre*/
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre 
			depressi_pre anxiety_pre antideprx_pre benzorx_pre
			teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	psclassvars=ga_quartile year_index2017 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre 
			depressi_pre anxiety_pre antideprx_pre benzorx_pre teratrx_pre, 
	doclassvars=ga_quartile year_index4 t2dmrx_pre t1t2dmrx_pre metforrx_pre diabetes_simp_pre
			nausea_pre recurlos_pre obesity_pre chronichypertension_pre
			depressi_pre anxiety_pre antideprx_pre benzorx_pre teratrx_pre rural2,
	trtvar=trt,
	numiterations=1000,
	initialseed=23244, 
	outds=ana.senspre_boot_ptb,
	outds_dist=ga_dist_senspre_boot_ptb,
	outds_ps = ps_senspre_boot_ptb,
	outds_surv = NA
	);
	
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.senspre_boot_ptb,
			inputDist= ana.ga_dist_senspre_boot_ptb,
			numStrata= 2, 
			output_stratified= ana.senspre_ptb_boot_strat,
			output_overall= ana.senspre_ptb_boot_overall);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.senspre_boot_ptb_1, pointdsn=ana.senspre_point_ptb_1, output=ana.senspre_boot_ptb_1_ci);
%strat_estimates_w_CI(bootdsn=ana.senspre_boot_ptb_2, pointdsn=ana.senspre_point_ptb_2, output=ana.senspre_boot_ptb_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.senspre_ptb_boot_overall, pointdsn=ana.senspre_point_ptb_overall, output=ana.senspre_overall_ptb_boot_ci);

	
%count_missing_zero(inds1=ana.senspre_boot_ptb_1, inds2=ana.senspre_boot_ptb_2);













/********************************************************************************************************************************************

															05 - AD-HOC DESCRIPTIVES

********************************************************************************************************************************************/

******Identify the number of pregnancies with UNK outcomes;

*Stack the trimmed cohorts;
proc sql;
	create table trimmed as
	select *
	from ana.primary_cohort_trim_1
	union corr
	select *
	from ana.primary_cohort_trim_2
	;
	quit;


**********************************
	Look at the pregnancies with UNK outcomes;

*****PREGANNCY LOSS;

data enroll_loss;
set ana.primary_cohort;

	days_disenroll = dt_disenroll_post_any - dt_gapreg;
	daysg31 = days_disenroll > 31;

	ga_at_end = (dt_gapreg -  dt_lmp)/7;
	ga_at_index = (dt_index - dt_lmp)/7;

run;

****UNK Counts;

proc sort data=enroll_loss;
	by trt;
run;

*Overall cohort;
proc freq data=enroll_loss;
	by trt;
	table daysg31*preg_outcome_clean;
run;

*GA at index: <14 w;
proc freq data=enroll_loss (where = (ga_index_cat = 1));
	by trt;
	table daysg31*preg_outcome_clean;
run;

*GA at index: >=14 w;
proc freq data=enroll_loss (where = (ga_index_cat = 2));
	by trt;
	table daysg31*preg_outcome_clean;
run;

****GA stats;

*Overall cohort;
proc means data=enroll_loss (where = (preg_outcome_clean = 'UNK')) median p25 p75;
	class trt daysg31;
	var ga_at_index ga_at_end;
run;

*Stratified by gestational age at index;

*Overall cohort;
proc means data=enroll_loss (where = (preg_outcome_clean = 'UNK')) median p25 p75;
	class trt ga_index_cat daysg31;
	var ga_at_index ga_at_end;
run;



**Make figures;

* Gestational age at index among nifedipine initiators in the overall cohort;
proc sort data=enroll_loss (where=(trt = 1 and preg_outcome_clean = 'UNK')) out=sort_nif; 
	by daysg31; 
run;
ods graphics on / reset imagename="GA at Index - Nifedipine - Overall";
proc sgplot data=sort_nif;
    title "Gestational age at index among censored nifedipine initiators, stratified by disenrollment >31 days after censoring";
    histogram ga_at_index / group=daysg31 transparency=0.3;
    keylegend / title="Disenrolled >31d After Index"; 
    xaxis label="Gestational Age at Index (Weeks)";
    styleattrs datacolors=(navy darkgoldenrod) datacontrastcolors=(navy darkgoldenrod);
run;

* Gestational age at outcome among nifedipine initiators in the overall cohort;
ods graphics on / reset imagename="GA at Censoring - Nifedipine - Overall";
proc sgplot data=sort_nif;
    title "Gestational age at censoring among censored nifedipine initiators, stratified by disenrollment >31 days after censoring";
    histogram ga_at_end / group=daysg31 transparency=0.3;
    keylegend / title="Disenrolled >31d After Index"; 
    xaxis label="Gestational Age at Index (Weeks)";
    styleattrs datacolors=(navy darkgoldenrod) datacontrastcolors=(navy darkgoldenrod);
run;





* Gestational age at index among labetalol initiators in the overall cohort;
proc sort data=enroll_loss (where=(trt = 0 and preg_outcome_clean = 'UNK')) out=sort_lab; 
	by daysg31; 
run;
ods graphics on / reset imagename="GA at Index - Labetalol - Overall";
proc sgplot data=sort_lab;
    title "Gestational age at index among censored labetalol initiators, stratified by disenrollment >31 days after censoring";
    histogram ga_at_index / group=daysg31 transparency=0.3;
    keylegend / title="Disenrolled >31d After Index"; 
    xaxis label="Gestational Age at Index (Weeks)";
    styleattrs datacolors=(navy darkgoldenrod) datacontrastcolors=(navy darkgoldenrod);
run;

* Gestational age at outcome among nifedipine initiators in the overall cohort;
ods graphics on / reset imagename="GA at Censoring - Labetalol - Overall";
proc sgplot data=sort_lab;
    title "Gestational age at censoring among censored labetalol initiators, stratified by disenrollment >31 days after censoring";
    histogram ga_at_end / group=daysg31 transparency=0.3;
    keylegend / title="Disenrolled >31d After Index"; 
    xaxis label="Gestational Age at Index (Weeks)";
    styleattrs datacolors=(navy darkgoldenrod) datacontrastcolors=(navy darkgoldenrod);
run;






*****PRETERM BIRTH;

data enroll_ptb;
set ana.primary_cohort;

	days_disenroll = dt_disenroll_post_any - dt_gapreg;
	daysg31 = days_disenroll > 31;

	ga_at_end = (dt_gapreg_ptb -  dt_lmp)/7;
	ga_at_index = (dt_index - dt_lmp)/7;

run;

****UNK Counts;

proc sort data=enroll_ptb;
	by trt;
run;

*Overall cohort;
proc freq data=enroll_ptb;
	by trt;
	table daysg31*preg_outcome_ptb;
run;

*GA at index: <14 w;
proc freq data=enroll_ptb (where = (ga_index_cat = 1));
	by trt;
	table daysg31*preg_outcome_ptb;
run;

*GA at index: >=14 w;
proc freq data=enroll_ptb (where = (ga_index_cat = 2));
	by trt;
	table daysg31*preg_outcome_ptb;
run;

****GA stats;

*Overall cohort;
proc means data=enroll_ptb (where = (preg_outcome_ptb = 'UNK')) median p25 p75;
	class trt daysg31;
	var ga_at_index ga_at_end;
run;

*Stratified by gestational age at index;

*Overall cohort;
proc means data=enroll_ptb (where = (preg_outcome_ptb = 'UNK')) median p25 p75;
	class trt ga_index_cat daysg31;
	var ga_at_index ga_at_end;
run;





