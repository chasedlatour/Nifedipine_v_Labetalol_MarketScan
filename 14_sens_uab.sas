/********************************************************************************************************************************************
PROGRAM: 14_sens_uab.sas
PROGRAMMER: Chase Latour
PURPOSE: To conduct a sensitivity analysis where we assume that unspecified abortions (UABs) are induced abortions (IABs).
	
Goal: 
Output data: 

Date: 12.19.2024
********************************************************************************************************************************************/






/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - CONDUCT PREGNANCY LOSS ANALYSIS
	- 02 - CONDUCT PRETERM BIRTH ANALYSIS

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
%setup(sample=full, programname=14_sens_uab, savelog=N);

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

														03 - CONDUCT PREGNANCY LOSS ANALYSIS

********************************************************************************************************************************************/

*Use ana.primary_cohort as the base dataset;

*********Get point estimates;
	
*Now incorporate the IPCW;
%competing2risk_weights(
	boot=0, 
	inds=ana.primary_cohort, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_clean,
	eventDT=dt_GApreg, 
	event = 'SAB' 'SB' 'MLS',
	cr1='IAB' 'UAB', 
	cr2='LBM' 'LBS' 'UDL',
	psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			year_index2017*t2dmrx_post diabetes_simp*t2dmrx_post
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars=ga_quartile age_at_index year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre depressi_post anxiety_post antideprx_post 
			benzorx_post teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	dovarsmodel = ga_quartile age_at_index year_index4
			t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp year_index4*t2dmrx_post 
			diabetes_simp*t2dmrx_post 
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
	outds=ana.sensuab_point,
	outds_dist=ga_dist_sensuab_point,
	outds_ps = ps_sensuab_point,
	outds_surv = ana.surv_sensuab_point
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_sensuab_point,
	numgastrat = 2,
	est = ana.sensuab_point_,
	outds = ana.sensuab_point_overall
	);

*********Finally, conduct the bootstrap for SE;
	
*Now incorporate the IPCW;
%competing2risk_weights(
	boot=1, 
	inds=ana.primary_cohort, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_clean,
	eventDT=dt_GApreg, 
	event = 'SAB' 'SB' 'MLS',
	cr1='IAB' 'UAB', 
	cr2='LBM' 'LBS' 'UDL',
	psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			year_index2017*t2dmrx_post diabetes_simp*t2dmrx_post
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars=ga_quartile age_at_index year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre depressi_post anxiety_post antideprx_post 
			benzorx_post teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	dovarsmodel = ga_quartile age_at_index year_index4
			t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp year_index4*t2dmrx_post 
			diabetes_simp*t2dmrx_post 
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
	outds=ana.sensuab_boot,
	outds_dist=ga_dist_sensuab_boot,
	outds_ps = ps_sensuab_boot,
	outds_surv = NA
	);
	
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.sensuab_boot, 
			inputDist= ana.ga_dist_sensuab_boot,
			numStrata= 2, 
			output_stratified= ana.sensuab_boot_strat,
			output_overall= ana.sensuab_boot_overall);


*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.sensuab_boot_1, pointdsn=ana.sensuab_point_1, output=ana.sensuab_boot_1_ci);
%strat_estimates_w_CI(bootdsn=ana.sensuab_boot_2, pointdsn=ana.sensuab_point_2, output=ana.sensuab_boot_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.sensuab_boot_overall, pointdsn=ana.sensuab_point_overall, output=ana.sensuab_overall_boot_ci);
	


%count_missing_zero(inds1=ana.sensuab_boot_1, inds2=ana.sensuab_boot_2);

	
	
/********************************************************************************************************************************************

														04 - CONDUCT PRETERM BIRTH ANALYSIS

********************************************************************************************************************************************/
	
*********Get the point estimates;
	
*Now incorporate the IPCW;
%competing3risk_weights(
	boot=0, 
	inds=ana.primary_cohort, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_ptb,
	eventDT=dt_GApreg_ptb, 
	event = 'PTB',
	cr1='IAB' 'UAB', 
	cr2='SAB' 'SB' 'MLS',
	cr3='TB',
	psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			year_index2017*t2dmrx_post diabetes_simp*t2dmrx_post
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars=ga_quartile age_at_index year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre depressi_post anxiety_post antideprx_post 
			benzorx_post teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	dovarsmodel = ga_quartile age_at_index year_index4
			t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp year_index4*t2dmrx_post 
			diabetes_simp*t2dmrx_post 
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
	outds=ana.sensuab_point_ptb,
	outds_dist=ga_dist_sensuab_point_ptb,
	outds_ps = NA,
	outds_surv = NA
	);

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_sensuab_point_ptb,
	numgastrat = 2,
	est = ana.sensuab_point_ptb_,
	outds = ana.sensuab_point_ptb_overall
	);


*********Finally, conduct the bootstrap for SE;
	
%competing3risk_weights(
	boot=1, 
	inds=ana.primary_cohort, 
	gacatvar = ga_index_cat,
	startDT=dt_index, 
	outcomevar = preg_outcome_ptb,
	eventDT=dt_GApreg_ptb, 
	event = 'PTB',
	cr1='IAB' 'UAB', 
	cr2='SAB' 'SB' 'MLS',
	cr3='TB',
	psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			year_index2017*t2dmrx_post diabetes_simp*t2dmrx_post
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
	dovars=ga_quartile age_at_index year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
			nausea_pre recurlos_pre obesity_post chronichypertension_pre depressi_post anxiety_post antideprx_post 
			benzorx_post teratrx_pre num_outptpnc num_outptpnc_2 rural2,
	dovarsmodel = ga_quartile age_at_index year_index4
			t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp year_index4*t2dmrx_post 
			diabetes_simp*t2dmrx_post 
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
	outds=ana.sensuab_boot_ptb,
	outds_dist=ga_dist_sensuab_boot_ptb,
	outds_ps = NA,
	outds_surv = NA
	);
	
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.sensuab_boot_ptb,
			inputDist= ana.ga_dist_sensuab_boot_ptb,
			numStrata= 2, 
			output_stratified= ana.sensuab_ptb_boot_strat,
			output_overall= ana.sensuab_ptb_boot_overall);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.sensuab_boot_ptb_1, pointdsn=ana.sensuab_point_ptb_1, output=ana.sensuab_boot_ptb_1_ci);
%strat_estimates_w_CI(bootdsn=ana.sensuab_boot_ptb_2, pointdsn=ana.sensuab_point_ptb_2, output=ana.sensuab_boot_ptb_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.sensuab_ptb_boot_overall, pointdsn=ana.sensuab_point_ptb_overall, output=ana.sensuab_overall_ptb_boot_ci);

	
%count_missing_zero(inds1=ana.sensuab_boot_ptb_1, inds2=ana.sensuab_boot_ptb_2);





