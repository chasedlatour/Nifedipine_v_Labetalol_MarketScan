/********************************************************************************************************************************************
PROGRAM: 06_enrl_bounds_loss.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to conduct a sensitivity analysis where we leverage enrollment information to determine which
pregnancies were lost to follow-up due to their outcome versus disenrolling (i.e., measurable causes).
	
Goal: 
Output data: 

Date: 12.19.2024
********************************************************************************************************************************************/






/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - REVISE OUTCOMES AND DESCRIBE
	- 02 - CONDUCT PREGNANCY LOSS ANALYSIS

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
%setup(sample=full, programname=06_enrl_bounds_preec_any, savelog=Y);

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

													01 - REVISE OUTCOMES AND DESCRIBE

********************************************************************************************************************************************/	

data partial_bounds;
set ana.preeclampsia_outc;

	length preg_outcome_clean_none $4 preg_outcome_clean_both $4 preg_outcome_clean_nif $4 preg_outcome_clean_lab $4;

	*Determine which pregnancies were LTFU due to disenrollment <= 31 days after their preg outcome;
	days_disenroll = dt_disenroll_post_any - preec_any_dt; *Now use the preeclmapsia outcome date here;
	daysle31 = days_disenroll <= 31; *If equal to zero, means that MNAR;

	*****Impute outcomes for those pregnancies where preg_outcome_clean = UNK and daysle31 = 1;
	
	***Assumption: Both nifedipine and labetalol initiators IMMEDIATELY experienced the study outcome upon censoring.
	If censored prior to 20 weeks gestation, then we assume that they have preeclampsia at 20 weeks of gestation;
	
	*Event indicator;
	if preec_any_outc = 'UNK' and daysle31 = 0 then preg_outcome_clean_both = 'PRE';
		else preg_outcome_clean_both = preec_any_outc;

	*New date indicator ;
	if preec_any_outc = 'UNK' and daysle31 = 0 then dt_gapreg_both = max(dt_lmp + 140, preec_any_dt); *-- The maximum of 20w gestation or their outcome date;
		else dt_gapreg_both = preec_any_dt;


	***Assumption: Neither nifedipine nor labetalol initiators expreienced the study outcome and survived the full follow-up
		Had a healthy live birth at 39 weeks gestation and survived 14 days without the outcome. Or delivered at
		last PNC if after 39 weeks;

	*Event indicator;
	if preec_any_outc = 'UNK' and daysle31 = 0 then preg_outcome_clean_none = 'LBS';
		else preg_outcome_clean_none = preec_any_outc;

	*New date indicator;
	if preec_any_outc = 'UNK' and daysle31 = 0 then dt_gapreg_none = max(dt_lmp + 273, preec_any_dt) + 14;
		else dt_gapreg_none = preec_any_dt;


	***Assumption: Nifedipine initiators IMMEDIATELY experienced the study outcome while labetalol initiators NEVER experienced the study outcome;

	*Event indicator;
	if trt = 1 and preec_any_outc = 'UNK' and daysle31 = 0 then preg_outcome_clean_nif = 'PRE';
		else if trt = 0 and preec_any_outc = 'UNK' and daysle31 = 0 then preg_outcome_clean_nif = 'LBS';
		else preg_outcome_clean_nif = preec_any_outc;

	*New event date;
	if trt = 1 and preec_any_outc = 'UNK' and daysle31 = 0 then dt_gapreg_nif = max(dt_lmp + 140, preec_any_dt);
		else if trt = 0 and preec_any_outc = 'UNK' and daysle31 = 0 then dt_gapreg_nif = max(dt_lmp + 273, preec_any_dt) + 14;
		else dt_gapreg_nif = preec_any_dt;


	***Assumption: Labetalol initiators IMMEDIATELY experienced the study outcome while labetalol initiators NEVER experienced the study outcome;

	*Event indicator;
	if trt = 1 and preec_any_outc = 'UNK' and daysle31 = 0 then preg_outcome_clean_lab = 'LBS';
		else if trt = 0 and preec_any_outc = 'UNK' and daysle31 = 0 then preg_outcome_clean_lab = 'PRE';
		else preg_outcome_clean_lab = preec_any_outc;

	*New event date;
	if trt = 1 and preec_any_outc = 'UNK' and daysle31 = 0 then dt_gapreg_lab = max(dt_lmp + 273, preec_any_dt) + 14;
		else if trt = 0 and preec_any_outc = 'UNK' and daysle31 = 0 then dt_gapreg_lab = max(dt_lmp + 140, preec_any_dt);
		else dt_gapreg_lab = preec_any_dt;


	drop days_disenroll;

run;



*Get descriptives of pregnancies that disenrolled <=31 days after the pregnancy outcome;
proc sort data=partial_bounds; by trt; run;
proc freq data=partial_bounds noprint;
	where preg_outcome_clean = "UNK";
	table daysle31 / out=ana.partial_bounds_count_pany;
run;



/**Confirm that there is no unexpected missingness;*/
proc freq data=partial_bounds;
	table preg_outcome_clean_none preg_outcome_clean_both preg_outcome_clean_nif preg_outcome_clean_lab / missing;
run;
proc means data=partial_bounds nmiss;
	var dt_gapreg_none dt_gapreg_nif dt_gapreg_lab dt_gapreg_both;
run;






	
	
/********************************************************************************************************************************************

														02 - CONDUCT PRETERM BIRTH ANALYSIS

********************************************************************************************************************************************/


%macro partial_bounds_pany(
		allassumptions= none both nif lab,
		numboot = 2,
		psvars=ga_quartile age_at_index age_at_index_2
					year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
					nausea_pre recurlos_pre obesity_post chronichypertension_pre
					depressi_post anxiety_post antideprx_post benzorx_post 
					teratrx_pre num_outptpnc num_outptpnc_2,
		dovars=ga_quartile age_at_index year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
					nausea_pre recurlos_pre obesity_post chronichypertension_pre depressi_post anxiety_post antideprx_post 
					benzorx_post teratrx_pre num_outptpnc num_outptpnc_2 rural2,
		dovarsmodel=ga_quartile age_at_index year_index4
					t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
					nausea_pre recurlos_pre obesity_post chronichypertension_pre 
					depressi_post anxiety_post antideprx_post benzorx_post
					teratrx_pre num_outptpnc num_outptpnc_2 rural2,
		psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
					nausea_pre recurlos_pre obesity_post chronichypertension_pre 
					depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre,
		doclassvars=ga_quartile year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
					nausea_pre recurlos_pre obesity_post chronichypertension_pre
					depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre rural2);

	%*Count the number of assumptions;
	%let numassumptions = %sysfunc(countw(&allassumptions));

	%*Repeat the analysis over each assumption;

	%do n=1 %to &numassumptions;

		%let assumption = %scan(&allassumptions, &n);

		%*Gather the point estimate;
		%competing3risk_weights(
			boot=0, 
			inds=partial_bounds, 
			gacatvar = ga_index_cat,
			startDT=dt_index, 
			outcomevar = preg_outcome_clean_&assumption,
			eventDT=dt_GApreg_&assumption, 
			event = 'PRE',
			cr1='IAB', 
			cr2= 'SAB' 'UAB' 'SB' 'MLS',
			cr3= 'LBS' 'LBM' 'UDL',
			psvars=&psvars,
			dovars=&dovars,
			dovarsmodel = &dovarsmodel,
			psclassvars=&psclassvars, 
			doclassvars=&doclassvars,
			trtvar=trt,
			numiterations=1,
			initialseed=23244, 
			outds=ana.partialb_point_pany_&assumption,
			outds_dist=gadist_partialb_point_pany_&assumption,
			outds_ps = NA,
			outds_surv = NA
			);

		%*Combine the stratified estimates according to the distribution of GA at index among
		those in the treated group;
		%combine_point_estimates(
			gadist = ana.gadist_partialb_point_pany_&assumption,
			numgastrat = 2,
			est = ana.partialb_point_pany_&assumption._,
			outds = ana.partialb_point_pany_overall_&assumption
			);

		%*********Finally, conduct the bootstrap;
		%competing3risk_weights(
			boot=1, 
			inds=partial_bounds, 
			gacatvar = ga_index_cat,
			startDT=dt_index, 
			outcomevar = preg_outcome_clean_&assumption,
			eventDT=dt_GApreg_&assumption, 
			event = 'PRE',
			cr1='IAB', 
			cr2='SAB' 'UAB' 'SB' 'MLS',
			cr3='LBS' 'LBM' 'UDL',
			psvars=&psvars,
			dovars=&dovars,
			dovarsmodel = &dovarsmodel,
			psclassvars= &psclassvars, 
			doclassvars= &doclassvars,
			trtvar=trt, 
			numiterations=&numboot,
			initialseed=23244, 
			outds=ana.partialb_boot_pany_&assumption,
			outds_dist=gadist_partialb_boot_pany_&assumption,
			outds_ps = NA,
			outds_surv = NA
			);
			
			
		options mlogic mprint symbolgen notes;
		*Combine all the bootstrapped estimates;
		%combine_boot_estimates(
					inputEst= ana.partialb_boot_pany_&assumption, 
					inputDist= ana.gadist_partialb_boot_pany_&assumption,
					numStrata= 2, 
					output_stratified= ana.partialb_boot_pany_strat_&assumption,
					output_overall= ana.partialb_boot_pany_overall_&assumption);


		*Stratified estimates with confidence intervals.;
		%strat_estimates_w_CI(bootdsn=ana.partialb_boot_pany_&assumption._1, pointdsn=ana.partialb_point_pany_&assumption._1, output=ana.partialb_boot_pany_&assumption._1_ci);
		%strat_estimates_w_CI(bootdsn=ana.partialb_boot_pany_&assumption._2, pointdsn=ana.partialb_point_pany_&assumption._2, output=ana.partialb_boot_pany_&assumption._2_ci);

		*OVerall estimates with confidence interval;
		%overall_estimates_w_CI(stderrdsn=ana.partialb_boot_pany_overall_&assumption, pointdsn=ana.partialb_point_pany_overall_&assumption, output=ana.partialb_all_boot_pany_&assumption._ci);

	%end;

%mend;
	

%partial_bounds_pany(numboot=1000); 





/****Assumption: Neither had an outcome;*/
%count_missing_zero(inds1=ana.partialb_boot_pany_none_1, inds2=ana.partialb_boot_pany_none_2);

***Both had an outcome;
%count_missing_zero(inds1=ana.partialb_boot_pany_both_1, inds2=ana.partialb_boot_pany_both_2);

/****Assumption: Nifedipine had an outcome;*/
%count_missing_zero(inds1=ana.partialb_boot_pany_nif_1, inds2=ana.partialb_boot_pany_nif_2);

/****Assumption: Labetalol had an outcome;*/
%count_missing_zero(inds1=ana.partialb_boot_pany_lab_1, inds2=ana.partialb_boot_pany_lab_2);
