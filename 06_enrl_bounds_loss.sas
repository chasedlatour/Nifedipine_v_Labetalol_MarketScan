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
%setup(sample=full, programname=06_enrl_bounds_loss, savelog=Y);

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

													01 - REVISE OUTCOMES AND DESCRIBE

********************************************************************************************************************************************/

*For term live births, we assume that they occur at 39w0d gestation or the LTFU date if after that.;

data partial_bounds;
set ana.primary_cohort;

	*Determine which pregnancies were LTFU due to disenrollment <= 31 days after their preg outcome;
	days_disenroll = dt_disenroll_post_any - dt_gapreg;
	daysle31 = days_disenroll <= 31;

	*****Impute outcomes for those pregnancies where preg_outcome_clean = UNK and daysle31 = 1;
	


	***Assumption: Both nifedipine and labetalol initiators IMMEDIATELY experienced the study outcome upon censoring;

	*Event indicator;
	if preg_outcome_clean = 'UNK' and daysle31 = 0 then preg_outcome_clean_both = 'SAB';
		else preg_outcome_clean_both = preg_outcome_clean;

	*New date indicator -- the same but created for easier macro;
	dt_gapreg_both = dt_gapreg;



	***Assumption: Both nifedipine and labetalol initiators NEVER experienced the study outcome and survived the full follow-up;

	*Event indicator;
	if preg_outcome_clean = 'UNK' and daysle31 = 0 then preg_outcome_clean_none = 'LBS'; /*Assume essentially that they had a live birth*/
		else preg_outcome_clean_none = preg_outcome_clean;

	*Create a new outcome date variable;
	if preg_outcome_clean = 'UNK' and daysle31 = 0 then dt_gapreg_none = max(dt_lmp + 273, dt_gapreg); /*Assume that the live birth occurred at last possible week of gestation*/
		else dt_gapreg_none = dt_gapreg;




	***Assumption: Nifedipine initiators IMMEDIATELY experienced the study outcome while labetalol initiators NEVER experienced the study outcome;

	*Event indicator;
	if trt = 1 and preg_outcome_clean = 'UNK' and daysle31 = 0 then preg_outcome_clean_nif = 'SAB';
		else if trt = 0 and preg_outcome_clean = 'UNK' and daysle31 = 0 then preg_outcome_clean_nif = 'LBS';
		else preg_outcome_clean_nif = preg_outcome_clean;

	*New event date variable;
	if trt = 0 and preg_outcome_clean = 'UNK' and daysle31 = 0 then dt_gapreg_nif = max(dt_lmp + 273, dt_gapreg);
		else dt_gapreg_nif = dt_gapreg;




	***Assumption: Labetalol initiators IMMEDIATELY experienced the study outcome while labetalol initiators NEVER experienced the study outcome;

	*Event indicator;
	if trt = 1 and preg_outcome_clean = 'UNK' and daysle31 = 0 then preg_outcome_clean_lab = 'LBS';
		else if trt = 0 and preg_outcome_clean = 'UNK' and daysle31 = 0 then preg_outcome_clean_lab = 'SAB';
		else preg_outcome_clean_lab = preg_outcome_clean;

	*New event date;
	if trt = 1 and preg_outcome_clean = 'UNK' and daysle31 = 0 then dt_gapreg_lab = max(dt_lmp + 273, dt_gapreg);
		else dt_gapreg_lab = dt_gapreg;



	drop days_disenroll;

run;

*Get descriptives of pregnancies that disenrolled <=31 days after the pregnancy outcome;
proc sort data=partial_bounds; by trt; run;
proc freq data=partial_bounds noprint;
	where preg_outcome_clean = "UNK";
	table daysle31 / out=ana.partial_bounds_count_loss;
run;

*Get counts of pregnancies that need imputed outcomes.;
proc freq data=partial_bounds;
	where preg_outcome_clean = "UNK" and daysle31 = 0;
	table trt / missing;
run;


/**Confirm that there is no unexpected missingness;*/
/*proc freq data=partial_bounds;*/
/*	table preg_outcome_clean_none preg_outcome_clean_both preg_outcome_clean_nif preg_outcome_clean_lab / missing;*/
/*run;*/
/*proc means data=partial_bounds nmiss;*/
/*	var dt_gapreg_none dt_gapreg_nif dt_gapreg_lab dt_gapreg_both;*/
/*run;*/






	
	
/********************************************************************************************************************************************

														02 - CONDUCT PREGNANCY LOSS ANALYSIS

********************************************************************************************************************************************/


%macro partial_bounds_loss(
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
		%competing2risk_weights(
			boot=0, 
			inds=partial_bounds, 
			gacatvar = ga_index_cat,
			startDT=dt_index, 
			outcomevar = preg_outcome_clean_&assumption,
			eventDT=dt_GApreg_&assumption, 
			event = 'SAB' 'UAB' 'SB' 'MLS',
			cr1='IAB', 
			cr2='LBM' 'LBS' 'UDL',
			psvars=&psvars,
			dovars=&dovars,
			dovarsmodel = &dovarsmodel,
			psclassvars=&psclassvars, 
			doclassvars=&doclassvars,
			trtvar=trt,
			numiterations=1,
			initialseed=23244, 
			outds=ana.partialb_point_&assumption,
			outds_dist=gadist_partialb_point_&assumption,
			outds_ps = NA,
			outds_surv = NA
			);

		%*Combine the stratified estimates according to the distribution of GA at index among
		those in the treated group;
		%combine_point_estimates(
			gadist = ana.gadist_partialb_point_&assumption,
			numgastrat = 2,
			est = ana.partialb_point_&assumption._,
			outds = ana.partialb_point_overall_&assumption
			);

		%*********Finally, conduct the bootstrap;
		%competing2risk_weights(
			boot=1, 
			inds=partial_bounds, 
			gacatvar = ga_index_cat,
			startDT=dt_index, 
			outcomevar = preg_outcome_clean_&assumption,
			eventDT=dt_GApreg_&assumption, 
			event = 'SAB' 'UAB' 'SB' 'MLS',
			cr1='IAB', 
			cr2='LBM' 'LBS' 'UDL',
			psvars=&psvars,
			dovars=&dovars,
			dovarsmodel = &dovarsmodel,
			psclassvars= &psclassvars, 
			doclassvars= &doclassvars,
			trtvar=trt, 
			numiterations=&numboot,
			initialseed=23244, 
			outds=ana.partialb_boot_&assumption,
			outds_dist=gadist_partialb_boot_&assumption,
			outds_ps = NA,
			outds_surv = NA
			);
			
			
		options mlogic mprint symbolgen notes;
		*Combine all the bootstrapped estimates;
		%combine_boot_estimates(
					inputEst= ana.partialb_boot_&assumption, 
					inputDist= ana.gadist_partialb_boot_&assumption,
					numStrata= 2, 
					output_stratified= ana.partialb_boot_strat_&assumption,
					output_overall= ana.partialb_boot_overall_&assumption);

		*Stratified estimates with confidence intervals.;
		%strat_estimates_w_CI(bootdsn=ana.partialb_boot_&assumption._1, pointdsn=ana.partialb_point_&assumption._1, output=ana.partialb_boot_&assumption._1_ci);
		%strat_estimates_w_CI(bootdsn=ana.partialb_boot_&assumption._2, pointdsn=ana.partialb_point_&assumption._2, output=ana.partialb_boot_&assumption._2_ci);

		*OVerall estimates with confidence interval;
		%overall_estimates_w_CI(stderrdsn=ana.partialb_boot_overall_&assumption, pointdsn=ana.partialb_point_overall_&assumption, output=ana.partialb_overall_boot_&assumption._ci);

	%end;

%mend;
	

%partial_bounds_loss(numboot=1000);


/****Assumption: Neither had an outcome;*/
%count_missing_zero(inds1=ana.partialb_boot_none_1, inds2=ana.partialb_boot_none_2);

***Both had an outcome;
%count_missing_zero(inds1=ana.partialb_boot_both_1, inds2=ana.partialb_boot_both_2);

/****Assumption: Nifedipine had an outcome;*/
%count_missing_zero(inds1=ana.partialb_boot_nif_1, inds2=ana.partialb_boot_nif_2);

/****Assumption: Labetalol had an outcome;*/
%count_missing_zero(inds1=ana.partialb_boot_lab_1, inds2=ana.partialb_boot_lab_2);




