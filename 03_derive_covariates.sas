/********************************************************************************************************************************************
PROGRAM: 03_derive_covariates.sas
PROGRAMMER: Chase Latour
PURPOSE: To derive covariates necessary for the analyses

Date: 12.19.2024

Updates:
	- 12.4.2025 - CDL: Added additional covariates that derived for revisions (anemia, GI disease, other HTN disorders, uterine fibroids, 
					   aspirin, anti-inflammatory medications
********************************************************************************************************************************************/






/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - IDENTIFY VARIABLES - STILL NEED TO RUN SOME DESCRIPTIVES
	- 02 - SUBSET TO ANALYTIC COHORT WITH DESCRIPTIVES -- NEED TO DO

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
%setup(sample=full, programname=03_derive_covariates, savelog=N);

******************************************************************************************************************************************;
/*Create local mirros of the libraries from the set up macro - Run locally*/
/*libname lout slibref=out server=server; */
/*libname lana slibref=ana server=server;*/
/*libname lpreg slibref=preg server=server;*/
/*libname lwork slibref=work server=server;*/
/*libname ltemp slibref=temp server=server;*/
/*libname lraw slibref=raw server=server;*/
/*libname lder slibref=der server=server;*/
/*libname lcovref slibref=covref server=server;*/
/*libname lexpref slibref=expref server=server;*/
/*libname lrxcov slibref=rxcov server=server;*/






/********************************************************************************************************************************************

															01 - IDENTIFY VARIABLES

********************************************************************************************************************************************/

/*
MACRO: covar
PURPOSE: To identify covariates (not incl, excl criteria) based upon predefined definitions


INPUTS:
INPUT_DATA -- input dataset of pregnancies
OUTPUT_DATA -- output dataset of pregnancies
*/


%macro covar(INPUT_DATA=pregnancies2, OUTPUT_DATA = pregnancies3);

	data &OUTPUT_DATA;
	format year_index $4.;
	set &INPUT_DATA;
	
		%*Create a year indicator for index date;
		year_index = year(dt_index);
		year_index_numeric = year(dt_index);
		%*Binary indicator for index date;
		if year(dt_index) <= 2019 then year_le2019 = 1;
			else year_le2019 = 0;

		%*Rurality;
		if MSA_at_index = 0 then rural = 1;
			else if MSA_at_index = . then rural = .;
			else rural = 0;

		%*****Opioid use disorder;

		%*Use pre-index information only;
		if sum(dx_OudOutptpre, dx_OudInptpre) >= 1 then oud_pre = 1;
			else if pr_Oudpre >= 1 then oud_pre = 1;
			else if rx_oudpre >=1 then oud_pre = 1;
			else oud_pre = 0;

		%*Include post-index information;
		if oud_pre = 1 then oud_post = 1;
			else if sum(dx_OudOutptpost, dx_OudInptpost) >= 1 then oud_post = 1;
			else if pr_oudpos >= 1 then oud_post = 1;
			else if rx_oudpos >=1 then oud_post = 1;
			else oud_post = 0;

		%*********Smoking;

		%*Use pre-index information only;
		if sum(dx_SmkOutptPre, dx_SmkInptPre) >= 1 then Smk_pre = 1;
			else if rx_SmkPre >= 1 then Smk_pre = 1;
			else Smk_pre = 0;

		%*Include post-index information;
		if Smk_pre = 1 then Smk_post = 1;
			else if sum(dx_SmkOutptPost, dx_SmkInptPost) >= 1 then Smk_post = 1;
			else if rx_SmkPos >= 1 then Smk_post = 1;
			else Smk_post = 0;

		%**********Covariates defined by only one diagnosis code;

		%do d=1 %to &num_one_dx;

			%let loop&d = %scan(&one_dx, &d);

			%*Only use pre-index information;
			if dx_&&loop&d..Outptpre >= 1 then &&loop&d.._pre = 1;
				else if dx_&&loop&d..Inptpre >= 1 then &&loop&d.._pre = 1;
				else &&loop&d.._pre = 0;

			%*Include post-index information;
			if &&loop&d.._pre = 1 then &&loop&d.._post = 1;
				else if dx_&&loop&d..Outptpost >= 1 then &&loop&d.._post = 1;
				else if dx_&&loop&d..Inptpost >= 1 then &&loop&d.._post = 1;
				else &&loop&d.._post = 0;

		%end;

		%******Nausea and vomiting;

		%*Only use information up to index date;
		if dx_NauseaOutptpre >= 1 and rx_Nauseapre >=1 then nausea_pre = 1;
			else if dx_NauseaInptpre >= 1 and rx_Nauseapre >= 1 then nausea_pre = 1;
			else if dx_ps_NauseaInptpre >= 1 then nausea_pre = 1;
			else nausea_pre = 0;

		%*****Diabetes -- derive the initial indicator here;
		
		%*Using only pre-index information; 
		%*CDL: REVISED 11.16.2025 because indicator was supposed to be at least 2 codes of any types on different dates;
		if sum(dx_DiabDMge141idx, dx_DiabDMle140idx) >= 2 then any_diabetes_pre = 1;
/*		if dx_DiabetesOutptPre >= 2 then any_diabetes_pre = 1;*/
/*			else if dx_DiabetesInptpre >= 1 then any_diabetes_pre = 1;*/
		*CDL: FIXED below code on 1.15.2026 so that referencing the same counts.;
/*			else if sum(dx_DiabetesOutptPre, dx_DiabetesInptpre) >= 1 and */
			else if sum(dx_DiabDMge141idx, dx_DiabDMle140idx) >= 1 and 
					sum(rx_t2dmPre, rx_t1t2dmpre)>= 1 then any_diabetes_pre = 1;
			else if sum(rx_t2dmpre, rx_t1t2dmpre)>= 2 then any_diabetes_pre = 1;
			else any_diabetes_pre = 0;

		%*Including post-index information; 
		%*CDL: REVISED 11.16.2025 because indicator was supposed to be at least 2 codes of any types on different dates;
		if sum(dx_DiabDMge141idx30, dx_DiabDMle140idx30) >= 2 then any_diabetes_post = 1;
			else if sum(dx_DiabDMge141idx30, dx_DiabDMle140idx30) >= 1 and 
					sum(rx_t2dmPre, rx_t2dmPos, rx_t1t2dmpre, rx_t1t2dmpos) >= 1 then any_diabetes_post = 1;
/*		if any_diabetes_pre_dx = 1 then any_diabetes_post = 1;*/
/*			else if dx_DiabetesInptPost >= 1 then any_diabetes_post = 1;*/
/*			else if sum(dx_DiabetesOutptPre, dx_DiabetesOutptPost) >= 2 then any_diabetes_post = 1;*/
/*			else if sum(dx_DiabetesOutptPre, dx_DiabetesOutptPost) >= 1 and */
/*					sum(rx_t2dmPre, rx_t2dmPos, rx_t1t2dmpre, rx_t1t2dmpos) >= 1 then any_diabetes_post = 1;*/
			else if sum(rx_t2dmpre, rx_t2dmpos, rx_t1t2dmpre, rx_t1t2dmpos) >= 2 then any_diabetes_post = 1;
			else any_diabetes_post = 0;
			
		%******Diabetes -- derive the initial indicator but require at least 1 dx code;
			
		%*Using only pre-index information;
		if dx_DiabetesOutptPre >= 2 then any_diabetes_pre_dx = 1;
			else if dx_DiabetesInptpre >= 1 then any_diabetes_pre_dx = 1;
			else if (dx_DiabetesOutptPre >= 1) and 
					sum(rx_t2dmPre, rx_t1t2dmpre)>= 1 then any_diabetes_pre_dx = 1;
/*			else if rx_t1t2dmpre>= 2 then any_diabetes_pre_dx = 1;*/
			else if sum(rx_t2dmPre, rx_t1t2dmpre)>= 2 then any_diabetes_pre_dx = 1; %*CDL: CHANGED 2.18.25 based on MGP code review.;
			else any_diabetes_pre_dx = 0;


		%*Including post-index information;
		if any_diabetes_pre = 1 then any_diabetes_post_dx = 1;
			else if sum(dx_DiabetesOutptPre, dx_DiabetesOutptPost) >= 2 then any_diabetes_post_dx = 1;
			else if dx_DiabetesInptPost >= 1 then any_diabetes_post_dx = 1;
			else if sum(dx_DiabetesOutptPre, dx_DiabetesOutptPost) >= 1 and 
					sum(rx_t2dmPre, rx_t2dmPos, rx_t1t2dmpre, rx_t1t2dmpos) >= 1 then any_diabetes_post_dx = 1;
/*			else if sum(rx_t1t2dmpre, rx_t1t2dmpos) >= 2 then any_diabetes_post_dx = 1;*/
			else if sum(rx_t2dmPre, rx_t2dmPos, rx_t1t2dmpre, rx_t1t2dmpos) >= 2 then any_diabetes_post_dx = 1; %*CDL: CHANGED 2.18.25 based on MGP code review;
			else any_diabetes_post_dx = 0;


		%******Variables defined by ONE outpatient fill;
		%do r=1 %to &num_one_rx;

			%let loop&r = %scan(&one_rx, &r);

			%*Only use pre-index information;
			if rx_&&loop&r..Pre >= 1 then &&loop&r..rx_pre = 1;
				else &&loop&r..rx_pre = 0;

			%*Include post-index information;
			if rx_&&loop&r..Pre >= 1 then &&loop&r..rx_post = 1;
				else if rx_&&loop&r..Pos >= 1 then &&loop&r..rx_post = 1;
				else &&loop&r..rx_post = 0;

		%end;


		%******Migraine;

		%*Use only pre-index information;
		if dx_MigraineOutptPre >= 2 then migraine_pre = 1;
			else if dx_MigraineInptPre >= 1 then migraine_pre = 1;
			else if dx_MigraineOutptPre >= 1 and rx_MigrainePre >= 1 then migraine_pre = 1;
			else migraine_pre = 0;

		%*Include post-index information;
		if sum(dx_MigraineOutptPre, dx_MigraineInptPre) >= 2 then migraine_post = 1;
			else if dx_MigraineInptPre >= 1 then migraine_post = 1;
			else if dx_MigraineInptPost >= 1 then migraine_post = 1;
			else if (sum(dx_MigraineOutptPre, dx_MigraineOutptPost) >= 1) and 
					(sum(rx_MigrainePre, rx_MigrainePos) >= 1) then migraine_post = 1;
			else migraine_post = 0;



		%******Covariates defined by 2 outpatient or 1 inpatient dx code;

		%*Deal with all of those diagnoses defined by only one diagnosis code;
		%do o=1 %to &num_two_out_one_in;

			%let loop&o = %scan(&two_out_one_in, &o);

			%*Only use pre-index information;
			if dx_&&loop&o..Outptpre >= 2 then &&loop&o.._pre = 1;
				else if dx_&&loop&o..Inptpre >= 1 then &&loop&o.._pre = 1;
				else &&loop&o.._pre = 0;

			%*Include post-index information;
			if dx_&&loop&o..Outptpre >= 2 then &&loop&o.._post = 1;
				else if dx_&&loop&o..Outptpost >= 2 then &&loop&o.._post = 1;
				else if dx_&&loop&o..Outptpre >= 1 and dx_&&loop&o..Outptpost >= 1 then &&loop&o.._post = 1;
				else if dx_&&loop&o..Inptpre >= 1 then &&loop&o.._post = 1;
				else if dx_&&loop&o..Inptpost >= 1 then &&loop&o.._post = 1;
				else &&loop&o.._post = 0;

		%end;


		%******Anxiety;

		%*Use only information prior to index;
		if dx_anxietyOutptPre >= 2 then anxiety_pre = 1;
			else if dx_anxietyInptpre >= 1 then anxiety_pre = 1;
			else if dx_anxietyOutptPre >= 1 and rx_anxietyPre >= 2 then anxiety_pre = 1;
			else anxiety_pre = 0;

		%*Include information after index;
		if anxiety_pre = 1 then anxiety_post = 1;
			else if sum(dx_anxietyOutptPre, dx_anxietyOutptPost) >= 2 then anxiety_post = 1;
			else if dx_anxiety_inptPost >= 1 then anxiety_post = 1;
			else if sum(dx_anxietyOutptPre, dx_anxietyOutptPost) >=1 and 
				sum(rx_anxietyPre, rx_anxietyPos) >= 2 then anxiety_post = 1;
			else anxiety_post = 0;

		%********ADHD;

		%*Use only information prior to index;
		if sum(dx_adhdOutptpre, dx_adhdInptpre) >= 1 and rx_ADHDpre >=1 then adhd_pre = 1;
			else if sum(dx_adhdOutptpre, dx_adhdInptpre) >= 2 then adhd_pre = 1;
			else adhd_pre = 0;

		%*Include information after index;
		if adhd_pre = 1 then adhd_post = 1;
			else if sum(dx_adhdOutptpre, dx_adhdInptpre, dx_adhdOutptpost, dx_adhdInptpost) >= 1 and
					sum(rx_ADHDpre, rx_ADHDpos) >= 1 then adhd_post = 1;
			else if sum(dx_adhdOutptpre, dx_adhdInptpre, dx_adhdOutptpost, dx_adhdInptpost) >= 2 then adhd_post = 1;
			else adhd_post = 0;


		%********Covariates with >=1 inpt or >=2 outpt where >=1 outpt must occur prior to index.;
		%do p=1 %to &num_two_one_prior;

			%let loop&p = %scan(&two_one_prior, &p);

			%*Only use pre-index information;
			if dx_&&loop&p..Outptpre >= 2 then &&loop&p.._pre = 1;
				else if dx_&&loop&p..Inptpre >= 1 then &&loop&p.._pre = 1;
				else &&loop&p.._pre = 0;

			%*Include post-index information;
			if dx_&&loop&p..Outptpre >= 2 then &&loop&p.._post = 1;
				else if dx_&&loop&p..Outptpre >= 1 and dx_&&loop&p..Outptpost >= 1 then &&loop&p.._post = 1; /*one outpt must occur prior to index*/
				else if dx_&&loop&p..Inptpre >= 1 then &&loop&p.._post = 1; /*inpatent only considered before index*/
				else &&loop&p.._post = 0;

		%end;


		%******Hyperlipidemia;

		%*Use only information up to index;
		if dx_HyperLipOutptPre >= 2 then hyperlip_pre = 1;
			else if dx_HyperLipInptPre >= 1 then hyperlip_pre = 1;
			else if rx_statinpre >= 1 then hyperlip_pre = 1;
			else hyperlip_pre = 0;

		%*Include information after index;
		if hyperlip_pre = 1 then hyperlip_post = 1;
			else if sum(dx_HyperLipOutptPre, dx_HyperLipOutptPost) >= 2 then hyperlip_post = 1;
			else if dx_HyperLipInptPost >=1 then hyperlip_post = 1;
			else if rx_statinpos >=1 then hyperlip_post = 1;
			else hyperlip_post = 0;
			
		%******Superimposed preeclampsia - only use information prior to index_date;
		if dx_PreecOutptpre >= 1 then preeclampsia_pre = 1;
			else if dx_PreecInptpre >= 1 then preeclampsia_pre = 1;
			else preeclampsia_pre = 0;
			
		%******Bariatric surgery;
		
		%*Only pre-index information;
		if sum(pr_Bariatricpre, pr_BariatricRpre) >= 1 then bariatric_pre = 1;
			else if sum(dx_BariatricOutptpre) >= 2 then bariatric_pre = 1;
			else if dx_BariatricInptpre >= 1 then bariatric_pre = 1;
			else bariatric_pre = 0;
			
		%*Include infomration up to 30 days post-index;
		if bariatric_pre = 1 then bariatric_post = 1;
			else if sum(dx_BariatricOutptpre, dx_BariatricOutptpost) >= 2 then bariatric_post = 1;
			else if dx_BariatricInptpost = 1 then bariatric_post = 1;
			else bariatric_post = 0;

	run;

%mend;




/*
MACRO: identify_variables
PURPOSE: To identify key variables for inclusion and exclusion criteria as well as confounding using prespecified variable
definitions

INPUTS: 
LOOKBACKDT - Date from which we look back for prior covariate information: dt_index or dt_lmp
LOOKBACKDAYS - The number of days we look back from LOOKBACKDT
LMPINDEX - Assumed GA at the indexing prenatal encounter
*/

%macro identify_variables(LOOKBACKDT, LOOKBACKDAYS, LMPINDEX, SUFFIX);

	/*TESTING:
	%let lookbackdt = dt_index;
	%let lookbackdays = 270;
	%let lmpindex = 63;
	*/

	%*Read in the pregnancy cohort;
	data pregnancies;
	set 
		%if &suffix = NA %then %do;
			out.preg_cohort_&lmpindex._&lookbackdt._&lookbackdays;
		%end;
			%*CDL: ADDED;
			%else %if &suffix = DEL %then %do;
				out.preg_del_cohort_&lmpindex._&lookbackdt._&lookbackdays;
			%end;
			%else %do;
				out.cohort_&lmpindex._&lookbackdt._&lookbackdays._&suffix;
			%end;
	run;
	
/* 	%*Temporary; */
/* 	data pregnancies; */
/* 	set temp.preg_covar_&lmpindex._&lookbackdt._&lookbackdays; */
/* 	run; */


	%*Identify variables use for inclusion and exclusion criteria a priori;

	data pregnancies2;
	set pregnancies;

		%******Gestational age at index date;
		ga_at_index = dt_index - dt_lmp;

		%******Chronic Hypertension - Using only information pre-index;
		if dx_ChHtnOutptpre >= 1 then chronichypertension_pre = 1;
			else if dx_ChHtnInptpre >= 1 then chronichypertension_pre = 1;
			else chronichypertension_pre = 0;

		%******Cancer;
		if dx_CancerOutptpre >= 2 then cancer_pre = 1;
			else if dx_CancerInptpre >= 1 then cancer_pre = 1;
			else cancer_pre = 0;

		%******Pre-pregnancy asthma - Using only information pre-index;
		if dx_30_AsthmaOutptpre >= 2 then asthma_pre = 1;
			else if dx_AsthmaInptpre >= 1 then asthma_pre = 1;
			else if dx_AsthmaOutptpre >= 1 and Rx_AsthmaPre >= 1 then asthma_pre = 1;
			else asthma_pre = 0;

		%*******Coronary heart disease;
		if dx_CADOutptpre >= 2 then coronary_heart_disease_pre = 1;
			else if dx_CADInptpre >=1 then coronary_heart_disease_pre = 1;
			else coronary_heart_disease_pre = 0;

		%*******Arrhythmias;
		if dx_ArrhyOutptpre >= 1 then arrhythmia_pre = 1;
			else if dx_ArrhyInptpre >= 1 then arrhythmia_pre = 1;
			else arrhythmia_pre = 0;

		%*******Congenital heart defects;
		if first_congenitalheartdx ne . and (first_fetalecho = . or (first_fetalecho ne . and first_congenitalheartdx < first_fetalecho))
				then congenital_heart_pre = 1;
			else congenital_heart_pre = 0;

		%*******Endocarditis;
		if dx_EndocardOutptpre >= 2 then endocarditis_pre = 1;
			else if dx_ps_EndocardInptpre >= 1 then endocarditis_pre = 1;
			else endocarditis_pre = 0;

		%*******Myo-Pericardiditis;
		if dx_MyoPeriOutptpre >= 2 then myopericarditis_pre = 1;
			else if dx_ps_MyoPeriInptpre >= 1 then myopericarditis_pre = 1;
			else myopericarditis_pre = 0;

		%*******Heart Failure;
		if dx_HFOutptpre >= 2 then heartfailure_pre = 1;
			else if dx_HFInptpre >= 1 then heartfailure_pre = 1;
			else heartfailure_pre = 0;

		%*******Heart valve disease;
		if dx_HValveOutptpre >= 2 then heart_valve_disease_pre = 1;
			else if dx_HValveInptpre >= 1 then heart_valve_disease_pre = 1;
			else heart_valve_disease_pre = 0;

		%*******Cardiomyopathy;
		if dx_CMyoOutptpre >=2 then cardiomyopathy_pre = 1;
			else if dx_CMyoInptpre >=1 then cardiomyopathy_pre = 1;
			else cardiomyopathy_pre = 0;

		%*******Other heart diseases;
		if dx_OtherHrtOutptpre >=2 then other_heart_disease_pre = 1;
			else if dx_OtherHrtInptpre >=1 then other_heart_disease_pre = 1;
			else other_heart_disease_pre = 0;

	run;


	%let one_dx = alc othersud RecurLos ckd obesity Bipolar Schizo fibroids gidiseas pulmonar anemia;
	%let num_one_dx = %sysfunc(countw(&one_dx));

	%*CDL: ADDED simulant 4.24.2025 -- NOTE stimulant was mispelled as simulant in the original coding and so retain that error here.;
	%let one_rx = t2dm t1t2dm metfor HyperThy HypoThy terat statin benzo antidep anticonvul MoodStab adhd simulant PTSD Antipsy glp1wgt otherwgt lda nsaid aspirin;
	%let num_one_rx = %sysfunc(countw(&one_rx));

	%let two_out_one_in = retino Antiphos Lupus HyperThy HypoThy depressi PTSD Athero PeriVasc Anemia SickleT SickleD;
	%let num_two_out_one_in = %sysfunc(countw(&two_out_one_in));

	%let two_one_prior = stroke mi angina;
	%let num_two_one_prior = %sysfunc(countw(&two_one_prior));


	%*****Identify covariates
		All variable names possible are included even though not all may be there. may leave warnings.;
	options mprint;
	%covar(input_data=pregnancies2, output_data=pregnancies3);

	%*****Take additional steps to identify diabetes types from any_diabetes_pre and _post;

	%*Determine if the person has pregestational diabetes: if they satisfy the above conditions based upon information 
	prior to the LMP up to 90 days of gestation;
	data pregnancies3_pregest;
	set pregnancies3;

		%*Using only pre-index information;
		if any_diabetes_pre = 1 then do;
			if dx_DiabPreGOutptIdx >=2 then pregestation_diab_pre = 1;
				else if dx_DiabPreGInptIdx >= 1 then pregestation_diab_pre = 1;
				else if dx_DiabPreGOutptIdx >= 1 and rx_nonmet_antidiable90GAidx then pregestation_diab_pre = 1;
				else if rx_nonmet_antidiable90GAidx >= 2 then pregestation_diab_pre = 1;
				else pregestation_diab_pre = 0;
		end;
		else pregestation_diab_pre = 0;

		%*Including up to 30 days post-index;
		if any_diabetes_post = 1 then do;
			if pregestation_diab_pre = 1 then pregestation_diab_post = 1;
/*				else if sum(dx_DiabPreGOutptIdx, dx_DiabPreGOutptIdx30) >= 2 then pregestation_diab_post = 1;*/
				%*CDL: CHANGED 2.18.25 based on MGP code review;
				else if dx_DiabPreGOutptIdx30 >= 2 then pregestation_diab_post = 1;
				else if dx_DiabPreGInptIdx30 >= 1 then pregestation_diab_post = 1;
				%*CDL: CHANGED 2.18.25 based on MGP code review;
				else if dx_DiabPreGOutptIdx30 >= 1 and rx_nonmet_antidiable90GAidx30 >= 1 then pregestation_diab_post = 1;
/*				else if sum(dx_DiabPreGOutptIdx, dx_DiabPreGOutptIdx30) >= 1 and */
/*						sum(rx_nonmet_antidiable90GAidx, rx_nonmet_antidiable90GAidx30) >=1 then pregestation_diab_post = 1;*/
				else if sum(rx_nonmet_antidiable90GAidx, rx_nonmet_antidiable90GAidx30) >= 2 then pregestation_diab_post = 1;
				else pregestation_diab_post = 0;
		end;
		else pregestation_diab_post = 0;
		
		%*Require daignosis code; 
		
		%*Using only pre-index information;
		if any_diabetes_pre_dx = 1 then do;
			if dx_DiabPreGOutptIdx >=2 then pregestation_diab_pre_dx = 1;
				else if dx_DiabPreGInptIdx >= 1 then pregestation_diab_pre_dx = 1;
				else if dx_DiabPreGOutptIdx >= 1 and rx_nonmet_antidiable90GAidx >= 1 then pregestation_diab_pre_dx = 1; %*CDL: ADDED second >=1 based on MGP code review 2.18.25;
				else if rx_nonmet_antidiable90GAidx >= 2 then pregestation_diab_pre_dx = 1;
				else pregestation_diab_pre_dx = 0;
		end;
		else pregestation_diab_pre_dx = 0;

		%*Including up to 30 days post-index;
		if any_diabetes_post_dx = 1 then do;
			if pregestation_diab_pre_dx = 1 then pregestation_diab_post_dx = 1;
				%*CDL: CHANGED 2.18.25 based on MGP code review;
/*				else if sum(dx_DiabPreGOutptIdx, dx_DiabPreGOutptIdx30) >= 2 then pregestation_diab_post_dx = 1;*/
				else if dx_DiabPreGOutptIdx30 >= 2 then pregestation_diab_post_dx = 1;
				else if dx_DiabPreGInptIdx30 >= 1 then pregestation_diab_post_dx = 1;
				%*CDL: CHANGED 2.18.25 based on MGP code review;
/*				else if sum(dx_DiabPreGOutptIdx, dx_DiabPreGOutptIdx30) >= 1 and */
/*						sum(rx_nonmet_antidiable90GAidx, rx_nonmet_antidiable90GAidx30) >=1 then pregestation_diab_post_dx = 1;*/
				else if dx_DiabPreGOutptIdx30 >= 1 and rx_nonmet_antidiable90GAidx30 >= 1 then pregestation_diab_post_dx = 1;
				%*CDL: CHANGED 2.18.25 based on MGP code review;
/*				else if sum(rx_nonmet_antidiable90GAidx, rx_nonmet_antidiable90GAidx30) >= 2 then pregestation_diab_post_dx = 1;*/
				else if rx_nonmet_antidiable90GAidx30 >= 2 then pregestation_diab_post_dx = 1;
				else pregestation_diab_post_dx = 0;
		end;
		else pregestation_diab_post_dx = 0;

	run;

	%*Distinguish the type of diabetes among those with pregestational diabetes;
	data pregnancies3_diab;
	set pregnancies3_pregest;
		length diabetes_type_pre $4 diabetes_type_post $4;
		
		%*******Look at only pre-index information;		%*Those without pregestational diabetes;
		if any_diabetes_pre = 0 then diabetes_type_pre = "NA";
			
			else if pregestation_diab_pre = 0 then diabetes_type_pre = "TBD";

			%*Those with pregestational diabetes;
			else do;
				
				%*T1DM;
				if dx_DiabPreGT1DMIdx >= 1 and
					dx_DiabPreGT2DMIdx < 2 and 
					rx_t2dm_antidiable90GAidx <= 0 and
					rx_insulinle90GAidx >= 1 then diabetes_type_pre = "T1DM";

				%*T2DM;
				else if rx_t2dm_antidiable90GAidx >= 2 or
							rx_insulinle90GAidx = 0 or
							(dx_DiabPreGT2DMIdx >= 2 and dx_DiabPreGT1DMIdx < 2 and rx_insulinle90GAidx >= 1) 
						then diabetes_type_pre = "T2DM";


				%*Not otherwise specified;
				else diabetes_type_pre = "NOS";
				
			end;

		%*Include up to 30d post-index information;
		if any_diabetes_post = 0 then diabetes_type_post = "NA";
			
			else if pregestation_diab_post = 0 then diabetes_type_post = "TBD";

			%*Those with pregestational diabetes;
			else do;
				
				%*T1DM;
				if dx_DiabPreGT1DMIdx30 >= 1 and
					dx_DiabPreGT2DMIdx30 < 2 and 
					rx_t2dm_antidiable90GAidx30 <= 0 and
					rx_insulinle90GAidx30 >= 1 then diabetes_type_post = "T1DM";

				%*T2DM;
				else if rx_t2dm_antidiable90GAidx30 >= 2 or
							rx_insulinle90GAidx30 = 0 or
							(dx_DiabPreGT2DMIdx30 >= 2 and dx_DiabPreGT1DMIdx30 < 2 and rx_insulinle90GAidx30 >= 1) 
						then diabetes_type_post = "T2DM";


				%*Not otherwise specified;
				else diabetes_type_post = "NOS";
				
			end;
			
		%*******Require diagnosis code;	
			
		%*******Look at only pre-index information;		%*Those without pregestational diabetes;
		if any_diabetes_pre_dx = 0 then diabetes_type_pre_dx = "NA";
			
			else if pregestation_diab_pre_dx = 0 then diabetes_type_pre_dx = "TBD";

			%*Those with pregestational diabetes;
			else do;
				
				%*T1DM;
				if dx_DiabPreGT1DMIdx >= 1 and
					dx_DiabPreGT2DMIdx < 2 and 
					rx_t2dm_antidiable90GAidx <= 0 and
					rx_insulinle90GAidx >= 1 then diabetes_type_pre_dx = "T1DM";

				%*T2DM;
				else if rx_t2dm_antidiable90GAidx >= 2 or
							rx_insulinle90GAidx = 0 or
							(dx_DiabPreGT2DMIdx >= 2 and dx_DiabPreGT1DMIdx < 2 and rx_insulinle90GAidx >= 1) 
						then diabetes_type_pre_dx = "T2DM";


				%*Not otherwise specified;
				else diabetes_type_pre_dx = "NOS";
				
			end;

		%*Include up to 30d post-index information;
		if any_diabetes_post_dx = 0 then diabetes_type_post_dx = "NA";
			
			else if pregestation_diab_post_dx = 0 then diabetes_type_post_dx = "TBD";

			%*Those with pregestational diabetes;
			else do;
				
				%*T1DM;
				if dx_DiabPreGT1DMIdx30 >= 1 and
					dx_DiabPreGT2DMIdx30 < 2 and 
					rx_t2dm_antidiable90GAidx30 <= 0 and
					rx_insulinle90GAidx30 >= 1 then diabetes_type_post_dx = "T1DM";

				%*T2DM;
				else if rx_t2dm_antidiable90GAidx30 >= 2 or
							rx_insulinle90GAidx30 = 0 or
							(dx_DiabPreGT2DMIdx30 >= 2 and dx_DiabPreGT1DMIdx30 < 2 and rx_insulinle90GAidx30 >= 1) 
						then diabetes_type_post_dx = "T2DM";


				%*Not otherwise specified;
				else diabetes_type_post_dx = "NOS";
				
			end;
		
	run;

/*	proc freq data=pregnancies3_diab;*/
/*		table diabetes_type_pre / missing;*/
/*	run;*/


	%*Distinguish the type of diabetes among those with pregestational diabetes;
	data pregnancies3_diab2;
	set pregnancies3_diab;

		%*Use only pre-index information;
		if any_diabetes_pre = 0 then diabetes_type_pre = diabetes_type_pre;
			
			else if pregestation_diab_pre = 1 then diabetes_type_pre = diabetes_type_pre;

			else do;

				%*Early onset GDM on newly diagnosed pregetational OR no claims in the past months for pregestational;
				if dt_gtt91to140Idx ne . and
					dx_gtt91to140idx >= 1 and
					sum(dx_DiabPreGOutptIdx, dx_DiabPreGInptIdx) = 0 and
					rx_nonmet_antidiable90GAidx = 0 then diabetes_type_pre = "EOND";

					else if (dx_DiabGDMge141idx >= 1 or
								(pr_gttge141idx = 1 and dx_DiabDMge141idx >= 1)) and
							dx_DiabDMle140Idx = 0 and
							rx_nonmet_antidiable140GAidx = 0 then diabetes_type_pre = "GDM";

					else diabetes_type_pre = "UNSP";

			end;
			
		%*Include information up to 30 days post-index;
		if any_diabetes_post = 0 then diabetes_type_post = diabetes_type_post;
			
			else if pregestation_diab_post = 1 then diabetes_type_post = diabetes_type_post;

			else do;

				%*Early onset GDM on newly diagnosed pregetational OR no claims in the past months for pregestational;
				if dt_gtt91to140Idx30 ne . and
					dx_gtt91to140idx30 >= 1 and
					sum(dx_DiabPreGOutptIdx30, dx_DiabPreGInptIdx30) = 0 and
					rx_nonmet_antidiable90GAidx30 = 0 then diabetes_type_post = "EOND";

					else if (dx_DiabGDMge141idx30 >= 1 or
								(pr_gttge141idx30 = 1 and dx_DiabDMge141idx30 >= 1)) and
							dx_DiabDMle140Idx30 = 0 and
							rx_nonmet_antidiable140GAidx30 = 0 then diabetes_type_post = "GDM";

					else diabetes_type_post = "UNSP";

			end;
			
			
		%*****require a diagnosis code;
			
		%*Use only pre-index information;
		if any_diabetes_pre_dx = 0 then diabetes_type_pre_dx = diabetes_type_pre_dx;
			
			else if pregestation_diab_pre_dx = 1 then diabetes_type_pre_dx = diabetes_type_pre_dx;

			else do;

				%*Early onset GDM on newly diagnosed pregetational OR no claims in the past months for pregestational;
				if dt_gtt91to140Idx ne . and
					dx_gtt91to140idx >= 1 and
					sum(dx_DiabPreGOutptIdx, dx_DiabPreGInptIdx) = 0 and
					rx_nonmet_antidiable90GAidx = 0 then diabetes_type_pre_dx = "EOND";

					else if (dx_DiabGDMge141idx >= 1 or
								(pr_gttge141idx = 1 and dx_DiabDMge141idx >= 1)) and
							dx_DiabDMle140Idx = 0 and
							rx_nonmet_antidiable140GAidx = 0 then diabetes_type_pre_dx = "GDM";

					else diabetes_type_pre_dx = "UNSP";

			end;
			
		%*Include information up to 30 days post-index;
		if any_diabetes_post_dx = 0 then diabetes_type_post_dx = diabetes_type_post_dx;
			
			else if pregestation_diab_post_dx = 1 then diabetes_type_post_dx = diabetes_type_post_dx;

			else do;

				%*Early onset GDM on newly diagnosed pregetational OR no claims in the past months for pregestational;
				if dt_gtt91to140Idx30 ne . and
					dx_gtt91to140idx30 >= 1 and
					sum(dx_DiabPreGOutptIdx30, dx_DiabPreGInptIdx30) = 0 and
					rx_nonmet_antidiable90GAidx30 = 0 then diabetes_type_post_dx = "EOND";

					else if (dx_DiabGDMge141idx30 >= 1 or
								(pr_gttge141idx30 = 1 and dx_DiabDMge141idx30 >= 1)) and
							dx_DiabDMle140Idx30 = 0 and
							rx_nonmet_antidiable140GAidx30 = 0 then diabetes_type_post_dx = "GDM";

					else diabetes_type_post_dx = "UNSP";

			end;

	run;

/*	proc freq data=pregnancies3_diab2;*/
/*		table diabetes_type_pre*pregestation_diab_pre / missing;*/
/*	run;*/

/*	%*Revise the numOutptPNC and InptPNC;*/
/*	proc sql;*/
/*		create table pregnancies_counts as*/
/*		select distinct a.*, case when sum(b.pren_outpatient) = . then 0 else sum(b.pren_outpatient) end as num_OutptPNC,*/
/*				case when sum(b.pren_inpatient) = . then 0 else sum(b.pren_inpatient) end as num_InptPNC*/
/*		from pregnancies3_diab2 (drop = num_OutptPNC num_InptPNC) as a*/
/*		left join (select distinct * from temp.codeprenatal_meg1_dts) as b*/
/*		on a.patient_deid = b.patient_deid and a.dt_indexprenatal <= b.enc_date <= a.dt_index*/
/*		group by a.enrolid, a.idxpren*/
/*		;*/
/*		quit;*/

	%*Output the final dataset;
	data 
		%if &suffix = NA %then %do;
			ana.preg_covar_&lmpindex._&lookbackdt._&lookbackdays;
		%end;
		%else %do;
			ana.preg_covar_&lmpindex._&lookbackdt._&lookbackdays._&suffix;
		%end;
	set pregnancies3_diab2;

		%*****Create a final treatment variable that is a binary 1 or 0;
		if exposure = 'NIFEDIPINE' then trt = 1;
			else trt = 0;
			
		%*****Combine all substance use disorders;
		if oud_pre = 1 or alc_pre = 1 or othersud_pre = 1 or smk_pre = 1 then substance_use_pre = 1;
			else substance_use_pre = 0;
			
		%*****Combine hyper and hypothyroid disorders ;
	
		%*Using pre-index information;
		if hyperthy_pre = 1 or hypothy_pre = 1 then thyroid_disorder_pre = 1;
			else thyroid_disorder_pre = 0;
		%*Include an indicator for receiving pharmacotherapy for this thyroid disease.;
		if thyroid_disorder_pre = 1 and (hyperthyrx_pre = 1 or hypothyrx_pre = 1) then thyroid_disorder_trt_pre = 2;
			else if thyroid_disorder_pre = 1 then thyroid_disorder_trt_pre = 1;
			else thyroid_disorder_trt_pre = 0;
		
		%*Using post-index information;
		if hyperthy_post = 1 or hypothy_post = 1 then thyroid_disorder_post = 1;
			else thyroid_disorder_post = 0;
		%*Include an indicator for receiving pharmacotherapy for this thyroid disease.;
		if thyroid_disorder_post = 1 and (hyperthyrx_post = 1 or hypothyrx_post = 1) then thyroid_disorder_trt_post = 2;
			else if thyroid_disorder_post = 1 then thyroid_disorder_trt_post = 1;
			else thyroid_disorder_trt_post = 0;
			
		%********Combine medications for hyper and hypothyroid disorders for its own separate indicator;
		if hyperthyrx_pre = 1 or hypothyrx_pre = 1 then thyroidrx_pre = 2;
			else thyroidrx_pre = 0;
		if hyperthyrx_post = 1 or hypothyrx_post = 1 then thyroidrx_post = 1;
			else thyroidrx_post = 0;
			
		%********Indicators for specific gestational age codes;
		
		%*prior to index date;
		specific_ga_preidx = num_specificga_preidx >= 1;
		%*throughout the entire pregnancy;
		specific_ga_preg = num_specificga_preg >= 1;
		
		%*******Bipolar disorder - incorporate medication treatment;
		
		%*Only information up to index;
		if bipolar_pre = 1 and MoodStabrx_pre = 1 then bipolar_trt_pre = 2;
			else if bipolar_pre = 1 then bipolar_trt_pre = 1;
			else bipolar_trt_pre = 0;
			
		%*Include information up to 30 days post-index;
		if bipolar_post = 1 and MoodStabrx_post = 1 then bipolar_trt_post = 2;
			else if bipolar_post = 1 then bipolar_trt_post = 1;
			else bipolar_trt_post = 0;
			
		%******Schizophrenia or schizoaffective disorder - incorporate medication treatment;
		
		%*Include only information up to index;
		if schizo_pre = 1 and antipsyrx_pre = 1 then schizo_trt_pre = 2;
			else if schizo_pre = 1 then schizo_trt_pre = 1;
			else schizo_trt_pre = 0;
			
		%*Include information up to 30 days post-index;
		if schizo_post = 1 and antipsyrx_post = 1 then schizo_trt_post = 2;
			else if schizo_post = 1 then schizo_trt_post = 1;
			else schizo_trt_post = 0;
		
	run;

%mend;



%*Run the analysis;
%identify_variables(LOOKBACKDT=dt_index, LOOKBACKDAYS=270, LMPINDEX=63, SUFFIX=NA);
%identify_variables(LOOKBACKDT=dt_lmp, LOOKBACKDAYS=180, LMPINDEX=63, SUFFIX=NA);
%identify_variables(LOOKBACKDT=dt_index, LOOKBACKDAYS=270, LMPINDEX=42, SUFFIX=NA);
%identify_variables(LOOKBACKDT=dt_index, LOOKBACKDAYS=270, LMPINDEX=84, SUFFIX=NA);
%identify_variables(LOOKBACKDT=dt_lmp, LOOKBACKDAYS=180, LMPINDEX=63, SUFFIX=DEL);

%*Applying the 3-week window around the LMP;
%identify_variables(LOOKBACKDT=dt_index, LOOKBACKDAYS=270, LMPINDEX=63, SUFFIX=m21);
%identify_variables(LOOKBACKDT=dt_index, LOOKBACKDAYS=270, LMPINDEX=63, SUFFIX=p21);



















