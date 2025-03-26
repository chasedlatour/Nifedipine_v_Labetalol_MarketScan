/********************************************************************************************************************************************
PROGRAM: 04_primary_analyses.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to clean the pregnancy cohort that we derived from the MarketScan claims data
and prepare it for analyses.
	
Goal: 
Output data: 

Date: 12.19.2024
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
%setup(sample=full, programname=04_primary_analyses, savelog=N);

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

										01 - IMPLEMENT EXCLUSION CRITERIA AND GET RELEVANT COUNTS

********************************************************************************************************************************************/

%*First create working directory copy of pregnancies.;

data pregnancies (where = (ga_index_days < 161)); %*Index had to occur prior to 23w of gestation;
set ana.preg_covar_63_dt_index_270;

	%*Derive the gestational age at the index date;
	ga_index_days = dt_index - dt_lmp;

run;

%*Output count;
proc sql noprint;
	select count(distinct idxpren) into :num_preg from pregnancies;
	quit;
%put Number of new users of nifedipine or labetalol: &num_preg;



%*Restrict to those pregnancies that had >= 270 days of continuous enrollment in a plan with RX coverage
prior to the index date;
data pregnancies_enr;
set pregnancies;
	where prior_enrollment_rx >= 270;
run;

%*Output counts;
proc sql noprint;
	select count(distinct idxpren) into :num_preg_enr from pregnancies_enr;
	quit;
%put Number of pregnancies with qualifying enrollment: &num_preg_enr ;
%put Number of pregnancies without qualifying enrolllment: %eval(&num_preg - &num_preg_enr );



%*Now get counts for the rest of the exclusion criteria;
proc sql noprint;
	select sum(age_at_index < 18), sum(asthma_pre = 1), 
		sum(coronary_heart_disease_pre = 1 or arrhythmia_pre = 1 or congenital_heart_pre = 1 or endocarditis_pre = 1 or
				myopericarditis_pre = 1 or heartfailure_pre = 1 or heart_valve_disease_pre = 1 or cardiomyopathy_pre = 1 or
				other_heart_disease_pre = 1), sum(cancer_pre = 1), sum(preg_outcome_clean = 'EM'),
				sum(stroke_pre = 1 or mi_pre = 1 or angina_pre = 1 or athero_pre = 1 or perivasc_pre = 1),
				sum(retino_pre = 1), sum(antiphos_pre = 1), sum(lupus_pre = 1), sum(sickled_pre = 1)
	into :num_preg_less18, :num_preg_asthma, :num_heart_disease, :num_cancer, :num_em, :num_cvd, :num_retino,
			:num_antiphos, :num_lupus, :num_sickled
	from pregnancies_enr;
	quit;
%put Number of pregnancies where person was <18 years of age: &num_preg_less18;
%put Number of pregnancies where person had asthma: &num_preg_asthma;
%put Number of pregnancies where person had heart disease: &num_heart_disease;
%put Number of pregnancies where person had cancer: &num_cancer;
%put Number of pregnancies where outcome was ectopic or molar: &num_em;
%put Number of pregnancies where person had cardiovascular disease: &num_cvd;
%put Number of pregnancies where person had diabetic retinopathy: &num_retino;
%put Number of pregnancies where person had antiphospholipid syndrome: &num_antiphos;
%put Number of pregnancies where person had lupus: &num_lupus;
%put Number of pregnancies where person had sickle cell disease: &num_sickled;


%*Remove individuals that meet any of those criteria;
data pregnancies_excl;
set pregnancies_enr;
	if age_at_index < 18 or 
			sum(asthma_pre, coronary_heart_disease_pre, arrhythmia_pre, congenital_heart_pre, endocarditis_pre,
			myopericarditis_pre, heartfailure_pre, heart_valve_disease_pre, cardiomyopathy_pre, other_heart_disease_pre, cancer_pre,
			stroke_pre, mi_pre, angina_pre, athero_pre, perivasc_pre, retino_pre, antiphos_pre, lupus_pre,
			sickled_pre = 1) >0 or
			preg_outcome_clean = 'EM' then delete;
run;

%*Output counts;
proc sql noprint;
	select count(distinct idxpren), count(distinct enrolid)
	into :num_preg_excl, :num_ppl
	from pregnancies_excl;
	quit;

%put Number of pregnancies after implementing exclusion criteria: &num_preg_excl;
%put Number of pregnancies removed after retaining only first: %eval(&num_preg_excl - &num_ppl);



%*Now, actually remove those that are after the first pregnancy;
proc sort data=pregnancies_excl;
	by enrolid dt_lmp;
run;
%*output the first pregnancy - this is the primary analysis cohort;
data ana.primary_cohort;
format year_index $4.;
set pregnancies_excl;
	by enrolid dt_lmp;
	
	%*Create gestational age categorical variable that aligns with CHAP study;
	if . < ga_index_days < 98 then ga_index_lt14 = 1;
		else ga_index_lt14 = 0;
		
	%*Create numeric gestational age indicator;
	if . < ga_index_days < 98 then ga_index_cat = 1;
		else ga_index_cat = 2;
		
	numOutptPNC_lt3 = num_OutptPNC < 3;
	numInptADM_gt1 = num_InptAdm > 1;
	
	%*Make a simplified categorical diabetes variable;
	if diabetes_type_post in ("T1DM" "T2DM") then diabetes_simp = "PRE";
		else if diabetes_type_post = "NA" then diabetes_simp = "NA";
		else diabetes_simp = "OTH";
		

		
	%*Output only the first pregnancy to a person;
	if first.enrolid then output;

run;

*Get quartiles of gestational age within ga strata;
proc sort data=ana.primary_cohort; by ga_index_cat; run;
proc means data=ana.primary_cohort p25 median p75;
	by ga_index_cat;
	var ga_index_days;
run;

*Conduct some final data cleaning for input into the PS and censoring models;
data ana.primary_cohort;
set ana.primary_cohort;
	format preg_outcome_clean $pregoutc.
			thyroid_disorder_trt_post treated.
			bipolar_trt_post treated.
			schizo_trt_post treated.
			rural unknown.
			diabetes_type_post $diabetes.
			diabetes_simp $simpdiab.
			dt_gapreg_ptb MMDDYY10.;
			
	label ga_index_days = "Gestational age (days) at index date"
		ga_index_lt14 = "Gestational age <14wks at index date"
		age_at_index = "Maternal age at index"
		year_index = "Year of index date"
		year_le2019 = "Index prior to 2019"
		preg_outcome_clean = "Pregnancy outcome"
		substance_use_pre = "Substance Use"
		nausea_pre = "Nausea/Vomiting"
		any_diabetes_post = "Any diabetes"
		pregestation_diab_post = "Pregestational diabetes"
		diabetes_type_post = "Type of diabetes"
		t2dmrx_post = "Type-2 diabetes medications"
		t1t2dmrx_post = "Type-1 or 2 diabetes medications"
		metforrx_post = "Metformin"
		obesity_post = "Obesity diagnosis"
		glp1wgtrx_post = "GLP-1 weightloss medications"
		otherwgtrx_post = "Other weightloss medications"
		bariatric_post = "Bariatric surgery"
		migraine_pre = "Migraine"
		recurlos_pre = "Recurrent loss"
		ckd_post = "Chronic kidney disease"
		thyroid_disorder_post = "Thyroid disease diagnosis"
		thyroidrx_post = "Thyroid disease medication fill"
		thyroid_disorder_trt_post = "Thyroid disease diagnosis and medication"
		depressi_post = "Depression"
		anxiety_post = "Anxiety"
		antideprx_post = "Antidepressants"
		adhd_post = "ADHD"
		adhdrx_post = "ADHD medications"
		bipolar_post = "Bipolar disorder diagnosis"
		moodstabrx_post = "Mood stabilizer fills"
		bipolar_trt_post = "Bipolar disorder diagnosis and medication"
		ptsd_post = "PTSD"
		schizo_post = "Schizophrenia or schizoaffective disorder"
		antipsyrx_post = "Antipsychotic fills"
		schizo_trt_post = "Schizophrenia diagnosis and medication"
		rural = "Rurality"
		hyperlip_post = "Hyperlipidemia"
		teratrx_pre = "Teratogenic medications"
		benzorx_post = "Benzodiazepines"
		anticonvulrx_post = "Anticonvulsants"
		num_OutptPNC = "Number of outpatient prenatal encounters prior to index"
		numOutptPNC_lt3 = "Less than 1 or 2 prenatal encounters prior to index"
		num_InptAdm = "Number of inpatient admissions prior to index"
		numInptADM_gt1 = "More than 1 inpatient admission prior to index"
		diabetes_simp = "Type of diabetes"
		chronichypertension_pre = "Chronic hypertension"
		preeclampsia_pre = "Preeclampsia"
		;

	*First, just combine 2016 and 2017;
	if year_index in ("2016", "2017") then year_index2017 = "2017";
		else year_index2017 = year_index;
		
	*Slightly more flexible year variable;
	if year_index in ("2016", "2017") then year_index4 = "2017";
		else if year_index in ("2018", "2019") then year_index4 = "2019";
		else if year_index in ("2020") then year_index4 = "2020";
		else if year_index in ("2021", "2022") then year_index4 = "2022";
		
	*Now, combine pre 2019, 2019-2020, and 2021-2022;
	if year_index in ("2016", "2017", "2018") then year_index2021 = "2018";
		else if year_index in ("2019", "2020") then year_index2021 = "2020";
		else year_index2021 = "2022";
		
	*Create a squared version of the age, ga_index, and num_outptpnc variables;
	age_at_index_2 = age_at_index**2;
	ga_index_days_2 = ga_index_days**2;
	num_outptpnc_2 = num_outptpnc**2;
	
	*GA quartiles - quartile values were determined within ga_index_cat strata;
	if . < ga_index_days <= 56 then ga_quartile = 1;
		else if ga_index_days <= 72 then ga_quartile = 2;
		else if ga_index_days <= 86 then ga_quartile = 3;
		else if ga_index_days < 98 then ga_quartile = 4;
		else if ga_index_days <= 112 then ga_quartile = 1;
		else if ga_index_days <= 126 then ga_quartile = 2;
		else if ga_index_days <= 144 then ga_quartile = 3;
		else ga_quartile = 4;
		
	*Assume that missing MSA is urban -- Aligns with probability assumptions;
	if rural = 1 then rural2 = 1;
		else rural2 = 0;

	/*For preterm birth analyses, we want to know individuals experienced a live birth prior to 37 weeks of gestation. If
	they do not, then they receive an outcome of TB (term birth) that occurs at the end of the risk period (36w6d).
	We need to then reassign their event indicator and event time.*/

	*Create outcome variable for preterm birth analyses;
	if preg_outcome_clean in ('LBM' 'LBS' 'UDL') and dt_gapreg - dt_lmp < 259 then preg_outcome_ptb = 'PTB';
		else if dt_gapreg - dt_lmp >= 259 then preg_outcome_ptb = 'TB'; /*If the outcome occurs after 37 weeks, then they survived the risk period. */
		else preg_outcome_ptb = preg_outcome_clean;

	if preg_outcome_ptb = 'TB' then dt_gapreg_ptb = dt_lmp + 258; /*36*7+6 = 258 -- They have a non-event at the end of the risk period*/
		else dt_gapreg_ptb = dt_gapreg;
	
run;

*Create a categorical variable for maternal age at index;
data ana.primary_cohort;
set ana.primary_cohort;
	if  age_at_index <= 24 then age_at_index_cat = "<= 24";
		else if age_at_index <= 29 then age_at_index_cat = "25-29";
		else if age_at_index <= 34 then age_at_index_cat = "30-34";
		else if age_at_index <= 39 then age_at_index_cat = "35-39";
		else age_at_index_cat = ">= 40";
run;











/********************************************************************************************************************************************

												02 - INVESTIGATE COVARIATE DISTRIBUTIONS

********************************************************************************************************************************************/



*First, describe the cohort with all of the variables that we had proposed including;

%*Now do some table 1 descriptives of the cohort - ana.primary_cohort;
%table1(inds = ana.primary_cohort, colVar = exposure,
	rowVars = ga_index_days ga_index_lt14 age_at_index age_at_index_cat year_index year_le2019 chronichypertension_pre preeclampsia_pre
		substance_use_pre nausea_pre pregestation_diab_post t2dmrx_post t1t2dmrx_post metforrx_post
		obesity_post glp1wgtrx_post otherwgtrx_post bariatric_post migraine_pre recurlos_pre ckd_post 
		thyroid_disorder_post thyroidrx_post 
		depressi_post anxiety_post antideprx_post adhd_post adhdrx_post bipolar_post moodstabrx_post
		ptsd_post schizo_post antipsyrx_post rural hyperlip_post teratrx_pre benzorx_post anticonvulrx_post
		num_OutptPNC numOutptPNC_lt3 num_InptAdm  numInptADM_gt1,
	outfile = Table 1: Primary cohort overall without weights,
	title = Table 1: Primary cohort);
	

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
		substance_use_pre nausea_pre pregestation_diab_post t2dmrx_post t1t2dmrx_post metforrx_post
		obesity_post glp1wgtrx_post otherwgtrx_post bariatric_post migraine_pre recurlos_pre ckd_post 
		thyroid_disorder_post thyroidrx_post 
		depressi_post anxiety_post antideprx_post adhd_post adhdrx_post bipolar_post moodstabrx_post
		ptsd_post schizo_post antipsyrx_post rural hyperlip_post teratrx_pre benzorx_post anticonvulrx_post
		num_OutptPNC numOutptPNC_lt3 num_InptAdm  numInptADM_gt1,
	outfile = Table 1: Unweighted primary cohort GA lt 14wk,
	title = Table 1: Primary cohort where GA less than 14w at index);
	
%*Now do some table 1 descriptives of the cohort - ana.primary_cohort;
%table1(inds = primary_ge14, colVar = exposure,
	rowVars = ga_index_days age_at_index age_at_index_cat year_index year_le2019 chronichypertension_pre preeclampsia_pre
		substance_use_pre nausea_pre pregestation_diab_post t2dmrx_post t1t2dmrx_post metforrx_post
		obesity_post glp1wgtrx_post otherwgtrx_post bariatric_post migraine_pre recurlos_pre ckd_post 
		thyroid_disorder_post thyroidrx_post 
		depressi_post anxiety_post antideprx_post adhd_post adhdrx_post bipolar_post moodstabrx_post
		ptsd_post schizo_post antipsyrx_post rural hyperlip_post teratrx_pre benzorx_post anticonvulrx_post
		num_OutptPNC numOutptPNC_lt3 num_InptAdm  numInptADM_gt1,
	outfile = Table 1: Unweighted primary cohort GA ge 14wk,
	title = Table 1: Primary cohort where GA less than 14w at index);
	
	
	
*Fit PS model within each strata. Include variables for inverse probability of censoring weights;

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
			year_index2017*t2dmrx_post diabetes_simp*t2dmrx_post
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
	outds=test, /*ana.primary_point,*/
	outds_dist=NA, /*ga_dist_primary_point,*/
	outds_ps = ps_primary_point,
	outds_surv = NA /*ana.surv_primary_point*/
	);


*Put out the trimmed cohorts for at-risk counts in cumulative incidence figure;

data ana.primary_cohort_trim_1;
set _anacohort_trim1;
run;
title 'Population distribution by treatment after trimming: <14gw';
proc freq data=ana.primary_cohort_trim_1;
	table trt;
run;
data ana.primary_cohort_trim_2;
set _anacohort_trim2;
run;
title 'Population distribution by treatment after trimming: >=14gw';
proc freq data=ana.primary_cohort_trim_2;
	table trt;
run;
title ;
	
*********Look at the PS distribution by GA strata;

*Less than 14 weeks;
proc sort data=ana.ps_primary_point_num_1; by trt; run;
ods graphics on / reset imagename="PS primary cohort GA lt 14 weeks";
proc sgplot data=ana.ps_primary_point_num_1;
	title "Propensity scores by treatment: Less than 14 weeks at index";
	histogram d / group=trt transparency=0.2;
/* 	density d / group = trt; */
run;

*Greater than or equal to 14 weeks;
proc sort data=ana.ps_primary_point_num_2; by trt; run;
ods graphics on / reset imagename="PS primary cohort GA ge 14 weeks";
proc sgplot data=ana.ps_primary_point_num_2;
	title "Propensity scores by treatment: At least 14 weeks at index";
	histogram d / group=trt transparency=0.2;
/* 	density d / group = trt; */
run;

proc means data=ana.ps_primary_point_num_1 min mean max;
	title "Gestational age <14 weeks";
	class trt;
	var d;
run;
proc means data=ana.ps_primary_point_num_2 min mean max;
	title "Gestational age >=14 weeks";
	class trt;
	var d;
run;


*Attach these weights to peoples observations for a weighted Table 1 by strata;
proc sql;
	create table weighted as 
	select a.*, case when missing(c.expwgt) then 1 else 2 end as ga_strat,
			coalesce(b.expwgt, c.expwgt) as smrw
	from ana.primary_cohort as a
	left join ana.ps_primary_point_num_1 as b
	on a.idxpren=b.idxpren
	left join ana.ps_primary_point_num_2 as c
	on a.idxpren=c.idxpren
	/*Subset to those that remained in the trimmed sample*/
	having idxpren in (select distinct idxpren 
						from ana.ps_primary_point_num_1 
						union
						select distinct idxpren
						from ana.ps_primary_point_num_2)
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
		substance_use_pre nausea_pre pregestation_diab_post t2dmrx_post t1t2dmrx_post metforrx_post
		obesity_post glp1wgtrx_post otherwgtrx_post bariatric_post migraine_pre recurlos_pre ckd_post 
		thyroid_disorder_post thyroidrx_post 
		depressi_post anxiety_post antideprx_post adhd_post adhdrx_post bipolar_post moodstabrx_post
		ptsd_post schizo_post antipsyrx_post rural hyperlip_post teratrx_pre benzorx_post anticonvulrx_post
		num_OutptPNC numOutptPNC_lt3 num_InptAdm  numInptADM_gt1,
	outfile = Table 1: Weighted primary cohort lt 14w,
	title = Table 1: Primary cohort where GA less than 14w at index);
	
%*Now do some table 1 descriptives of the cohort - ana.primary_cohort;
%table1(inds = primary_ge14, colVar = exposure, wgtvar=smrw,
	rowVars = ga_index_days ga_index_lt14 age_at_index age_at_index_cat year_index year_le2019 chronichypertension_pre preeclampsia_pre
		substance_use_pre nausea_pre pregestation_diab_post t2dmrx_post t1t2dmrx_post metforrx_post
		obesity_post glp1wgtrx_post otherwgtrx_post bariatric_post migraine_pre recurlos_pre ckd_post 
		thyroid_disorder_post thyroidrx_post 
		depressi_post anxiety_post antideprx_post adhd_post adhdrx_post bipolar_post moodstabrx_post
		ptsd_post schizo_post antipsyrx_post rural hyperlip_post teratrx_pre benzorx_post anticonvulrx_post
		num_OutptPNC numOutptPNC_lt3 num_InptAdm  numInptADM_gt1,
	outfile = Table 1: Weighted primary cohort ge 14w,
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
	outds=ana.primary_point,
	outds_dist=ga_dist_primary_point,
	outds_ps = ps_primary_point,
	outds_surv = ana.surv_primary_point
	);

*Get the sample sizes across the two treatment groups;
proc freq data=ana.ps_primary_point_num_1;
	table trt / missing;
run;
proc freq data=ana.ps_primary_point_num_2;
	table trt / missing;
run;

*Look at the stratified estimates

Datasets:
- ana.primary_point_noipcw_1
- ana.primary_point_noipcw_2;

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_primary_point,
	numgastrat = 2,
	est = ana.primary_point_,
	outds = ana.primary_point_overall
	);

	
*******Output the cumulative incidence curves, stratified by treatment and Gestational age at index;
data surv_primary1;
set ana.surv_primary_point_1;
	cum_comprisk1_2 = sum(cum_outcome, cum_comprisk1);
	cum_comprisk2_2 = sum(cum_outcome, cum_comprisk1, cum_comprisk2);
run;

*GA <14 weeks at index and untreated;
ods graphics on / reset imagename="Untreated combined cumulative Incidence Pregnancy Loss GA lt 14 weeks";
proc sgplot data=surv_primary1 (where = (trt = 0));
    series x=days y=cum_outcome / lineattrs=(color=red thickness=2) legendlabel="Pregnancy Loss";
    series x=days y=cum_comprisk1_2 / lineattrs=(color=blue thickness=2) legendlabel="Induced Abortion";
    series x=days y=cum_comprisk2_2 / lineattrs=(color=green thickness=2) legendlabel="Live Birth";
    xaxis label="Days from Index Fill";
    yaxis label="Cumulative Incidence";
    keylegend / position=topright;
    title 'Combined Cumulative Incidence of Each Outcome Among the Untreated where GA at Index was < 14 Weeks';
run;
*GA <14 weeks at index and treated;
ods graphics on / reset imagename="Treated combined cumulative Incidence Pregnancy Loss GA lt 14 weeks";
proc sgplot data=surv_primary1 (where = (trt = 1));
    series x=days y=cum_outcome / lineattrs=(color=red thickness=2) legendlabel="Pregnancy Loss";
    series x=days y=cum_comprisk1_2 / lineattrs=(color=blue thickness=2) legendlabel="Induced Abortion";
    series x=days y=cum_comprisk2_2 / lineattrs=(color=green thickness=2) legendlabel="Live Birth";
    xaxis label="Days from Index Fill";
    yaxis label="Cumulative Incidence";
    keylegend / position=topright;
    title 'Combined Cumulative Incidence of Each Outcome Among the Treated where GA at Index was < 14 Weeks';
run;

data surv_primary2;
set ana.surv_primary_point_2;
	cum_comprisk1_2 = sum(cum_outcome, cum_comprisk1);
	cum_comprisk2_2 = sum(cum_outcome, cum_comprisk1, cum_comprisk2);
run;
*GA >=14 weeks at index and untreated;
ods graphics on / reset imagename="Untreated combined cumulative Incidence Pregnancy Loss GA ge 14 weeks";
proc sgplot data=surv_primary2 (where = (trt = 0));
    series x=days y=cum_outcome / lineattrs=(color=red thickness=2) legendlabel="Pregnancy Loss";
    series x=days y=cum_comprisk1_2 / lineattrs=(color=blue thickness=2) legendlabel="Induced Abortion";
    series x=days y=cum_comprisk2_2 / lineattrs=(color=green thickness=2) legendlabel="Live Birth";
    xaxis label="Days from Index Fill";
    yaxis label="Cumulative Incidence";
    keylegend / position=topright;
    title 'Combined Cumulative Incidence of Each Outcome Among the Untreated where GA at Index was >= 14 Weeks';
run;
*GA >=14 weeks at index and treated;
ods graphics on / reset imagename="Treated combined cumulative Incidence Pregnancy Loss GA lt 14 weeks";
proc sgplot data=surv_primary2 (where = (trt = 1));
    series x=days y=cum_outcome / lineattrs=(color=red thickness=2) legendlabel="Pregnancy Loss";
    series x=days y=cum_comprisk1_2 / lineattrs=(color=blue thickness=2) legendlabel="Induced Abortion";
    series x=days y=cum_comprisk2_2 / lineattrs=(color=green thickness=2) legendlabel="Live Birth";
    xaxis label="Days from Index Fill";
    yaxis label="Cumulative Incidence";
    keylegend / position=topright;
    title 'Combined Cumulative Incidence of Each Outcome Among the Treated where GA at Index was >= 14 Weeks';
run;



*******Output the cumulative incidence curves for the study outcome only, stratified by Gestational age at index;
*GA <14 weeks at index;
ods graphics on / reset imagename="Cumulative Incidence Pregnancy Loss GA lt 14 weeks";
proc sgplot data=surv_primary1;
    series x=days y=cum_outcome / group=trt;
    xaxis label="Days from Index Fill";
    yaxis label="Cumulative Incidence";
    keylegend / position=topright;
    title 'Combined Cumulative Incidence of Pregnancy Loss Among those with GA at Index < 14 Weeks';
run;
*GA >=14 weeks at index;
ods graphics on / reset imagename="Cumulative Incidence Pregnancy Loss GA ge 14 weeks";
proc sgplot data=surv_primary2;
    series x=days y=cum_outcome / group=trt;
    xaxis label="Days from Index Fill";
    yaxis label="Cumulative Incidence" min=0 max=0.12;
    keylegend / position=topright;
    title 'Cumulative Incidence Pregnancy Loss Among those with GA at Index >= 14 Weeks';
run;


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
	outds=ana.primary_boot,
	outds_dist=ga_dist_primary_boot,
	outds_ps = ps_primary_boot,
	outds_surv = NA
	);
	
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.primary_boot, 
			inputDist= ana.ga_dist_primary_boot,
			numStrata= 2, 
			output_stratified= ana.primary_boot_strat,
			output_overall= ana.primary_boot_overall);


*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.primary_boot_1, pointdsn=ana.primary_point_1, output=ana.primary_boot_1_ci);
%strat_estimates_w_CI(bootdsn=ana.primary_boot_2, pointdsn=ana.primary_point_2, output=ana.primary_boot_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.primary_boot_overall, pointdsn=ana.primary_point_overall, output=ana.primary_overall_boot_ci);
	
	
%count_missing_zero(inds1=ana.primary_boot_1, inds2=ana.primary_boot_2);
















	
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
	outds=ana.primary_point_ptb,
	outds_dist=ga_dist_primary_point_ptb,
	outds_ps = ps_primary_point_ptb,
	outds_surv = ana.surv_primary_point_ptb
	);

*Get the sample sizes across the two treatment groups;
proc freq data=ana.ps_primary_point_num_1;
	table trt / missing;
run;
proc freq data=ana.ps_primary_point_num_2;
	table trt / missing;
run;

*Look at the stratified estimates;

*Combine the stratified estimates according to the distribution of GA at index among
those in the treated group;
%combine_point_estimates(
	gadist = ana.ga_dist_primary_point_ptb,
	numgastrat = 2,
	est = ana.primary_point_ptb_,
	outds = ana.primary_point_ptb_overall
	);

	

/********Output the cumulative incidence curves for the study outcome only, stratified by Gestational age at index;*/
*GA <14 weeks at index;
ods graphics on / reset imagename="Cumulative Incidence of Preterm Birth_GA lt14 Wk at Index";
proc sgplot data=ana.surv_primary_point_ptb_1;
    series x=days y=cum_outcome / group=trt;
    xaxis label="Days from Index Fill";
    yaxis label="Cumulative Incidence" min=0 max=0.35;
    keylegend / position=topright;
    title 'Cumulative Incidence of Preterm Birth Among those with GA at Index < 14 Weeks';
run;
*GA >=14 weeks at index;
ods graphics on / reset imagename="Cumulative Incidence of Preterm Birth_GA ge14 Wk at Index";
proc sgplot data=ana.surv_primary_point_ptb_2;
    series x=days y=cum_outcome / group=trt;
    xaxis label="Days from Index Fill";
    yaxis label="Cumulative Incidence" min=0 max=0.35;
    keylegend / position=topright;
    title 'Cumulative Incidence of Preterm Birth Among those with GA at Index >= 14 Weeks';
run;


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
	outds=ana.primary_boot_ptb,
	outds_dist=ga_dist_primary_boot_ptb,
	outds_ps = ps_primary_boot_ptb,
	outds_surv = NA
	);
	
	
options mlogic mprint symbolgen notes;
*Combine all the bootstrapped estimates;
%combine_boot_estimates(
			inputEst= ana.primary_boot_ptb,
			inputDist= ana.ga_dist_primary_boot_ptb,
			numStrata= 2, 
			output_stratified= ana.primary_ptb_boot_strat,
			output_overall= ana.primary_ptb_boot_overall);

*Stratified estimates with confidence intervals.;
%strat_estimates_w_CI(bootdsn=ana.primary_boot_ptb_1, pointdsn=ana.primary_point_ptb_1, output=ana.primary_boot_ptb_1_ci);
%strat_estimates_w_CI(bootdsn=ana.primary_boot_ptb_2, pointdsn=ana.primary_point_ptb_2, output=ana.primary_boot_ptb_2_ci);

*OVerall estimates with confidence interval;
%overall_estimates_w_CI(stderrdsn=ana.primary_ptb_boot_overall, pointdsn=ana.primary_point_ptb_overall, output=ana.primary_overall_ptb_boot_ci);

	
%count_missing_zero(inds1=ana.primary_boot_ptb_1, inds2=ana.primary_boot_ptb_2);













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



******
	Explore how many preterm births have ptb codes;

data ptb;
set ana.primary_cohort;
run;

*Merge the gestational age codes onto the pregnancies;
proc sql;
	create table ptb_ga as
	select a.*, b.code, b.code_type, b.description, b.preg_outcome, b.code_hierarchy, b.gestational_age_days
	from ptb as a
	left join temp.gestage as b
	on a.patient_deid=b.patient_deid and a.dt_gapreg-7 <= b.enc_date <= a.dt_gapreg+7
	;
	quit;

*Assign a hierarchy based on our algorithms hierarchy;

proc sql;
	create table ptb_ga2 as
	select distinct enrolid, idxpren, dt_lmp, dt_gapreg, dt_indexprenatal, preg_outcome_ptb, 
			min(case when code_hierarchy = "Specific gestational age" then 1
					when code_hierarchy = "Extreme prematurity" then 2
					when code_hierarchy = "Other preterm" then 3
					else 9 end) as hierarchy,
			max(case when code_hierarchy in ("Extreme prematurity", "Other preterm") then 1 else 0 end) as preterm
	from ptb_ga
	group by idxpren
	;
	quit;

proc freq data=ptb_ga2 (where = ( preg_outcome_ptb = "PTB"));
	table hierarchy preterm / missing;
run;

*Look at the gestational ages represented by specific GA codes for preterm births;
data ptb_ga3;
set ptb_ga;
	where preg_outcome_ptb = "PTB" and code_hierarchy = "Specific gestational age";

	gest_age_weeks = gestational_age_days / 7;
run;


proc means data=ptb_ga3 min p10 p25 p50 p75 p90 max;
	var gest_age_weeks;
run;
	

