/*
MACRO: implement_exclusions
PROGRAMMER: Chase Latour

PURPOSE: To implement the primary cohort exclusion criteria in another subset.

INPUT:
- inds = input dataset in which you want to implement the exclusion criteria
- outds = name of the output dataset.

*/



%macro implement_exclusions(inds=, outds=);

	%*First create working directory copy of pregnancies.;
	data pregnancies (where = (ga_index_days < 161)); %*Index had to occur prior to 23w of gestation;
	set &inds;

		%*Derive the gestational age at the index date;
		ga_index_days = dt_index - dt_lmp;

	run;

	%*Restrict to those pregnancies that had >= 270 days of continuous enrollment in a plan with RX coverage
	prior to the index date;
	data pregnancies_enr;
	set pregnancies;
		where prior_enrollment_rx >= 270;
	run;


	%*Remove individuals that meet any of the exclusion criteria;
	data pregnancies_excl;
	set pregnancies_enr;
		if age_at_index < 18 or 
				sum(asthma_pre, coronary_heart_disease_pre, arrhythmia_pre, congenital_heart_pre, endocarditis_pre,
				myopericarditis_pre, heartfailure_pre, heart_valve_disease_pre, cardiomyopathy_pre, other_heart_disease_pre, cancer_pre,
				stroke_pre, mi_pre, angina_pre, athero_pre, perivasc_pre, retino_pre, antiphos_pre, lupus_pre,
				sickled_pre = 1) >0 or
				preg_outcome_clean = 'EM' then delete;
	run;

	%*Now, remove those that are after the first pregnancy;
	proc sort data=pregnancies_excl;
		by enrolid dt_lmp;
	run;
	%*output the first pregnancy - this is the primary analysis cohort;
	data first_preg;
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

	*Conduct some final data cleaning for input into the PS and censoring models;
	data clean1;
	set first_preg;
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
	data &outds;
	set clean1;
		if  age_at_index <= 24 then age_at_index_cat = "<= 24";
			else if age_at_index <= 29 then age_at_index_cat = "25-29";
			else if age_at_index <= 34 then age_at_index_cat = "30-34";
			else if age_at_index <= 39 then age_at_index_cat = "35-39";
			else age_at_index_cat = ">= 40";
	run;


%mend;
