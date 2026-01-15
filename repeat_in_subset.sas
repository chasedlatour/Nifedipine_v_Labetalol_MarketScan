/*
MACRO: repeat_in_subset
PROGRAMMER: Chase Latour
PURPOSE: To re-run the steps of the primary analysis among some prespecified subset of the population.
*/

%macro repeat_in_subset(
		inds=ana.primary_cohort,
		sens=sens6w,
		where=NA,
		gacatvar = ga_index_cat,
		psvars=ga_quartile age_at_index age_at_index_2
				year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
				nausea_pre recurlos_pre obesity_post chronichypertension_pre
				depressi_post anxiety_post antideprx_post benzorx_post 
				teratrx_pre num_outptpnc num_outptpnc_2,
		dovars= ga_quartile age_at_index year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
				nausea_pre recurlos_pre obesity_post chronichypertension_pre depressi_post anxiety_post antideprx_post 
				benzorx_post teratrx_pre num_outptpnc num_outptpnc_2 rural2,
		dovarsmodel =  ga_quartile age_at_index year_index4
				t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
				nausea_pre recurlos_pre obesity_post chronichypertension_pre 
				depressi_post anxiety_post antideprx_post benzorx_post
				teratrx_pre num_outptpnc num_outptpnc_2 rural2,
		psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
				nausea_pre recurlos_pre obesity_post chronichypertension_pre  
				depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre, 
		doclassvars= ga_quartile year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
				nausea_pre recurlos_pre obesity_post chronichypertension_pre
				depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre rural2,
		trtvar=trt,
		numboot=2
	);

	%*****01 - IMPLEMENT EXCLUSION CRITERIA AND GET RELEVANT COUNTS;
	data pregnancies (where = (ga_index_days < 161)); %*Index had to occur prior to 23w of gestation;
	set &inds;
		%*Derive the gestational age at the index date;
		ga_index_days = dt_index - dt_lmp;
	run;
	%*Output count;
	proc sql noprint; select count(distinct idxpren) into :num_preg from pregnancies; quit;
	%put Number of new users of nifedipine or labetalol: &num_preg;

	%*Restrict to those pregnancies that had >= 270 days of continuous enrollment in a plan with RX coverage
	prior to the index date;
	data pregnancies_enr;
	set pregnancies;
		where prior_enrollment_rx >= 270;
	run;
	%*Output counts;
	proc sql noprint; select count(distinct idxpren) into :num_preg_enr from pregnancies_enr; quit;
	%put Number of pregnancies with qualifying enrollment: &num_preg_enr ;
	%put Number of pregnancies without qualifying enrolllment: %eval(&num_preg - &num_preg_enr );

	%*CDL: ADDED 11.16.2025
	Exclude those with observed outcomes that occurred on the index date;
	data pregnancies_outc;
	set pregnancies_enr;
		if preg_outcome_clean = "UNK" and dt_gapreg = dt_index then delete;
	run;
	%*Get count;
	proc sql notprint; select count(preg_outcome_clean = "UNK" and dt_gapreg = dt_index) into :num_outc from pregnancies_outc; quit;
	%put Number of pregnancies with observed outcome on the index date: &num_outc;

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
	set pregnancies_outc;
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

	%*Now, remove those that are after the first pregnancy;
	proc sort data=pregnancies_excl;
		by enrolid dt_lmp;
	run;
	%*output the first pregnancy - this is the primary analysis cohort;
	data &sens._cohort;
	format year_index $4.;
	set pregnancies_excl;
		by enrolid dt_lmp;

		%*Implement the restriction criteria;
		%if &where ne NA %then %do;
			where &where;
		%end;
		
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
	data &sens._cohort;
	set &sens._cohort;
		format preg_outcome_clean $pregoutc.
				thyroid_disorder_trt_post treated.
				bipolar_trt_post treated.
				schizo_trt_post treated.
				rural unknown.
				diabetes_type_post $diabetes.
				diabetes_simp $simpdiab.;
				
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



	%*************INVESTIGATE COVARIATE DISTRIBUTIONS;

	%*Now do some table 1 descriptives of the cohortt;
	%table1(inds = &sens._cohort, colVar = exposure,
		rowVars = ga_index_days ga_index_lt14 age_at_index year_index year_le2019 preg_outcome_clean chronichypertension_pre preeclampsia_pre
			substance_use_pre nausea_pre any_diabetes_post pregestation_diab_post diabetes_type_post 
			t2dmrx_post t1t2dmrx_post metforrx_post obesity_post glp1wgtrx_post otherwgtrx_post bariatric_post
			migraine_pre recurlos_pre ckd_post thyroid_disorder_post thyroidrx_post thyroid_disorder_trt_post
			depressi_post anxiety_post antideprx_post adhd_post adhdrx_post
			bipolar_post moodstabrx_post bipolar_trt_post ptsd_post schizo_post antipsyrx_post schizo_trt_post  
			rural hyperlip_post teratrx_pre benzorx_post anticonvulrx_post num_OutptPNC numOutptPNC_lt3 num_InptAdm  numInptADM_gt1,
		outfile = Table 1: Unweighted &sens ,
		title = Table 1: Sensitivity analysis cohort);
		
	%*Stratify by gestational age at index date;
	data primary_lt14 primary_ge14;
	set &sens._cohort;
		if ga_index_lt14 = 1 then output primary_lt14;
			else output primary_ge14;
	run;

	*Now output a table for each strata;
	%table1(inds = primary_lt14, colVar = exposure,
		rowVars = ga_index_days ga_index_lt14 age_at_index year_index year_le2019 preg_outcome_clean
			chronichypertension_pre preeclampsia_pre substance_use_pre nausea_pre diabetes_simp t2dmrx_post t1t2dmrx_post metforrx_post
			obesity_post migraine_pre recurlos_pre ckd_post thyroid_disorder_post thyroidrx_post thyroid_disorder_trt_post
			depressi_post anxiety_post antideprx_post adhd_post adhdrx_post rural hyperlip_post teratrx_pre benzorx_post anticonvulrx_post
			num_OutptPNC numOutptPNC_lt3,
		outfile = Table 1: Unweighted &sens GA lt 14wk,
		title = Table 1: Sensitivity analysis cohort where GA less than 14w at index);
		
	%*Now do some table 1 descriptives of the cohort - ana.primary_cohort;
	%table1(inds = primary_ge14, colVar = exposure,
		rowVars = ga_index_days ga_index_lt14 age_at_index year_index year_le2019 preg_outcome_clean
			chronichypertension_pre preeclampsia_pre substance_use_pre nausea_pre diabetes_simp t2dmrx_post t1t2dmrx_post metforrx_post
			obesity_post migraine_pre recurlos_pre ckd_post  thyroid_disorder_post thyroidrx_post thyroid_disorder_trt_post
			depressi_post anxiety_post antideprx_post adhd_post adhdrx_post rural hyperlip_post teratrx_pre benzorx_post anticonvulrx_post
			num_OutptPNC numOutptPNC_lt3,
		outfile = Table 1: Unweighted &sens GA ge 14wk,
		title = Table 1: Sensitivity analysis cohort where GA less than 14w at index);
		
	*Fit PS model within each strata to derive SMR weights;
	%competing2risk_weights(
		boot=0, 
		inds=&sens._cohort, 
		gacatvar = ga_index_cat,
		startDT=dt_index, 
		outcomevar = preg_outcome_clean,
		eventDT=dt_GApreg, 
		event = 'SAB' 'UAB' 'SB' 'MLS',
		cr1='IAB', 
		cr2='LBM' 'LBS' 'UDL',
		psvars=&psvars,
		dovars=&dovars,
		dovarsmodel = &dovarsmodel,
		psclassvars=&psclassvars, 
		doclassvars=&doclassvars,
		trtvar=&trtvar, 
		numiterations=1,
		initialseed=23244, 
		outds=&sens._point,
		outds_dist=ga_dist_&sens._point,
		outds_ps = ps_&sens._point,
		outds_surv = NA
		);

	*Attach these weights to peoples observations for a weighted Table 1 by strata;
	proc sql;
		create table weighted as 
		select a.*, case when missing(c.expwgt) then 1 else 2 end as ga_strat,
				coalesce(b.expwgt, c.expwgt) as smrw
		from &sens._cohort as a
		left join ana.ps_&sens._point_num_1 as b
		on a.idxpren=b.idxpren
		left join ana.ps_&sens._point_num_2 as c
		on a.idxpren=c.idxpren
		/*Subset to those that remained in the trimmed sample*/
		having idxpren in (select distinct idxpren 
							from ana.ps_&sens._point_num_1 
							union
							select distinct idxpren
							from ana.ps_&sens._point_num_2)
		;
		quit;

	%*Stratify by gestational age at index date;
	data primary_lt14 primary_ge14;
	set weighted;
		if ga_strat = 1 then output primary_lt14;
			else output primary_ge14;
	run;

	*Now output a table for each strata;
	%table1(inds = primary_lt14, colVar = exposure, wgtvar=smrw,
		rowVars = ga_index_days age_at_index year_index chronichypertension_pre preeclampsia_pre
			substance_use_pre nausea_pre diabetes_simp t2dmrx_post t1t2dmrx_post metforrx_post
			obesity_post migraine_pre recurlos_pre ckd_post thyroid_disorder_post thyroidrx_post thyroid_disorder_trt_post
			depressi_post anxiety_post antideprx_post adhd_post adhdrx_post
			hyperlip_post teratrx_pre benzorx_post anticonvulrx_post num_OutptPNC numOutptPNC_lt3,
		outfile = Table 1: Weighted &sens ga lt14,
		title = Table 1: Sensitivity analysis cohort where GA less than 14w at index);
		
	%*Now do some table 1 descriptives of the cohort - ana.primary_cohort;
	%table1(inds = primary_ge14, colVar = exposure, wgtvar=smrw,
		rowVars = ga_index_days age_at_index year_index chronichypertension_pre preeclampsia_pre
			substance_use_pre nausea_pre diabetes_simp t2dmrx_post t1t2dmrx_post metforrx_post
			obesity_post migraine_pre recurlos_pre ckd_post thyroid_disorder_post thyroidrx_post thyroid_disorder_trt_post
			depressi_post anxiety_post antideprx_post adhd_post adhdrx_post
			hyperlip_post teratrx_pre benzorx_post anticonvulrx_post num_OutptPNC numOutptPNC_lt3,
		outfile = Table 1: Weight &sens ga ge14,
		title = Table 1: Sensitivity analysis cohort where GA at least 14w at index);


	%*******03 - CONDUCT PREGNANCY LOSS ANALYSIS;
	
	%*Get the point estimate;
	%competing2risk_weights(
		boot=0, 
		inds=&sens._cohort, 
		gacatvar = ga_index_cat,
		startDT=dt_index, 
		outcomevar = preg_outcome_clean,
		eventDT=dt_GApreg, 
		event = 'SAB' 'UAB' 'SB' 'MLS',
		cr1='IAB', 
		cr2='LBM' 'LBS' 'UDL',
		psvars=&psvars,
		dovars=&dovars,
		dovarsmodel = &dovarsmodel,
		psclassvars=&psclassvars, 
		doclassvars=&doclassvars,
		trtvar=&trtvar,
		numiterations=1,
		initialseed=23244, 
		outds=ana.&sens._point,
		outds_dist=ga_dist_&sens._point,
		outds_ps = ps_&sens._point,
		outds_surv = NA
		);

	*Combine the stratified estimates according to the distribution of GA at index among those in the treated group;
	%combine_point_estimates(
		gadist = ana.ga_dist_&sens._point,
		numgastrat = 2,
		est = ana.&sens._point_,
		outds = ana.&sens._point_overall
		);

	*********Finally, conduct the bootstrap for standard errors;
	%competing2risk_weights(
		boot=1, 
		inds=&sens._cohort, 
		gacatvar = ga_index_cat,
		startDT=dt_index, 
		outcomevar = preg_outcome_clean,
		eventDT=dt_GApreg, 
		event = 'SAB' 'UAB' 'SB' 'MLS',
		cr1='IAB', 
		cr2='LBM' 'LBS' 'UDL',
		psvars=&psvars,
		dovars=&dovars,
		dovarsmodel = &dovarsmodel,
		psclassvars = &psclassvars, 
		doclassvars=&doclassvars,
		trtvar=&trtvar, 
		numiterations=&numboot,
		initialseed=23244, 
		outds=ana.&sens._boot,
		outds_dist=ga_dist_&sens._boot,
		outds_ps = ps_&sens._boot,
		outds_surv = NA
		);
		
	options mlogic mprint symbolgen notes;
	*Combine all the bootstrapped estimates;
	%combine_boot_estimates(
				inputEst= ana.&sens._boot, 
				inputDist= ana.ga_dist_&sens._boot,
				numStrata= 2, 
				output_stratified= ana.&sens._boot_strat,
				output_overall= ana.&sens._boot_overall);

	*Stratified estimates with confidence intervals.;
	%strat_estimates_w_CI(bootdsn=ana.&sens._boot_1, pointdsn=ana.&sens._point_1, output=ana.&sens._boot_1_ci);
	%strat_estimates_w_CI(bootdsn=ana.&sens._boot_2, pointdsn=ana.&sens._point_2, output=ana.&sens._boot_2_ci);

	*OVerall estimates with confidence interval;
	%overall_estimates_w_CI(stderrdsn=ana.&sens._boot_overall, pointdsn=ana.&sens._point_overall, output=ana.&sens._overall_boot_ci);


	%**********04 - CONDUCT PRETERM BIRTH ANALYSIS;

	%*get the point estimate;
	%competing3risk_weights(
		boot=0, 
		inds=&sens._cohort, 
		gacatvar = ga_index_cat,
		startDT=dt_index, 
		outcomevar = preg_outcome_ptb,
		eventDT=dt_GApreg_ptb, 
		event = 'PTB',
		cr1='IAB', 
		cr2='SAB' 'UAB' 'SB' 'MLS',
		cr3='TB',
		psvars=&psvars,
		dovars=&dovars,
		dovarsmodel = &dovarsmodel,
		psclassvars=&psclassvars, 
		doclassvars=&doclassvars,
		trtvar=&trtvar,
		numiterations=1,
		initialseed=23244, 
		outds=ana.&sens._point_ptb,
		outds_dist=ga_dist_&sens._point_ptb,
		outds_ps = ps_&sens._point_ptb,
		outds_surv = ana.surv_&sens._point_ptb
		);

	*Combine the stratified estimates according to the distribution of GA at index among those in the treated group;
	%combine_point_estimates(
		gadist = ana.ga_dist_&sens._point_ptb,
		numgastrat = 2,
		est = ana.&sens._point_ptb_,
		outds = ana.&sens._point_ptb_overall
		);

	%*Finally, conduct the bootstrap for standard errors;
	%competing3risk_weights(
		boot=1, 
		inds=&sens._cohort, 
		gacatvar = ga_index_cat,
		startDT=dt_index, 
		outcomevar = preg_outcome_ptb,
		eventDT=dt_GApreg_ptb, 
		event = 'PTB',
		cr1='IAB', 
		cr2='SAB' 'UAB' 'SB' 'MLS',
		cr3='TB',
		psvars=&psvars,
		dovars=&dovars,
		dovarsmodel = &dovarsmodel,
		psclassvars=&psclassvars, 
		doclassvars=&doclassvars,
		trtvar=&trtvar,
		numiterations=&numboot,
		initialseed=23244, 
		outds=ana.&sens._boot_ptb,
		outds_dist=ga_dist_&sens._boot_ptb,
		outds_ps = ps_&sens._boot_ptb,
		outds_surv = NA
		);	
		
	options mlogic mprint symbolgen notes;
	*Combine all the bootstrapped estimates;
	%combine_boot_estimates(
				inputEst= ana.&sens._boot_ptb, 
				inputDist= ana.ga_dist_&sens._boot_ptb,
				numStrata= 2, 
				output_stratified= ana.&sens._ptb_boot_strat,
				output_overall= ana.&sens._ptb_boot_overall);

	*Stratified estimates with confidence intervals.;
	%strat_estimates_w_CI(bootdsn=ana.&sens._boot_ptb_1, pointdsn=ana.&sens._point_ptb_1, output=ana.&sens._boot_ptb_1_ci);
	%strat_estimates_w_CI(bootdsn=ana.&sens._boot_ptb_2, pointdsn=ana.&sens._point_ptb_2, output=ana.&sens._boot_ptb_2_ci);

	*OVerall estimates with confidence interval;
	%overall_estimates_w_CI(stderrdsn=ana.&sens._ptb_boot_overall, pointdsn=ana.&sens._point_ptb_overall, output=ana.&sens._overall_ptb_boot_ci);

%mend;




	
