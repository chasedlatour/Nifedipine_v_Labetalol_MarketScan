/********************************************************************************************************************************************
PROGRAM: 19_descr_antihypertensives.sas
PROGRAMMER: Chase Latour
PURPOSE: To describe the use of antihypertensives after the index date.
	
Goal: 
Output data: 

Date: 12.19.2024
********************************************************************************************************************************************/






/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - DESCRIBE ANTIHYPERTENSIVES

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
%setup(sample=full, programname=19_descr_antihypertensives, savelog=N);

******************************************************************************************************************************************;
/*Create local mirros of the libraries from the set up macro - Run locally*/
/*libname lout slibref=out server=server; */
/*libname lana slibref=ana server=server;*/
/*libname lwork slibref=work server=server;*/
/*libname ltemp slibref=temp server=server;*/
/*libname lrxcov slibref=rxcov server=server;*/
/*libname lraw slibref=raw server=server;*/
/*libname lred slibref=red server=server;*/
/*libname lder slibref=der server=server;*/
/*libname lexpref slibref=expref server=server;*/

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

													01 - DESCRIBE ANTIHYPERTENSIVES

********************************************************************************************************************************************/

*The file temp.postidx_antihypertensives_63 contains information for all the antihypertensive fills for a pregnancy
between their index date and the outcome date. The 63 (9 weeks) represents the gestational age that we assumed at index for pregnancies
with UNK outcomes and no gestational age information. This should index lookback based on the first fill date (i.e., initiation), not
the LMP.;


*First, merge the treatment variable onto the dataset;
proc sql;
	create table pregnancies as
	select a.trt, a.dt_gapreg, b.*
	from ana.primary_cohort as a
	left join temp.postidx_antihypertensives_63 as b
	on  a.idxpren=b.idxpren
	;
	quit;


/**Look at the distribution of ATC labels within this dataset of antihypertensives between*/
/*the index date and pregnancy outcome or LTFU date;*/
/*proc freq data=pregnancies;*/
/*	table atc_label / missing;*/
/*run;*/



*Create counts of the number of people, by treatment, that filled a medication other than 
their initiation medication after their index date and prior to their pregnancy outcome or LTFU
date;
proc sql;
	create table non_initiation_counts as
	select distinct enrolid, idxpren, trt, sum(case when atc_label ne exposure then 1 else 0 end) > 0 as non_exposure_fill
	from pregnancies
	group by enrolid, idxpren, trt
	;
	quit;

proc freq data=non_initiation_counts;
	table trt*non_exposure_fill / missing;
run;
	



*Now, we want to know information about when the person discontinued their initial antihypertensive, if they did.

For this analysis, we define discontinuation as experiencing a gap in medications at home, based upon days supply
for the medication, allowing for a 30-day gap in fills.;

*Sort by svcdate for the antihypertensive fill. Subset to those rows where the antihypertensive filled
matches the initial antihypertensive;
proc sort data=pregnancies out=cohort ( where = (exposure = atc_label));
	by enrolid idxpren svcdate;
run;

*Confirm that everyone is still represented - Everyone should at least have thier first indexing fill;
proc sql;
	select count(distinct enrolid) as num_people, count(distinct idxpren) as num_preg
	from cohort;
	quit;

*Now use the days supply to evaluate continual usage. We want to be sure that if a fill occurs prior to the
	last fills day supply end, that we incorporate that overlapping time;
data cohort2;
set cohort;
	format next_fill MMDDYY10. last_date MMDDYY10. last_dt MMDDYY10.;
	by enrolid idxpren svcdate;
	retain last_date count30;

	last_id = lag1(idxpren);

	if first.idxpren then do;
		days_between = 0;
		next_fill = svcdate + daysupp;
		last_date = next_fill;
		end;
		else do;
			days_between = svcdate - last_date;

			if svcdate < last_date then do;
				next_fill = svcdate + daysupp + abs(days_between); *abs(days_between) + daysupp;
				last_date = next_fill;
			end;
				else do;
					next_fill = svcdate + daysupp;	
					last_date = next_fill;
				end; 
		end;

/*	last_dt = lag1(next_fill);*/

	drop last_date;

	if first.idxpren then count30 = 1;
		else if days_between > 30 then count30 = count30+1;
		else count30 = count30;

run;


%*Now the date of the first gap will be the next_fill date for the last row where count = 1;
data exposure;
set cohort2 (where = (count30 = 1));
	by enrolid idxpren svcdate;
	if last.idxpren then output;
run;

*Determine if the end of that gap occurred prior to the end of the pregnancy;
data exposure2;
set exposure;
	if next_fill < dt_gapreg then gap = 1;
		else gap = 0;
run;

proc freq data=exposure2;
	table trt*gap / missing;
run;




*How many of these discontinuers filled another antihypertensive prior to that date?;

*First, subset to the discontinuers -- those w gap in their initial exposure;
data discontinuers;
set exposure2;
	where gap = 1;

	drop last_date last_dt count30 last_id days_between;
run;

*Identify the first fill for an antihypertensive other than the exposure;
proc sql;
	create table non_exposure as
	select distinct enrolid, idxpren, min(svcdate) as first_non_exposure format=MMDDYY10.
	from pregnancies (where = (exposure ne atc_label))
	group by enrolid, idxpren
	;
	quit;

*Left join that information onto the discontinuers;
proc sql;
	create table discontinuers2 as
	select a.*, b.first_non_exposure format=MMDDYY10.
	from discontinuers as a
	left join non_exposure as b
	on a.enrolid=b.enrolid and a.idxpren=b.idxpren;
	quit;

*Determine the following among discontinuers
	- How many initiated a different antihyperensive at any point post-index
	- How many initiated a different antihypertensive prior to discontinuation of the first med;

data discontinuers3;
set discontinuers2;
	if first_non_exposure ne . then new_med = 1;
		else new_med = 0;

	if first_non_exposure < next_fill then new_med_before_gap  = 1;
		else new_med_before_gap = 0;

run;


proc freq data=discontinuers3;
	table trt*new_med trt*new_med_before_gap;
run;







/**I think that the original descriptive variables pulled for the medication usage (discontinuation dates) are not */
/*quite right.;*/
/**/
/**Count the number of rows in the antihypertensives rx dataset;*/
/*proc sql;*/
/*	select count(*) as n_row*/
/*	from temp.postidx_antihypertensives_63;*/
/*	quit;*/
/**/
/**/
/**/
/*proc sql;*/
/*	select count(distinct enrolid) as num_people, count(distinct idxpren) as num_preg*/
/*	from cohort;*/
/*	quit;*/
/**/
/**/
/*proc means data=ana.primary_cohort nmiss; */
/*	var dt_lastfill_gap0;*/
/*run;*/
/**/
/*	*/
/**/
/**/
/**/
/**/
/*proc freq data=rxcov.antihypertensives_rx;*/
/*	table atc_label;*/
/*run;*/
/**/
/**/
/**See if any of these IDs link in the one percent sample.;*/
/**/
/*data noantihyperensives;*/
/*set cohort;*/
/*	where idxpren = .;*/
/*run;*/
/**/
/*proc sql;*/
/*	create table noantihypertensives_1pct as*/
/*	select a.**/
/*	from der.enrlper as a*/
/*	inner join noantihyperensives as b*/
/*	on a.enrolid=b.enrolid;*/
/*	quit;*/
/**/
/**Now get that pregnancy row for the person;*/
/*proc sql;*/
/*	create table preg_in_1pct as*/
/*	select a.**/
/*	from ana.primary_cohort as a*/
/*	inner join noantihypertensives_1pct as b*/
/*	on a.enrolid=b.enrolid*/
/*	;*/
/*	quit;*/
/**/
/**/
/**/
/**Look at their fills in 2019 and make sure that they actually have an indexing fill.;*/
/*proc sql;*/
/*	create table fills as*/
/*	select a.idxpren, a.dt_indexprenatal, a.dt_index, a.trt, b.**/
/*	from preg_in_1pct as a*/
/*	left join raw.outptdrug2019 as b*/
/*	on a.enrolid=b.enrolid*/
/*	;*/
/*	quit;*/
/**/
/**/
/**/
/**/
/**/
/*proc sql;*/
/*	create table test as*/
/*	select a.**/
/*	from ana.primary_cohort as a*/
/*	inner join noantihyperensives as b*/
/*	on a.enrolid=b.enrolid;*/
/*	quit;*/
/**/
/**/
/*proc sql;*/
/*	create table any_fills_noantihypertensives as*/
/*	select b.**/
/*	from noantihypertensives_1pct as a*/
/*	inner join _postidx_antihypertensives as b*/
/*	on a.enrolid=b.enrolid*/
/*	;*/
/*	quit;*/
/**/
/**/
/**/
/*data postindex;*/
/*set temp.postidx_antihypertensives_63;*/
/*	where enrolid = 3438211802 or idxpren = 1348894;*/
/*run;*/
/**/
/**/
/*proc sql;*/
/*	select count(*) as n_row*/
/*	from _postidx_antihypertensives*/
/*	;*/
/*	quit;*/
/**/
/*data postindex;*/
/*set _postidx_antihypertensives;*/
/*	where enrolid = 3438211802 or idxpren = 1348894;*/
/*run;*/
/**/
/**/
/*proc sort data=_postidx_antihypertensives nodup; by enrolid idxpren; run;*/
/**/
/*data postindex2;*/
/*set _postidx_antihypertensives;*/
/*	where enrolid = 3438211802 or idxpren = 1348894;*/
/*run;*/
/**/
/**/
/*proc compare base=temp.preg_newusers_63_dt_index_270 compare=temp.preg_newusers_63_dt_lmp_180;*/
/*run;*/



















/*01 - DESCRIBE ANTIDIABETICS*/

/**First, look at the code list for metformin;*/
/*proc freq data=rxcov.metformin_rx;*/
/*	table atc_label / missing;*/
/*run;*/
/**/
/**Look at code list for T2DM meds;*/
/*proc freq data=rxcov.t2dm_antidiabetics_rx;*/
/*	table atc_label / missing;*/
/*run;*/
/**/
/**/
/*proc freq data=temp.cov_meds_63_dt_index_270;*/
/*	table medication;*/
/*run;*/
/**/
/*data antidiabetics;*/
/*set temp.cov_meds_63_dt_index_270;*/
/*	where medication in ('t2dm' 'Metfor' 't1t2dm');*/
/**/
/*	dt_index_year = year(dt_index);*/
/**/
/*	keep enrolid idxpren dt_index dt_index_year atc_label medication;*/
/*run;*/
/**/
/**Assign a class to the T2DM medications;*/
/**/
/*data antidiabetics2;*/
/*length class $45;*/
/*set antidiabetics;*/
/*	*/
/*	if medication = 't1t2dm' then class  = "insulin";*/
/*		else if medication = 'Metfor' then class = "metformin";*/
/*		else do;*/
/**/
/*			if atc_label = "ACARBOSE" then class = "Alpha glucosidase inhibitors";*/
/*				else if atc_label = "ACETOHEXAMIDE" then class = "Sulfonylureas";*/
/*				else if atc_label = "ALBIGLUTIDE" then class = "GLP-1";*/
/*				else if atc_label in ("ALOGLIPTIN BENZOATE", "ALOGLIPTIN BENZOATE/METFORMIN HYDROCHLORIDE") then class = "DPP4i";*/
/*				else if atc_label = "ALOGLIPTIN BENZOATE/PIOGLITAZONE HYDROCHLORIDE" then class = "DPP4i/Thiazolidinediones";*/
/*				else if atc_label in ("CANAGLIFLOZIN","CANAGLIFLOZIN/METFORMIN HYDROCHLORIDE") then class = "SGLT2";*/
/*				else if atc_label = "CHLORPROPAMIDE" then class = "Sulfonylureas";*/
/*				else if atc_label in ("DAPAGLIFLOZIN PROPANEDIOL", "DAPAGLIFLOZIN PROPANEDIOL/METFORMIN HYDROCHLORIDE") then class = "SGLT2";*/
/*				else if atc_label = "DAPAGLIFLOZIN/SAXAGLIPTIN" then class = "SGLT2/DPP4i";*/
/*				else if atc_label = "DULAGLUTIDE" then class = "GLP-1";*/
/*				else if atc_label in ("EMPAGLIFLOZIN", "EMPAGLIFLOZIN/METFORMIN HYDROCHLORIDE") then class = "SGLT2";*/
/*				else if atc_label in ("EMPAGLIFLOZIN/LINAGLIPTIN", "EMPAGLIFLOZIN/LINAGLIPTIN/METFORMIN HYDROCHLORIDE") then class= "SGLT2/DPP4i";*/
/*				else if atc_label in ("ERTUGLIFLOZIN", "ERTUGLIFLOZIN/METFORMIN HYDROCHLORIDE") then class = "SGLT2";*/
/*				else if atc_label = "ERTUGLIFLOZIN/SITAGLIPTIN" then class = "SGLT2/DPP4i";*/
/*				else if atc_label = "EXENATIDE" then class = "GLP-1";*/
/*				else if atc_label = "GLIMEPIRIDE" then class = "Sulfonylureas";*/
/*				else if atc_label in ("GLIMEPIRIDE/PIOGLITAZONE HYDROCHLORIDE", "GLIMEPIRIDE/ROSIGLITAZONE MALEATE") then class = "Sulfonylureas/Thiazolidinediones";*/
/*				else if atc_label in ("GLIPIZIDE", "GLIPIZIDE/METFORMIN HYDROCHLORIDE") then class = "Sulfonylureas";*/
/*				else if atc_label in ("GLYBURIDE", "GLYBURIDE, MICRONIZED", "GLYBURIDE/METFORMIN HYDROCHLORIDE", "TOLAZAMIDE") then class = "Sulfonylureas";*/
/*				else if atc_label in ("INSULIN DEGLUDEC/LIRAGLUTIDE","LIRAGLUTIDE", "SEMAGLUTIDE") then class = "GLP-1";*/
/*				else if atc_label in ("LINAGLIPTIN","LINAGLIPTIN/METFORMIN HYDROCHLORIDE") then class = "DPP4i";*/
/*				else if atc_label in ("PIOGLITAZONE HYDROCHLORIDE", "METFORMIN HYDROCHLORIDE/PIOGLITAZONE HYDROCHLORIDE") then class = "Thiazolidinediones";*/
/*				else if atc_label in ("METFORMIN HYDROCHLORIDE/REPAGLINIDE", "REPAGLINIDE", "NATEGLINIDE") then class = "Meglinitide";*/
/*				else if atc_label in ("SAXAGLIPTIN HYDROCHLORIDE", "METFORMIN HYDROCHLORIDE/SAXAGLIPTIN HYDROCHLORIDE",*/
/*										"SITAGLIPTIN PHOSPHATE" , "METFORMIN HYDROCHLORIDE/SITAGLIPTIN PHOSPHATE",*/
/*										"SIMVASTATIN/SITAGLIPTIN PHOSPHATE") then class = "DPP4i";*/
/*				else if atc_label in ("TOLBUTAMIDE", "TOLBUTAMIDE SODIUM") then class = "Sulfonylureas";*/
/*				else if atc_label in ("METFORMIN HYDROCHLORIDE/ROSIGLITAZONE MALEATE", "ROSIGLITAZONE MALEATE") then class = "Thiazolidinediones";*/
/*				else if atc_label = "MIGLITOL" then class = "Alpha glucosidase inhibitors";*/
/*		end;*/
/**/
/*run;*/
/**/
/*proc freq data=antidiabetics;*/
/*	table atc_label;*/
/*run;*/
/**/
/*proc sql;*/
/*	create table distinct_antidiabetics as*/
/*	select distinct enrolid, idxpren, atc_label*/
/*	from antidiabetics*/
/*	;*/
/*	quit;*/
/**/
/*proc freq data=distinct_antidiabetics;*/
/*	table atc_label;*/
/*run;*/
/**/
/**First just look at the distribution by treatment to see that we can get the same thing.;*/
/**/
/*proc sql;*/
/*	create table linked as*/
/*	select distinct a.trt, a.enrolid, a.idxpren, a.dt_index, year(a.dt_index) as dt_index_year, b.medication, b.class*/
/*	from ana.primary_cohort as a*/
/*	left join antidiabetics2 as b*/
/*	on a.enrolid = b.enrolid and a.dt_index = b.dt_index*/
/*	;*/
/*	quit;*/
/**/
/*proc freq data=linked;*/
/*	table class / missing;*/
/*run;*/
/**/
/**Create a dataset with the number of medications per person;*/
/*proc sql;*/
/*	create table med_person as*/
/*	select distinct enrolid, idxpren, trt, dt_index, dt_index_year, sum(medication ne "") as num_med,*/
/*		sum(class = "DPP4i") as num_dpp4, sum(class = "GLP-1") as num_glp1, sum(class = "SGLT2") as num_sglt2,*/
/*		sum(class = "Sulfonylureas") as num_sulfonylureas, sum(class = "Thiazolidinediones") as num_thiazolidinediones,*/
/*		sum(class = "insulin") as num_insulin, sum(class = "metformin") as num_metformin*/
/*	from linked */
/*	group by enrolid, idxpren, trt, dt_index, dt_index_year*/
/*	;*/
/*	quit;*/
/**/
/**/
/**Create a table that contains counts for each antidiabetic type;*/
/*proc sql;*/
/*	create table any_med as*/
/*	select trt, dt_index_year, count(*) as n_people, sum(num_med > 0) as any_med,*/
/*		sum(num_dpp4 > 0) as any_dpp4i, sum(num_glp1 > 0) as any_glp1,*/
/*		sum(num_sglt2 > 0) as any_sglt2, sum(num_sulfonylureas > 0) as any_sulfonylureas,*/
/*		sum(num_thiazolidinediones > 0) as any_thiazolidinediones,*/
/*		sum(num_insulin > 0) as any_insulin, sum(num_metformin) as any_metformin*/
/*	from med_person*/
/*	group by trt, dt_index_year;*/
/*	quit;*/
/**/
/*data any_med2;*/
/*set any_med;*/
/*	any_med_perc = cat(any_med, " (", round(100*any_med/n_people, 0.1), "%)");*/
/*	any_dpp4i_perc = cat(any_dpp4i, " (", round(100*any_dpp4i/n_people, 0.1), "%)");*/
/*	any_glp1_perc = cat(any_glp1, " (", round(100*any_glp1/n_people, 0.1), "%)");*/
/*	any_sglt2_perc = cat(any_sglt2, " (", round(100*any_sglt2/n_people, 0.1), "%)");*/
/*	any_sulfonylureas_perc = cat(any_sulfonylureas, " (", round(100*any_sulfonylureas/n_people, 0.1), "%)");*/
/*	any_thiazolidinediones_perc = cat(any_thiazolidinediones, " (", round(100*any_thiazolidinediones/n_people, 0.1), "%)");*/
/*	any_insulin_perc = cat(any_insulin, " (", round(100*any_insulin/n_people, 0.1), "%)");*/
/*	any_metformin_perc = cat(any_metformin, " (", round(100*any_metformin/n_people, 0.1), "%)");*/
/*	drop any_med any_dpp4i any_glp1 any_sglt2 any_sulfonylureas any_thiazolidinediones any_insulin any_metformin;*/
/*run;*/
/**/
/**/
/**/
/*data tolazamide;*/
/*set rxcov.t2dm_antidiabetics_rx;*/
/*	where atc_label = "TOLAZAMIDE";*/
/*run;*/
/**/
/*proc freq data=tolazamide;*/
/*	table ndcnum ndc9;*/
/*run;*/
/**/
/*proc sql;*/
/*	*/
/*	create table rx2020 as*/
/*	select a.*, b.**/
/*	from raw.outptdrug2020 as a*/
/*	inner join tolazamide as b*/
/*	on substr(a.ndcnum, 1, 9) = substr(b.ndcnum, 1, 9)*/
/*	;*/
/**/
/*	create table rx2021 as*/
/*	select a.*, b.**/
/*	from raw.outptdrug2021 as a*/
/*	inner join tolazamide as b*/
/*	on substr(a.ndcnum, 1, 9) = substr(b.ndcnum, 1, 9)*/
/*	;*/
/*	quit;*/
/**/
/*proc freq data=rx2020;*/
/*	table ndcnum / out=rx2020_table;*/
/*run;*/
/*proc freq data=rx2021;*/
/*	table ndcnum / out=rx2021_table;*/
/*run;*/
/**/
/*data red;*/
/*set red.redbook;*/
/*	where substr(ndcnum, 1, 9) = "498840122";*/
/*run;*/
/**/
/**At this point, discovered that a retired NDC code for TOLAZAMIDE was reassigned to labetalol.*/
/*This was causing huge differences in these variables between the two treatment groups.;*/
