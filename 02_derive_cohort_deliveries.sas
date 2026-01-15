/********************************************************************************************************************************************
PROGRAM: 02_derive_cohort_deliveries.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to extract covariate information for the pregnancy cohort we derived from marketScan claims data.
However, this analysis is restricted to deliveries only.
	
Goal: 
Output data: 

Date: 11.16.2025


********************************************************************************************************************************************/








/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - IDENTIFY PREGNANCIES WITH QUALIFYING FILL
	- 02 - IDENTIFY INCLUSION AND EXCLUSION CRITERIA
	- 03 - DEFINE COVARIATES

********************************************************************************************************************************************/










/********************************************************************************************************************************************

															00 - SET UP LIBRARIES

********************************************************************************************************************************************/

/*run local:*/
/*SIGNOFF;*/
/*%LET server=n2.schsr.unc.edu 1234; */
/*options comamid=tcp remote=server; */
/*signon username=_prompt_;*/

*No errors;;
*Run setup macro and define libnames;
options sasautos=(SASAUTOS "/local/projects/marketscan_preg/Latour_23_2322/programs/macros");
/*options mprint;*/

/*change "saveLog=" to "Y" when program is closer to complete*/
%setup(sample=full, programname=02_derive_cohort_deliveries, savelog=Y);



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

*Run this format file locally (Ctrl+A on the file with local run) if you want to view datasets;
/*%inc "/local/projects/marketscan_preg/Latour_23_2322/programs/FormatStatements_CDWH.sas";*/









*No errors;

/********************************************************************************************************************************************

											01 - IDENTIFY PREGNANCIES WITH QUALIFYING FILL

It is possible that we will want to consider different iterations of inclusion and exclusion criteria. However, in all cases, we want to
include pregnancies that have a qualifying nifedipine or labetalol fill. We first create that dataset and then we will add columns
that provide information on inclusion and exclusion criteria.

********************************************************************************************************************************************/


/*
MACRO: get_exposure
PURPOSE: Grab all the outpatient medication fills to identify those for nifedipine and labetalol.

This is written so that it will create the dataset _exp_subset, which will be appended (i.e., proc append) onto the initially 
empty dataset pregnancies_exp. Pregnancies_exp will be defined prior to applying the get_exposure macro.
*/

%macro get_exposure;

	%let numYr = %sysfunc(countw(&years));

	%do d=1 %to &numYr;

		%let loop&d = %scan(&years, &d);

		/*Create the dataset _exp_subset which contains all outpatient presription fills for nifedipine
		or labetalol which occur between 4 and 23 weeks of gestation*/
		proc sql;
			create table _exp_subset as
			select distinct a.enrolid as enrolid format=best12. length=8, a.ndcnum as ndcnum format=$12. length=12, 
				a.svcdate as svcdate format=MMDDYY10. length=8, a.age as age format=BEST12. length=8, 
				a.msa as msa format=Z5. length=4, a.region as region format=$3. length=3,
				a.egeoloc as egeoloc format=$3. length=3, a.eestatu as eestatu format=$3. length=3,
				a.metqty as metqty format=best12. length=4, a.daysupp format=best12. length=3,
				b.idxpren as idxpren format=BEST12. length=8, c.atc_label as atc_label format=$56. length=56,
				c.route as route format=$40. length=40, c.strength as strength format=$40. length=40,
				c.strength_uom format=$40. length=40, c.form as form format=$40. length=40
			from (select * from raw.outptdrug&&loop&d where metqty > 0) as a
			inner join _pregnancies as b
			on a.enrolid = b.enrolid and b.dt_lmp + 28 <= a.svcdate and a.svcdate <= b.max_fill
			inner join (select distinct atc_label, ndc9, route, strength, strength_uom, form from exposure) as c
			on substr(a.ndcnum, 1, 9) = c.ndc9
			;
			quit;

		%*Append the subset to the pregnancies dataset. This works even if 0 are identified.;
		proc append base=_pregnancies_exp data=_exp_subset; run;

	%end;

%mend;





/*
MACRO: derive_initial_cohort
PURPOSE: This macro identifies any pregnancies with at least 1 fill for nifedipine or labetalol between 4 and 28 weeks of gestation. 
This macro does not consider whether this represents initiation. That is an issue
accomplished in a subsequent macro.

LMPINDEX - The gestational age assigned to pregnancies with unknown outcomes without gestational age infomration
MAXFILL - The maximum gestational age at which an indexing fill can occur.
*/

%macro derive_initial_cohort(lmpindex);

	/*TESTING:
	%let lmpindex=63;
	%let maxfill=258;
	*/

	%*The clean_preg_cohort file limited to pregnancies in the study frame. Output the number of rows in that dataset
	EXCLUDE THOSE PREGNANCIES THAT DO NOT END IN DELIVERY;
	proc sql noprint;
	    select count(*) 
	    into :num_pregnancies_&lmpindex
	    from ana.preg_cohort_&lmpindex (where = (preg_outcome_clean in ('UDL' 'SB' 'LBS' 'LBM' 'MLS')));
		quit;

	%put Number of deliveries in the cohort before requiring medication fills before 22w0d days gestation: &&num_pregnancies_&lmpindex;



	/*******************
		Now, identify those individuals with a nifedipine or labetalol fill between their indexing prenatal encounter and 23_0 weeks gestation
	*******************/

	%*Conduct some initial data cleaning of the pregnancy dataset so that we can link outpatient drug claims to the relevant
	pregnancies. We want to subset, for now, to those pregnancies that have a defined index date that occurs prior to maxfill;
	data _pregnancies;
	set ana.preg_cohort_&lmpindex;
		%*Limit to those pregnancies that ended in delivery;
		where preg_outcome_clean in ('UDL' 'SB' 'LBS' 'LBM' 'MLS');

		%*Establish the date that the person was at the maximum gestational age for an index date;
		max_index = dt_lmp + 160; *22w 6d of gestation;


		%*Establish the latest date that someone can have a fill. This fill must occur prior to the delivery date;
		max_fill = min(dt_gapreg - 1, dt_lmp + 160);

	run;
	%*Make sure that there are no duplicates;
	proc sort data=_pregnancies nodup;
		by enrolid idxpren;
	run;


	%*Now, we want to get all fill records for the retained pregnancies that
		occurred between 4 and 22 weeks of gestation.
		For now, we focus only on those with fills for nifedipine or labetalol.;

	%*stack the labetalol and nifedipine datasets;
	data exposure;
	set expref.labetalol expref.nifedipine;
	run;

	%*Create an empty dataset that we are going to append the records too. We need to manually set the length
	and formats so that we avoid errors.; 
	data _pregnancies_exp;
	    length enrolid 8 ndcnum $12 svcdate 8 age 8 msa 4 region $3 egeoloc $3 eestatu $3 daysupp 3 metqty 4 idxpren 8 atc_label $56 route $40 strength $40 strength_uom $40 form $40;
	    format enrolid best12. ndcnum $12. svcdate MMDDYY10. age best12. msa Z5. region $3. egeoloc $3. eestatu $3. daysupp best12. metqty best12. idxpren best12. atc_label $56. route $40. strength $40. strength_uom $40. form $40.;
	    stop;
	run;

	%let years = 2016 2017 2018 2019 2020 2021 2022;
	%get_exposure;

	%*pregancies_exp now contains a dataset with all pregnancies that had a nifedipine or labetalol fill between
	their indexing prenatal encounter and 22w6d gestation.;

	%*Output the number of distinct pregnancies with indexing fill;
	proc sql noprint;
		select count(distinct idxpren) as num_pregnancies_nif_or_lab
		into :num_pregnancies_nif_or_lab
		from _pregnancies_exp;
		quit;

	%put Number of pregnancies with an indexing fill: &num_pregnancies_nif_or_lab;

	%*Now, subset to the first fill for each of the pregnancies and output that first fill, as this is now their index date;
	proc sort data=_pregnancies_exp;
		by idxpren svcdate;
	run;
	data _index;
	set _pregnancies_exp;
		by idxpren svcdate;
		if first.idxpren then output;
	run;

	%*Link this information back onto the pregnancy dataset;
	proc sql;
		create table temp.pregnancies_del_with_index_&lmpindex as
		select distinct a.*, b.svcdate as dt_index, b.ndcnum as ndc11, b.atc_label as exposure, b.age as age_at_index,
			b.msa as msa_at_index format=Z5. length=4, b.region as region_at_index format=$3. length=3,
			b.egeoloc as egeoloc_at_index format=$3. length=3, b.eestatu as eestatu_at_index format=$3. length=3,
			b.metqty as metqty_at_index format=best12. length=4, b.route as exposure_route format=$40. length=40,
			b.strength as exposure_strength format=$40. length=40, b.strength_uom as exposure_strength_uom format=$40. length=40,
			b.form as exposure_form format=$40. length=40, b.daysupp format=best12. length=3
		from _pregnancies as a
		inner join _index as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		having a.dt_lmp + 28 <= b.svcdate <= a.max_fill /*Only include those individuals whose fills occurred in gestational age range*/
		;
		quit;

%mend;










**************************************************************************************************
Now, subset those datasets to new users. This is the final dataset to which we will be adding
the information on exclusion and inclusion criteria.;


/*
MACRO: get_antihypertensives
PURPOSE: Go through all outpatient fills to identify those fills for antihypertensives
that are either nifedipine, labetalol, or neither. This loops through all the separate
year files from the raw data. In particular, we want to identify all those fills that occurred
between the lookback date and the index date.
*/
%macro get_antihypertensives;

	%let numYr = %sysfunc(countw(&years));

	%do d=1 %to &numYr;

		%let loop&d = %scan(&years, &d);

		%*This is a subset of all the antihypertensive fills where a dispensing quantity was greater than 0. This is
		created for each year and then appended (i.e., proc append) onto the _prior_antihypertensives;
		proc sql;
			create table _antihtn_subset as
			select a.enrolid as enrolid format=best12. length=8, a.svcdate as svcdate format=MMDDYY10. length=8,
				b.idxpren as idxpren format=BEST12. length=8, b.dt_index as dt_index format=MMDDYY10. length=8, 
				b.exposure as exposure format=$56. length=56, b.dt_lookback as dt_lookback format=MMDDYY10. length=8,
				c.atc_label as atc_label format=$56. length=56
			from (select distinct * from raw.outptdrug&&loop&d where metqty > 0) as a
			inner join _pregnancies_index as b /*Only grab those fills for pregnancies in the current dataset*/
			on a.enrolid = b.enrolid and b.dt_lookback <= a.svcdate <= b.dt_index /*Fills must occur between the lookback date and the index date*/
			inner join (select distinct atc_label, ndc9 from rxcov.antihypertensives_rx) as c /*Link based upon NDC-9 code, as there may be missed packagings with NDC-11*/
			on substr(a.ndcnum, 1, 9) = c.ndc9
			;
			quit;

		proc append base=_prior_antihypertensives data=_antihtn_subset; run;

		/*Output the number of pregnancies in that year*/
		proc sql noprint;
			select count(*) as num_pregnancies_loop_year
			into :num_pregnancies_loop_year
			from _antihtn_subset;
			quit;

		%put Number of pregnancies with qualifying antihypertensives in &&loop&d: &num_pregnancies_loop_year;

	%end;

%mend;


/*******************
MACRO: identify_new_users
PURPOSE: To identify those individuals that were new users of antihypertensives according to the presence
of any antihypertensive fill 270 days prior to the index date. Further, we want to remove those individuals
who have a fill for an antihypertensive other than their assigned exposure on the same date.

INPUTS:
LMPINDEX -- The gestatational age assumed at the index date for those pregnancies with UNK outcomes and no GA information
LOOKBACK_DT -- The date that we will use to anchor a lookback window. This is either dt_index (primary)
	or dt_lmp (sensitivity analysis)
LOOKBACK_DAYS -- The number of days to look back from the lookback_dt (270 for primary and 180 based on LMP)
******************/
%macro identify_new_users(lmpindex, lookbackdt, lookbackdays);

	/* TESTING:
	%let lmpindex=63;
	%let lookbackdt=dt_index;
	%let lookbackdays=270;
	*/

	%*Pull in the dataset that we are going to be using;
	data _pregnancies_index;
	set temp.pregnancies_with_index_&lmpindex;
		length dt_lookback 8;
		format dt_lookback MMDDYY10.;
		
		*Create the lookback date: This the maximum date from which we should look for covariate information;
		dt_lookback = &lookbackdt - &lookbackdays;
	run;

	%*Create an empty dataset that we are going to append the records too. We need to manually set the length
	and formats so that we avoid errors.;
	data _prior_antihypertensives;
	    length enrolid 8 svcdate 8 idxpren 8 dt_index 8 dt_lookback 8 exposure $56 atc_label $56;
	    format enrolid best12. svcdate MMDDYY10. idxpren best12. dt_index MMDDYY10. dt_lookback MMDDYY10. exposure $56. atc_label $56.;
	    stop;
	run;

	%let years = 2016 2017 2018 2019 2020 2021 2022;
	%get_antihypertensives;

	%*Now we want to identify the number that have a fill prior to their index date
	OR a fill for a medication OTHER THAN the exposure on their index date;
	proc sql;
		create table _prior_antihypertensives2 as
		select distinct enrolid, svcdate, idxpren, dt_index, dt_lookback, exposure, atc_label,
			sum(svcdate < dt_index) as num_pre_index,
			sum(svcdate = dt_index and atc_label ne exposure) as num_others_index_date
		from _prior_antihypertensives
		group by idxpren
		order by idxpren, svcdate
		;
		quit;

	%*****Provide counts in the log;
	proc sql noprint;
		select sum(num_pre_index > 0) as num_preg_pre_index_fill,
				sum(num_others_index_date > 0) as num_preg_other_fill_on_index,
				sum(num_pre_index > 0 and num_others_index_date > 0) as num_preg_both
		into :num_preg_pre_index_fill, :num_preg_other_fill_on_index, :num_preg_both
		from (
				select idxpren, min(num_pre_index) as num_pre_index, min(num_others_index_date) as num_others_index_date
				from _prior_antihypertensives2
				group by idxpren
			  )
		;
		quit;

	%put Number of pregnancies with an antihypertensive fill prior to the index date: &num_preg_pre_index_fill;
	%put Number of pregnancies with an antihypertensive fill other than the exposure on the index date: &num_preg_other_fill_on_index;
	%put Number of pregnancies with a pre-index antihypertensive fill and another antihypertensive fill on the index date: &num_preg_both;

	%*Subset to new-users;
	data _new_users;
	set  _prior_antihypertensives2;
		where num_pre_index = 0 and num_others_index_date = 0;
	run;

	%*Output a count of new-users;
	proc sql noprint;
		select count(distinct idxpren)
		into :num_new_users
		from _new_users;
		quit;
	%put Number of pregnancies identified as new users: &num_new_users;

	%*Now subset to those pregnancies that are identified as new-users using an inner join;
	proc sql;
		create table temp.preg_del_newusers_&lmpindex._&lookbackdt._&lookbackdays as
		select distinct a.*
		from _pregnancies_index as a
		inner join _new_users as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		;
		quit;

%mend;















/********************************************************************************************************************************************

											02 - IDENTIFY INCLUSION AND EXCLUSION CRITERIA

NOTE: Some of our covariates use information up to 30 days after the index date, particularly chronic conditions. However, we do not want
to capture information after the outcome or LTFU date nor the end of continuous enrollment. Thus, the maximum time at which covariate 
information can be captured post-index is the minimum of dt_index+30 and dt_gapreg, unless otherwise specified.

********************************************************************************************************************************************/


/*
MACRO: get_meds
PURPOSE: The purpose fo this macro is to grab the medication fills with NDC9s that match
those in a reference dataset.

INPUTS:
INPUT_PREGS -- Input dataset with all of the pregnancies that we are evaluating.
INPUT_MEDS -- Reference dataset with medication information, which contains NDC codes for the medications assocaited with a covariate. Derived
from the ATC-NDC file from the First Data Bank.
BASE -- Empty dataset to which the pulled dataset will be added (i.e., via proc append)
*/

%macro get_meds(input_pregs, input_meds, base);

	%let numYr = %sysfunc(countw(&years));

	%do d=1 %to &numYr;

		%let loop&d = %scan(&years, &d);

		/*Create the dataset _med_subset which contains all outpatient presription fills for the input meds
		between 270 days prior to the indexing prenatal encounter and 30 days after the indexing fill*/
		proc sql;
			create table _med_subset as
			select distinct a.enrolid as enrolid format=best12. length=8, a.svcdate as svcdate format=MMDDYY10. length=8, 
				b.idxpren as idxpren format=BEST12. length=8, b.dt_index as dt_index format=MMDDYY10. length=8,
				b.dt_lookback as dt_lookback format=MMDDYY10. length=8, 
				b.dt_lmp as dt_lmp format MMDDYY10. length=8,
				c.atc_label as atc_label format=$56. length=56, c.medication as medication format=$10. length=10,
				b.dt_gapreg as dt_gapreg format MMDDYY10. length=8
			from (select * from raw.outptdrug&&loop&d. where metqty > 0) as a
			inner join &input_pregs as b /*Only select those fills for those in the covariate dataset*/
			on a.enrolid = b.enrolid and b.dt_lookback <= a.svcdate <= min(b.dt_index+30, b.dt_gapreg)
			inner join (select distinct atc_label, ndc9, medication from &input_meds) as c
			on substr(a.ndcnum, 1, 9) = c.ndc9
			;
			quit;

		/*Append the pulled dataset to the base dataset*/
		proc append base=&base data=_med_subset; run;

		%*Output counts;
		proc sql noprint;
			select count(distinct idxpren) as num_preg_w_fill
			into :num_preg_w_fill
			from _med_subset;
			quit;
		%put Number of pregnancies with at least 1 fill from &&loop&d : &num_preg_w_fill;

	%end;

%mend;





/*
MACRO: get_dx_covariates
PURPOSE: The purpose of this macro is to get all the relevant diagnosis codes from the derived files. Relevant diagnosis codes correspond
to those on the reference file for at least one of the covariates. This pulls from the derived diagnosis code files created by Virginia
Pate. This is a file where each row is a diagnosis code from an inpatient service, outpatient service, or inpatient admission claim. The
type of claim is indicated in the dataset.

INPUTS:
INPUT_PREGS -- Input dataset with all of the pregnancies
INPUT_REF -- Input dataset with all of the reference codes for the defined diagnoses. Expects a column called diagnosis
BASE -- Empty dataset to which the pulled dataset will be added via proc append
*/

%macro get_dx_covariates(input_pregs, input_ref, base);

	%let numYr = %sysfunc(countw(&years));

	%do d=1 %to &numYr;

		%let loop&d = %scan(&years, &d);

		%*Create dataset with all relevant diagnosis codes for identified pregnancy for the year defined by &loop&d;
		proc sql;
			create table _diagnoses_subset as
			select distinct a.enrolid format=best12. length=8, a.svcdate format=MMDDYY10. length=8, 
				c.idxpren format=best12. length=8, c.dt_index format=MMDDYY10. length=8, 
				a.dxLoc format=$9. length=9, a.dxNum format=best12. length=8, b.code format=$7. length=7, 
				b.diagnosis format=$8. length=8
			from %if &&loop&d = 2016 %then %do; der.alldx102016 %end; %else %do; der.alldx&&loop&d. %end; as a
			inner join &input_pregs as c
			on a.enrolid=c.enrolid and c.dt_lookback <= a.svcdate <= min(c.dt_index + 30, c.dt_gapreg)
			inner join &input_ref as b
			on a.dx&&loop&d. = b.code
			;
			quit;

		%*Append the dataset onto the base dataset;
		proc append base=&base data=_diagnoses_subset; run;

		%*Output counts into the log;
		proc sql noprint;
			select count(distinct idxpren) as num_preg_w_dx
			into :num_preg_w_dx
			from _diagnoses_subset;
			quit;
		%put Number of pregnancies with at least 1 diagnosis code from &&loop&d : &num_preg_w_dx;

	%end;

%mend;




/*
MACRO: get_chd_procedures
PURPOSE: The purpose of this macro is to get all the procedure codes from the derived procedure code file (same as dx via Virginia Pate)
that correspond to fetal echocardiography. This will be used to derive whether a congenital heart diagnosis refers to 
the pregnant person or the fetus. Of note, we are only interested in fetal echocardiography procedures after the LMP
in case there are prior pregnancies in the lookback window.

INPUTS:
INPUT_PREGS -- Input dataset with all of the pregnancies
INPUT_PROC -- Input dataset with all of the reference codes for the defined procedures
BASE -- Empty dataset to which the pulled dataset will be added
*/

%macro get_chd_procedures(input_pregs, input_proc, base);

	%let numYr = %sysfunc(countw(&years));

	%do d=1 %to &numYr;

		%let loop&d = %scan(&years, &d);

		%*Create the dataset with fetal echocardiography procedure codes for the year defined by &&loop&d;
		proc sql;
			create table _proc_subset as
			select distinct a.enrolid as enrolid format=best12. length=8, a.svcdate as svcdate format=MMDDYY10. length=8,
					b.idxpren as idxpren format=best12. length=8, a.proc&&loop&d as code format=$7. length=7,
					b.dt_index as dt_index format=MMDDYY10. length=8
			from der.allproc&&loop&d as a
			inner join &input_pregs as b
			on a.enrolid=b.enrolid and b.dt_lmp <= a.svcdate <= min(b.dt_index + 30, b.dt_gapreg)
			inner join (select * from &input_proc where codetype="CPT") as c
			on a.proc&&loop&d = c.code
			;
			quit;

		proc append base=&base data=_proc_subset; run;

		%*Output counts;
		proc sql noprint;
			select count(distinct idxpren) as num_preg_fetal
			into :num_preg_fetal
			from _proc_subset;
			quit;
		%put Number of pregnancies with at least 1 fetal echocardiography for &&loop&d : &num_preg_fetal;

	%end;

%mend;




/*
MACRO: cardiovascular_diagnoses
PURPOSE: The purpose of this macro is to get all the inpatient and outpatient services records where the 
STDPROV variable is equal to 250, 440, 540, 580. These are cardiovascular specialties. This is particularly
important for defining maternal congenital heart defects. This is derived from the raw files provided by
MarketScan.

INPUTS:
INPUT_PREGS -- Input dataset with all of the pregnancies
BASE -- Empty dataset to which the pulled dataset will be added
*/

%macro cardiovascular_diagnoses(input_preg, base);

	%let numYr = %sysfunc(countw(&years));

	%do d=1 %to &numYr;

		%let loop&d = %scan(&years, &d);

		/*Grab all diagnosis codes for claims where the provider type was of those considered in our study*/
		proc sql;
			create table _cardiovascular_diagnoses_subset as
			/*First, grab the records from the outpatient services files*/
			select distinct *
			from (
				select a.enrolid as enrolid format=best12. length=8, a.svcdate as svcdate format=MMDDYY10. length=8,
					b.idxpren as idxpren format=best12. length=8, b.dt_index as dt_index format=MMDDYY10. length=8,
					/*WEnt through all files and confirmed only 4 diagnosis codes on each out patient service claim*/
					a.dx1 as dx1 format=$7. length=7, a.dx2 as dx2 format=$7. length=7,
					a.dx3 as dx3 format=$7. length=7, a.dx4 as dx4 format=$7. length=7
				from (select * from raw.outptserv&&loop&d (keep = enrolid svcdate stdprov dx:) where dxver = '0' and STDPROV in (250 440 540 580)) as a
					/*dxver = 0 is limiting to ICD-10 diagnosis codes*/
				inner join &input_preg as b
				on a.enrolid=b.enrolid and b.dt_lookback <= a.svcdate <= min(b.dt_index + 30, b.dt_gapreg)
			)
			union corr
			/*Union with those records from the inpatient services files*/
			select *
			from (
				select a.enrolid as enrolid format=best12. length=8, a.svcdate as svcdate format=MMDDYY10. length=8,
					b.idxpren as idxpren format=best12. length=8, b.dt_index as dt_index format=MMDDYY10. length=8,
					/*WEnt through all files and confirmed only 4 diagnosis codes on each out patient service claim*/
					a.dx1 as dx1 format=$7. length=7, a.dx2 as dx2 format=$7. length=7,
					a.dx3 as dx3 format=$7. length=7, a.dx4 as dx4 format=$7. length=7
				from (select * from raw.inptserv&&loop&d (keep = enrolid svcdate stdprov dx:) where dxver = '0' and STDPROV in (250 440 540 580)) as a
				inner join &input_preg as b
				on a.enrolid=b.enrolid and b.dt_lookback <= a.svcdate <= min(b.dt_index + 30, b.dt_gapreg)
			)
			;
			quit;

		proc append base=&base data=_cardiovascular_diagnoses_subset; run;

		%*Output counts;
		proc sql noprint;
			select count(distinct idxpren) as num_preg_cvdenc
			into :num_preg_cvdenc
			from _cardiovascular_diagnoses_subset;
			quit;
		%put Number of pregnancies with a cardiovascular health encounter in &&loop&d : &num_preg_cvdenc;

	%end;

%mend;






/*
MACRO: identify_dx_covariates
PURPOSE: TO identify the relevant covariates for defining inclusion and excusion criteria using only diagnosis codes.
Specifically, this macro creates variables for:
- The number of distinct svcdates for inpatient and outpatient diagnosis codes for each diagnosis, pre and post index
- The number of distinct svcdates where the diagnosis codes were in the primary or secondary position for an inpatient
	admission, pre and post index
- The number of distinct svcdates for inpatient and outpatient diagnosis codes that were at least 30 days apart, pre and
	post index. The 30-day requirement is implemented in the pre and post index period separately.

INPUT:
- REF: Reference code list with columns indicating the code and diagnosis associated with the code
- PREG_INPT: Input dataset with the base pregnancies that we want to append diagnosis code-based information to
- OUTPUT: Name for the desired name of the pregnancy-level dataset output by the macro
*/

%macro identify_dx_covariates(REF=, PREG_INPT=, OUTPUT=);

	%*Create an empty dataset that we are going to append the records too. We need to manually set the length
	and formats so that we avoid errors.;
	data _diagnoses_covariate;
		length enrolid 8 svcdate 8 idxpren 8 dt_index 8 dxLoc $9 dxNum 8 code $7 diagnosis $8;
		format enrolid best12. svcdate MMDDYY10. idxpren best12. dt_index MMDDYY10. dxLoc $9. dxnum best12. code $7. diagnosis $8.;
		stop;
	run;

	%get_dx_covariates(input_pregs=&PREG_INPT, input_ref = &REF, base = _diagnoses_covariate);

	%*Save this dataset for subsequent review - diabetes results are not as expected;
	data temp.covdiagnoses_&lmpindex._&lookbackdt._&lookbackdays; set _diagnoses_covariate; run;

	%*Do some data cleaning on the codes;
	data _diagnoses2;
	set _diagnoses_covariate;
		%*Determine if the diagnosis code was on an outpatient record (outpatient service claim) or inpatient record (inpatient service or
			inpatient admission claim);
		if dxLoc = "OutptServ" then location = "Outpt";
			else if dxLoc in ("InptServ" "InptAdm") then location = "Inpt";
			else location = "";
		%*Create a variable to indicate if the diagnosis code occurred prior to or after the index date;
		if svcdate > dt_index then timing = "post";
			else if svcdate <= dt_index then timing = "pre";
			else timing= "";

		if location = "" then delete; /*Should not be a problem*/
	run;



	%********************
		STEP 1: Identify counts for codes that occur on different service dates
	*********************;

	%*Count the number of diagnosis codes that occur on different dates for each diagnosis, stratified by code location (inpatient or outpatient) 
	and timing (pre or post index);
	proc sql;
		create table _diagnoses3 as
		select distinct enrolid, idxpren, diagnosis, location, timing, count(distinct svcdate) as numdx /*These codes need to be on different dates*/
		from _diagnoses2
		group by enrolid, idxpren, diagnosis, location, timing
		;
		quit;

	%*Transpose the dataset so that each diagnosis, location, and timing  (pre v post index) are their own column.
	These columns contain counts of the diagnoses that satisfy the requirement.
	NOTE: If no one has a code that meets a criteria, then that column will not be created.;
	proc transpose data=_diagnoses3 out=_diagnoses4 (drop = _NAME_) prefix=dx_;
		by enrolid idxpren;
		id diagnosis location timing;
		var numdx;
	run;

	%*Finally, join all of the dx variables onto the pregnancy dataset. This will produce a warning because enrolid and 
	idxpren are represented in both datasets. However, this is substantially simpler than referencing the dictionary
	to grab the names of all columns that start with dx_.
	IF A COLUMN IS NOT PRESENT, IT MEANS THAT NO RELEVANT DIAGNOSES WERE IDENTIFIED.;
	proc sql;
		create table _pregnancies3 as
		select distinct a.*, b.*
		from &PREG_INPT as a
		left join _diagnoses4 as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		;
		quit;

	%*All variables with missing values for the dx_ need to be replaced with 0.;
	data _pregnancies4_a;
	set _pregnancies3;
		array dx_vars{*} dx_:;
		do i=1 to dim(dx_vars);
			if missing(dx_vars[i]) then dx_vars[i]=0;
		end;
		drop i;
	run;




	%********************
		STEP 2: Identify codes that occur in the primary or secondary posiitons of inpatient admissions
	*********************;

	%*Count the number of diagnosis codes in the primary and secondary position of an inpt adm (different dates);
	proc sql;
		create table _diagnoses3_ps as
		select distinct enrolid, idxpren, diagnosis, location, timing, count(distinct svcdate) as numdx /*These codes need to be on different dates*/
		from (select * from _diagnoses2 where dxLoc = "InptAdm" and dxNum <= 2) /*This subsets to primary or secondary positions of codes*/
		group by enrolid, idxpren, diagnosis, location, timing
		;
		quit;

	%*Transpose the dataset so that each diagnosis, location, and timing  (pre v post index) are their own column.
	These columns contain counts of the diagnoses that satisfy the requirement.
	NOTE: If a code is not foundn that meets a criteria, that meets that NO ONE has information that supports it.;
	proc transpose data=_diagnoses3_ps out=_diagnoses4_ps (drop = _NAME_) prefix=dx_ps_;
		by enrolid idxpren;
		id diagnosis location timing;
		var numdx;
	run;

	%*Finally, join all of the dx variables onto the pregnancy dataset.
	This will produce a warning because enrolid and idxpren are represented
	in both datasets. However, this is substantially simpler than referencing the dictionary
	since we do not know what all columns there will be.
	IF A COLUMN IS NOT PRESENT, IT MEANS THAT NO RELEVANT DIAGNOSES WERE IDENTIFIED.;
	proc sql;
		create table _pregnancies3_ps as
		select distinct a.*, b.*
		from _pregnancies4_a as a
		left join _diagnoses4_ps as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		;
		quit;

	%*All variables with missing values for the dx_ need to be replaced with 0.;
	data _pregnancies4_b;
	set _pregnancies3_ps;
		array dx_vars{*} dx_:;
		do i=1 to dim(dx_vars);
			if missing(dx_vars[i]) then dx_vars[i]=0;
		end;
		drop i;
	run;







	%********************
		STEP 3: Indentify codes that occur at least 30 days apart
	*********************;

	%*Repeat the original diagnosis code process but identify those that are >= 30 days apart;
	
	%*Do some data cleaning to identify those codes that are at least 30d apart;
	proc sort data=_diagnoses2;
		by enrolid idxpren timing svcdate;
	run;
	data _diagnoses_30d;
	set _diagnoses2;
		by enrolid idxpren timing svcdate;
		retain count;

		last_id = lag1(idxpren);
		last_dt = lag1(svcdate);
		days_elapsed = svcdate - last_dt;

		if first.idxpren then count = 1;
			else if idxpren = last_id and days_elapsed > 30 then count = count + 1;

	run;

	%*Count the number of diagnosis that occur on different dates, stratified by inpatient and outpatient codes;
	proc sql;
		create table _diagnoses3 as
		/*CDL: 2.4.2025 -- Changed count(distinct cats(svcdate, count)) to count(distinct count) based on MPG note that gave incorrect results.*/
		select distinct enrolid, idxpren, diagnosis, location, timing, count(distinct count) as numdx /*These codes need to be on different dates at least 30 days apart*/
		from _diagnoses_30d
		group by enrolid, idxpren, diagnosis, location, timing
		;
		quit;

	%*Transpose the dataset so that each diagnosis, location, and timing  (pre v post index) are their own column.
	These columns contain counts of the diagnoses that satisfy the requirement.
	NOTE: If a column is not found, that meets that NO ONE has information that supports it.;
	proc transpose data=_diagnoses3 out=_diagnoses4 (drop = _NAME_) prefix=dx_30_;
		by enrolid idxpren;
		id diagnosis location timing;
		var numdx;
	run;

	%*Finally, join all of the dx variables onto the pregnancy dataset. This will produce a warning because enrolid and idxpren are represented
	in both datasets. However, this is substantially simpler than referencing the dictionary since we do not know what all columns there will be.
	IF A COLUMN IS NOT PRESENT, IT MEANS THAT NO RELEVANT DIAGNOSES WERE IDENTIFIED.;
	proc sql;
		create table _pregnancies30d as
		select distinct a.*, b.*
		from _pregnancies4_b as a
		left join _diagnoses4 as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		;
		quit;

	%*All variables with missing values for the dx_ need to be replaced with 0.;
	data &OUTPUT;
	set _pregnancies30d;
		array dx_vars{*} dx_:;
		do i=1 to dim(dx_vars);
			if missing(dx_vars[i]) then dx_vars[i]=0;
		end;
		drop i;
	run;

	%*Delete the unnecessary datasets;
	proc datasets gennum=all noprint;
		delete _diagnoses_covaraites _diagnoses3 _diagnoses4 _diagnoses_subset
			_pregnancies3 _pregnancies4_a 
			_diagnoses3_ps _diagnoses4_ps _pregnancies3_ps _pregnancies4_b
			_diagnoses_30d _diagnoses_ps _pregnancies_30d;
	run;


%mend;






/*
MACRO: define_incl_excl
PURPOSE: To define covariates required to apply the inclusion and exclusion criteria to the pregnancy cohort.

Inputs: 
lmpindex -- gestational age assumed at the index date for pregnancies with UNK outcomes and no GA information
gap -- Number of day gap period used to define a period of continuous enrollment prior to the index date
lookbackdt -- date being used to index the lookback period from: dt_index (primary) or dt_lmp (sensitivity analysis)
lookbackdays -- the number of days we are looking back from lookbackdt to assign covariates, etc.
*/


%macro define_incl_excl(lmpindex, gap, lookbackdt, lookbackdays);

	/*Testing
	%let lmpindex=63;
	%let gap = 31;
	%let lookbackdt = dt_index;
	%let lookbackdays = 270;
	*/

	%let years = 2016 2017 2018 2019 2020 2021 2022;

	/*******
	STEP 1: Call in the cohort that we are going to add the covariate values to
	********/

	data _pregnancies;
	set temp.preg_newusers_&lmpindex._&lookbackdt._&lookbackdays;
	run;


	/********
	STEP 2a: Add information  related to enrollment in a FFS plan with RX coverage.
	********/

	%*Subset to those continuous enrollment periods for those individuals in our pregnancy dataset via inner join.
	This is using the continuous enrollment derived files created by Virginia Pate. For the exclusoin criteria, 
	we are interested in enrollment in a Fee-for-service (FFS) plan with prescription (rx) coverage.;
	proc sql;
		create table _enrollment as
		select distinct b.patient_deid, a.enrolid, a.start, a.end, b.idxpren, b.dt_index, b.dt_lookback, b.dt_gapreg
		from der.rxenrlper_ffs as a
		inner join _pregnancies as b
		on a.enrolid = b.enrolid
		;
		quit;

	%*Output some relevant counts;
	proc sql noprint;
		/*Number of pregnancies prior to implementing any continous enrollment requirements.*/
		select count(distinct idxpren) as num_preg_pre_enrl
		into :num_preg_pre_enrl
		from _pregnancies;
		/*Number of pregnancies that had no enrollment information*/
		select count(distinct idxpren) as num_preg_any_enrl_rxffs
		into :num_preg_any_enrl_rxffs
		from _enrollment;
		quit;
	%put Number of pregnancies in our initial dataset: &num_preg_pre_enrl;
	%put Number of pregnancies without any FFS Rx enrollment information: &num_preg_any_enrl_rxffs;

	%*****Implement allowable gap;

	%*Calculate the time elapsed from the last continuous enrollment period;
	proc sort data=_enrollment;
		by patient_deid enrolid idxpren start;
	run;
	data _enrollment2;
	set _enrollment;
		by patient_deid enrolid idxpren start;
		retain enrl_count;
		
		last_id = lag1(idxpren);
		elapsed = start - lag1(end);

		%*Carry down the information through a count variable;
		if first.idxpren then enrl_count = 1;
			else if idxpren = last_id and elapsed > &gap then enrl_count = enrl_count + 1;
	run;

	%*Output the continuous enrollment rows, revised with the gap;
	proc sql;
		create table _enrollment3 as
		select distinct patient_deid, enrolid, idxpren, dt_index, dt_lookback, enrl_count, 
			min(start) as start_ffsrx format=DATE9., max(end) as end_ffsrx format=DATE9. 
		from _enrollment2
		group by patient_deid, enrolid, idxpren, dt_index, dt_lookback, enrl_count
		;
		quit;

	
	%*Subset to the row that contains the index date, if one contains the index date
	AND
	Calculate the number of days the index date occurs after the start date;
	data _enrollment4;
	set _enrollment3;
		where start_ffsrx <= dt_index <= end_ffsrx;

		prior_enrollment = dt_index - start_ffsrx;
	run;



	/********
	STEP 2b: Add information  related to enrollment in any plan with RX coverage.
	********/

	%*Subset to those continuous enrollment periods for those individuals in our pregnancy dataset via inner join.
	This is using the continuous enrollment derived files created by Virginia Pate. For the exclusoin criteria, 
	we are interested in enrollment in a Fee-for-service (FFS) plan with prescription (rx) coverage.;
	proc sql;
		create table _enrollment_rx as
		select distinct b.patient_deid, a.enrolid, a.start, a.end, b.idxpren, b.dt_index, b.dt_lookback, b.dt_gapreg
		from der.rxenrlper as a
		inner join _pregnancies as b
		on a.enrolid = b.enrolid
		;
		quit;

	%*Output some relevant counts;
	proc sql noprint;
		/*Number of pregnancies that had no enrollment information*/
		select count(distinct idxpren) as num_preg_any_enrl_rxffs
		into :num_preg_any_enrl_rx
		from _enrollment_rx;
		quit;
	%put Number of pregnancies without any Rx enrollment information: &num_preg_any_enrl_rx;

	%*****Implement allowable gap;

	%*Calculate the time elapsed from the last continuous enrollment period;
	proc sort data=_enrollment_rx;
		by patient_deid enrolid idxpren start;
	run;
	data _enrollment2_rx;
	set _enrollment_rx;
		by patient_deid enrolid idxpren start;
		retain enrl_count;
		
		last_id = lag1(idxpren);
		elapsed = start - lag1(end);

		%*Carry down the information through a count variable;
		if first.idxpren then enrl_count = 1;
			else if idxpren = last_id and elapsed > &gap then enrl_count = enrl_count + 1;
	run;

	%*Output the continuous enrollment rows, revised with the gap;
	proc sql;
		create table _enrollment3_rx as
		select distinct patient_deid, enrolid, idxpren, dt_index, dt_lookback, enrl_count, 
			min(start) as start_rx format=DATE9., max(end) as end_rx format=DATE9. 
		from _enrollment2_rx
		group by patient_deid, enrolid, idxpren, dt_index, dt_lookback, enrl_count
		;
		quit;

	
	%*Subset to the row that contains the index date, if one contains the index date
	AND
	Calculate the number of days the index date occurs after the start date;
	data _enrollment4_rx;
	set _enrollment3_rx;
		where start_rx <= dt_index <= end_rx;

		prior_enrollment = dt_index - start_rx;
	run;









	%*link the enrollment information back onto the pregnancies. 
	Further, we want to identify:
		(1) the first disenrollment date from ANY plan after the indexing prenatal encounter, and
		(2) the first disenrollment date from a FFS plan after the indexing prenatal encounter, ;
	proc sql;
		create table _pregnancies2a as
		select distinct a.*, 
			b.start_ffsrx as cont_enrl_start_ffsrx format=DATE9., b.end_ffsrx as cont_enrl_end_ffsrx format=DATE9.,
			bb.start_rx as cont_enrl_start_rx format=DATE9., bb.end_rx as cont_enrl_end_rx format=DATE9.,
			b.prior_enrollment as prior_enrollment_ffsrx, bb.prior_enrollment as prior_enrollment_rx,
			c.end as cont_enrl_end_any format=DATE9., /*End of enrollment from any plan*/
			d.end as cont_enrl_end_ffs format=DATE9. /*End of enrollment from FFS plan*/
		from _pregnancies as a
		left join _enrollment4 as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		left join _enrollment4_rx as bb
		on a.enrolid=bb.enrolid and a.idxpren = bb.idxpren
		left join der.enrlper as c
		on a.enrolid=c.enrolid and c.start <= a.dt_index <= c.end
		left join der.enrlper_ffs as d
		on a.enrolid=d.enrolid and d.start <= a.dt_index <= d.end
		;
		quit;

	%*We want to calcualte the date of the last PNC encounter PRIOR to disenrollment from ANY plan 
	after the indexing prenatal encounter;
	proc sql;
		create table _pregnancies2b as /*Previously _pregnancies2*/
		select distinct a.*, 
			max(case when a.dt_index <= b.enc_date <= a.cont_enrl_end_any then b.enc_date else . end) as cont_enrl_end_pnc_any format=DATE9.,
			max(case when a.dt_index <= b.enc_date <= a.cont_enrl_end_ffs then b.enc_date else . end) as cont_enrl_end_pnc_ffs format=DATE9.,
			max(case when a.dt_index <= b.enc_date <= a.cont_enrl_end_ffsrx then b.enc_date else . end) as cont_enrl_end_pnc_ffsrx format=DATE9.,
			max(case when a.dt_index <= b.enc_date <= a.cont_enrl_end_rx then b.enc_date else . end) as cont_enrl_end_pnc_rx format=DATE9.
		from _pregnancies2a as a
		left join temp.codeprenatal_meg1_dts as b
		on a.patient_deid=b.patient_deid and a.dt_index <= b.enc_date <= a.cont_enrl_end_any
		group by a.enrolid, a.idxpren
		;
		quit;

	%*Finally, we want to know the first date of disenrollment after the outcome date;
	proc sql;
		create table _pregnancies2 as
		select distinct a.*, b.end as dt_disenroll_post_any format=MMDDYY10.,
			c.end as dt_disenroll_post_ffs format=MMDDYY10.,
			d.end as dt_disenroll_post_rx format=MMDDYY10.,
			e.end as dt_disenroll_post_ffsrx format=MMDDYY10.
		from _pregnancies2b as a
		left join der.enrlper as b
		on a.enrolid=b.enrolid and b.start <= a.dt_gapreg <= b.end
		left join der.enrlper_ffs as c
		on a.enrolid=c.enrolid and c.start <= a.dt_gapreg <= c.end
		left join der.rxenrlper as d
		on a.enrolid=d.enrolid and d.start <= a.dt_gapreg <= d.end
		left join der.rxenrlper_ffs as e
		on a.enrolid=e.enrolid and e.start <= a.dt_gapreg <= e.end
		;
		quit;


	%*Delete unnecessary datasets to retain sufficient working memory;
	proc datasets gennum = all noprint;
		delete _enrollment:;
	run;



	/*******************
		STEP 3: Implement the diagnosis-code based exclusion criteria, bulleted below.
		- Chronic hypertension
		- Asthma -- This also has a med component, which we deal with later
		- Preexisting heart disease, excluding congenital heart defects, dealt with later becuase it requires information on provider type.
		- Any prevalent cancer diagnosis.
	******************/

	%*Stack all of the code lists. All included only have diagnosis codes;
	proc sql;
		create table _exclusions as
		select distinct code, "ChHtn" as diagnosis from covref.chronic_htn_dx
		union corr
		select distinct code, "Asthma" as diagnosis from covref.asthma_dx
		union corr
		select distinct code, "CAD" as diagnosis from covref.cad_dx
		union corr
		select distinct code, "Arrhy" as diagnosis from covref.arrhythmia_dx
		union corr
		select distinct code, "Endocard" as diagnosis from covref.endocarditis_dx
		union corr
		select distinct code, "MyoPeri" as diagnosis from covref.myo_pericarditis_dx
		union corr
		select distinct code, "HF" as diagnosis from covref.hf_dx
		union corr
		select distinct code, "HValve" as diagnosis from covref.heart_valve_dx
		union corr
		select distinct code, "CMyo" as diagnosis from covref.cardiomyopathy_dx
		union corr
		select distinct code, "OtherHrt" as diagnosis from covref.other_heart_disease_dx
		union corr
		/*Breast cancer*/
		select distinct code, "Cancer" as diagnosis from covref.breast_cancer_dx
		union corr
		/*Colorectal cancer*/
		select distinct code, "Cancer" as diagnosis from covref.colorectal_cancer_dx
		union corr
		/*Endometrial cancer*/
		select distinct code, "Cancer" as diagnosis from covref.endometrial_cancer_dx
		union corr
		/*Lung cancer*/
		select distinct code, "Cancer" as diagnosis from covref.lung_cancer_dx
		union corr
		/*Urologic cancer*/
		select distinct code, "Cancer" as diagnosis from covref.urologic_cancer_dx
		;
		quit;

	%identify_dx_covariates(REF=_exclusions, PREG_INPT=_pregnancies2, OUTPUT=_pregnancies4a)


	%*Get infomration on chronic hypertension diagnoses prior to 20 weeks gestation;

	%*First, subset to the relevant diagnosis codes;
	data _diagnoses_chtn;
	set _diagnoses2;
		where diagnosis = 'ChHtn';
	run;

	%*Now add on LMP;
	proc sql;
		create table _diagnoses_chtn2 as
		select a.*, b.dt_lmp
		from _diagnoses_chtn as a
		left join _pregnancies4a as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		;
		quit;

	%*Calculate counts for those diagnoses prior to 20w;
	proc sql;
		create table _diagnses_chtn_count as
		select enrolid, idxpren, location, timing, count(distinct svcdate) as numChtn
		from _diagnoses_chtn2 (where = (svcdate < dt_lmp + 140))
		group by enrolid, idxpren, location, timing
		;
		quit;

	%*Transpose the variables wide;
	proc transpose data=_diagnses_chtn_count out=_diagnses_chtn_count_wide (drop = _NAME_) prefix = dx_pre20_CHtn_;
		by enrolid idxpren;
		var numChtn;
		id location timing;
	run;

	%*Join those variables onto the pregnancy dataset;
	proc sql;
		create table _pregnancies4b as
		select a.*, b.*
		from _pregnancies4a as a
		left join _diagnses_chtn_count_wide as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		;
		quit;

	%*All variables with missing values for the dx_ need to be replaced with 0.;
	data _pregnancies4;
	set _pregnancies4b;
		array dx_vars{*} dx_:;
		do i=1 to dim(dx_vars);
			if missing(dx_vars[i]) then dx_vars[i]=0;
		end;
		drop i;
	run;

	/*********
		STEP 4: Get information for asthma medications. These are required to identify individuals with asthma.
	*********/

	%*Create an empty dataset that we are going to append the records too. We need to manually set the length
	and formats so that we avoid errors.;
	data _asthma_meds;
		length enrolid 8 svcdate 8 idxpren 8 dt_index 8 dt_lookback 8 dt_lmp 8 atc_label $56 medication $10 dt_gapreg 8;
		format enrolid best12. svcdate MMDDYY10. idxpren best12. dt_index MMDDYY10. dt_lookback MMDDYY10. dt_lmp MMDDYY10. atc_label $56. medication $10. dt_GApreg MMDDYY10.;
		stop;
	run;

	%*Add medication column to covref.asthma_rx;
	data asthma_rx;
	set rxcov.asthma_rx;
		medication = "Asthma";
	run;
	
	%*Now grab the medications for asthma;
	%get_meds(input_pregs=_pregnancies4, input_meds=asthma_rx, base=_asthma_meds);

	%*Now do some data cleaning - Identify those fills that occur on or prior to the index date and those
	fills that occur within 30 days after the index date;
	data _asthma_meds2;
	set _asthma_meds;
		if svcdate <= dt_index then timing = "Pre";
			else if svcdate > dt_index then timing = "Pos";
			else timing = "";
	run;

	%*Count the number of unique fills that occur within the two periods;
	proc sql;
		create table _asthma_meds3 as
		select distinct enrolid, idxpren, dt_index, timing, count(distinct cats(svcdate, atc_label)) as numRx /*Make sure different ATCs being dispensed if same date*/
		from _asthma_meds2
		group by enrolid, idxpren, dt_index, timing
		;
		quit;

	%*Transpose the dataset so that the distinct periods have two separate columns;
	proc transpose data=_asthma_meds3 out=_asthma_meds4 (drop = _NAME_) prefix = Rx_Asthma;
		by enrolid idxpren;
		var numRx;
		id timing;
	run;

	%*Merge this information back onto the pregnancies dataset;
	proc sql;
		create table _pregnancies5 as 
		select distinct a.*, case
						when b.rx_asthmapre ne . then b.rx_asthmapre
						else 0 end as Rx_AsthmaPre,
					case
						when b.rx_asthmapos ne . then b.rx_asthmapos
						else 0 end as Rx_AsthmaPos
		from _pregnancies4 as a
		left join _asthma_meds4 as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		;
		quit;

	%*Delete the unnecessary datasets;
	proc datasets gennum = all noprint;
		delete _asthma:;
	run;





	/************
		STEP 5: Collect the necessary information to identify individuals with confenital heart defects.

		Algorithm: They need at least 1 diagnosis code where the provider type is Cardiology, Pediatric Cardiology, 
		Cardiovascular Surgery, or Cardiothoracic surgery. Further,this must occur prior to a procedure code for a fetal echocardiogram.
	************/

	%*First step is to grab all procedure codes for fetal echocardiography;

	%*Create an empty file;
	data _congenital_proc;
		length enrolid 8 svcdate 8 idxpren 8 code $7 dt_index 8;
		format enrolid best12. svcdate MMDDYY10. idxpren best12. code $7. dt_index MMDDYY10.;
		stop;
	run;

/*	%let years = 2016 2017 2018 2019 2020 2021 2022;*/
	%get_chd_procedures(input_pregs = _pregnancies5, input_proc = covref.congenital_heart_defect_dx, base=_congenital_proc);


	%*Now grab all diagnosis codes;
	*Create an empty file;
	data _congenital_dx;
		length enrolid 8 svcdate 8 idxpren 8 dt_index 8 dx1 $7 dx2 $7 dx3 $7 dx4 $7;
		format enrolid best12. svcdate MMDDYY10. idxpren best12. dt_index MMDDYY10. dx1 $7. dx2 $7. dx3 $7. dx4 $7.;
		stop;
	run;

	%cardiovascular_diagnoses(input_preg = _pregnancies5, base=_congenital_dx);

	%*The diagnosis codes were in a different format such that each  diagnosis code on a claim is a column.
	Below, we transpose the datset such that all diagnoses are contained in one column: code.
	This will allow us to compare the diagnoses with those from the codelist for maternal congenital heart defects;

	%*Remove duplicates if any exist;
	proc sort data= _congenital_dx nodup;
		by enrolid svcdate idxpren dt_index;
	run;
	%*There are some duplicate diagnosis across the same svcdate for a person, so add this unique row indicator so that 
	we can move from wide to long format with one column;
	data _congenital_dx2;
	set _congenital_dx;
		row = _n_;
	run;
	proc transpose data=_congenital_dx2 out=_congenital_dx_long(rename=(col1=code)) name=source;
	    by row enrolid svcdate idxpren dt_index; /* Variables to keep */
	    var dx1 dx2 dx3 dx4; /* Variables to transpose */
	run;

	%*Remove empty diagnosis codes and the unnecessary viables;
	data _congenital_dx_long2;
	set _congenital_dx_long;
		if code = "" then delete;
		drop row source;
	run;
	%*Now remove duplicate rows;
	proc sort data=_congenital_dx_long2 nodup;
		by enrolid svcdate idxpren dt_index;
	run;

	%*Inner join with the congenital heart defects codes to ensure that
	we only have the relevant diagnosis codes;
	proc sql;
		create table _congenital_dx3 as
		select distinct a.*
		from _congenital_dx_long2 as a
		inner join (select * from covref.congenital_heart_defect_dx where codetype = "DX10") as b
		on a.code = b.code
		;
		quit; /*0 met this in the 1pct sample*/

	%*For each pregnancy, we want their first svcdate with one of these qualifying codes.
	This is because the first congenital heart defect diagnosis must occur PRIOR to the first fetal echocardiography procedure code;
	proc sort data=_congenital_dx3;
		by enrolid idxpren svcdate ;
	run;
	data _congenital_dx4;
	set _congenital_dx3;
		by enrolid idxpren svcdate ;
		if first.idxpren then output;
	run;

	%*Repeat this process for the procedure codes (i.e., grab the date of the first fetal echocardiography procedure code).;
	proc sort data=_congenital_proc;
		by enrolid idxpren svcdate ;
	run;
	data _congenital_proc2;
	set _congenital_proc;
		by enrolid idxpren svcdate ;
		if first.idxpren then output;
	run;

	%*Now, conceptually, we have the first date of the diagnosis code with a congenital heart defect and
	a procedure code for a fetal echocardiography.

	Of note, in the 1pct sample, there were 0 people with the diagnosis codes for congenital heart defects,
	so I could not test that code. However, it is the same set-up as the procedure codes which works 
	as expected.;

	proc sql;
		create table _pregnancies6 as 
		select distinct a.*, b.svcdate as first_fetalecho format=MMDDYY10. length=8,
			c.svcdate as first_congenitalheartdx format=MMDDYY10. length=8
		from _pregnancies5 as a
		left join _congenital_proc2 as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		left join _congenital_dx4 as c
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		;
		quit;

	%*Delete unnecessary datasets;
	proc datasets gennum=all noprint;
		delete _congenital:;
	run;



	/************
		STEP 6: Output the final dataset
	*************/
	data temp.preg_del_incl_excl_&lmpindex._&lookbackdt._&lookbackdays;
	set _pregnancies6;
	run;


	%*Delete the rest of the datasets that output in the macro;
	proc datasets gennum = all noprint;
		delete _:;
	run;


%mend;

























/********************************************************************************************************************************************

														03 - DEFINE COVARIATES

Last key step to deriving the cohort data set involved identifying the key baseline covariates necessary for the analysis.

********************************************************************************************************************************************/


/*
MACRO: combine_codes
PURPOSE: The goal of this program is to combine the codelists from different diagnoses
so that they can easily be put into a macro to identify the count of the codes present.

INPUT:
output_dsn -- Name of the dataset you want output
codetype -- DX10 or PR10 CPT HCPCS
*/
%macro combine_codes(output_dsn, codetype);

	proc sql;
		create table &output_dsn as
		/*secondary hypertension*/
		select distinct code, "SecHtn" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.second_htn_dx where codetype in (&codetype)
		union corr
		/*gestational hypertension*/
		select distinct code, "GHtn" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.gest_htn_dx where codetype in (&codetype)
		union corr
		/*preeclampsia*/
		select distinct code, "Preec" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.preeclampsia_dx where codetype in (&codetype)
		union corr
		/*unspecified hypertension*/
		select distinct code, "UHtn" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.unsp_htn_dx where codetype in (&codetype)
		union corr
		/*HELLP syndrom*/
		select distinct code, "HELLP" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.hellp_dx where codetype in (&codetype)
		union corr
		/*Eclampsia*/
		select distinct code, "Eclamp" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.eclampsia_dx where codetype in (&codetype)
		union corr
		/*opioid use disorder*/
		select distinct code, "Oud" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.oud_dx where codetype in (&codetype)
		union corr
		/*smoking or tobacco use*/
		select distinct code, "Smk" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.smoking_dx where codetype in (&codetype)
		union corr
		/*Alcohol use*/
		select distinct code, "Alc" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.alcohol_dx where codetype in (&codetype)
		union corr
		/*other substance use disorder*/
		select distinct code, "OtherSUD" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.other_sud_dx where codetype in (&codetype)
		union corr
		/*nausea or vomiting*/
		select distinct code, "Nausea" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.nausea_vomiting_dx where codetype in (&codetype) 
		union corr
		/*diabetes - leaving becuase need some of these*/
		select distinct code, "Diabetes" as diagnosis, 0 as exclusion_criteria, outcome, "" as algorithm_step format=$9. length=9
		from covref.diabetes_dx where codetype in (&codetype)
		union corr
		/*migraine*/
		select distinct code, "Migraine" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.migraine_dx where codetype in (&codetype)
		union corr
		/*Recurrent pregnancy loss*/
		select distinct code, "RecurLoss" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.recur_preg_loss_dx where codetype in (&codetype)
		union corr
		/*Chronic kidney disease*/
		select distinct code, "Ckd" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.ckd_dx where codetype in (&codetype) 
		union corr
		/*Obesity*/
		select distinct code, "Obesity" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.obese_dx where codetype in (&codetype) 
		union corr
		/*Diabetic retinopathy*/
		select distinct code, "Retino" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.retinopathy_dx where codetype in (&codetype)
		union corr
		/*Antiphospholipid syndrome - keep but likely will not use these variables*/
		select distinct code, "Antiphos" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, algorithm_step
		from covref.antiphospholipid_dx where codetype in (&codetype)
		union corr
		/*Lupus*/
		select distinct code, "Lupus" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.lupus_dx where codetype in (&codetype)
		union corr
		/*Hyperthyroid disease*/
		select distinct code, "HyperThy" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.hyperthyroid_dx where codetype in (&codetype)
		union corr
		/*Hypothyroid disease*/
		select distinct code, "HypoThy" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.hypothyroid_dx where codetype in (&codetype)
		union corr
		/*Depression*/
		select distinct code, "Depression" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.depression_dx where codetype in (&codetype)
		union corr
		/*Anxiety*/
		select distinct code, "Anxiety" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.anxiety_dx where codetype in (&codetype) 
		union corr
		/*Bipolar disorder*/
		select distinct code, "Bipolar" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.bipolar_disorder_dx where codetype in (&codetype)
		union corr
		/*ADHD*/
		select distinct code, "ADHD" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.adhd_dx where codetype in (&codetype)
		/*PTSD*/
		union corr
		select distinct code, "PTSD" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.ptsd_dx where codetype in (&codetype)
		union corr 
		/*Schizophrenia*/
		select distinct code, "Schizo" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.schizophrenia_dx where codetype in (&codetype)
		union corr
		/*Stroke or TIA*/
		select distinct code, "Stroke" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.stroke_tia_dx where codetype in (&codetype)
		union corr
		/*Myocardial infarction*/
		select distinct code, "MI" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.mi_dx where codetype in (&codetype)
		union corr
		/*Angina*/
		select distinct code, "Angina" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.angina_dx where codetype in (&codetype)
		union corr
		/*Atherosclerosis*/
		select distinct code, "Athero" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.atherosclerosis_dx where codetype in (&codetype)
		union corr
		/*Peripheral vascular disease*/
		select distinct code, "PeriVasc" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.peripheral_vasc_disease_dx where codetype in (&codetype)
		union corr
		/*Hyperlipidemia*/
		select distinct code, "HyperLip" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.hyperlipidemia_dx where codetype in (&codetype)
		union corr
		/*Anemia*/
		select distinct code, "Anemia" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.anemia_dx where codetype in (&codetype)
		union corr
		/*Sickle cell trait*/
		select distinct code, "SickleT" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.sickle_trait_dx where codetype in (&codetype)
		union corr
		/*Sickle cell disease*/
		select distinct code, "SickleD" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.sickle_disease_dx where codetype in (&codetype)
		union corr
		/*Bariatric surgery*/
		select distinct code, "Bariatric" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.bariatric_surgery_dx where codetype in (&codetype)
		union corr
		/*Bariatric surgery revision*/
		select distinct code, "BariatricR" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.rev_bariatric_dx where codetype in (&codetype)
		union corr
		/*CDL: 11.16.2025 ADDED the below code to extract additional diagnoses*/
		/*Anemia*/
		select distinct code, "Anemia" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.anemia_dx where codetype in (&codetype)
		union corr
		/*Gastrointestinal disease*/
		select distinct code, "GIDisease" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.gi_disease_dx where codetype in (&codetype)
		union corr
		/*Pulmonary hypertension*/
		select distinct code, "PulmonaryHTN" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.Pulmonary_htn_dx where codetype in (&codetype)
		union corr
		/*Uterine fibroids*/
		select distinct code, "Fibroids" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.fibroids_dx where codetype in (&codetype)
		;
		quit;

%mend;





/*
MACRO: identify_diabetes_dx
PURPOSE: The purpose of this macro is to identify the additional variables, based on diagnosis codes,
that are necessary to define diabetes types in the dataset.

These include:
- The number of diabetes diagnosis codes (overall and by type) that occur prior to their LMP through LMP +90 
	(within <= 270d prior to or <= 30 d after the index date)
	- We need information by location type (outpt or inpt and overall)
- The number of diabetes diagnosis codes (overall and by type) that occur >= 141d of gestation
- The number of diabetes diagnosis codes (any) <= 140d gestation

INPUTs:
- INPUT_DATA = input pregnancy datset
- OUTPUT_DATA = output pregnancy datset
- DIAGNOSES_DATA = dataset with diagnosis codes from the claims data
*/

%macro identify_diabetes_dx(INPUT_DATA, OUTPUT_DATA, DIAGNOSES_DATA=_diagnoses2);

	/*TESTING:
	%let diagnoses_data = _diagnoses2;
	*/

	%*Import the dataset with all claims with diagnoses from the referent file and subset to those with diabetes diagnosis codes.
		Each row is a diagnosis code.;
	data _diabetes_diagnoses;
	set &DIAGNOSES_DATA;
		where diagnosis = "Diabetes";
	run;

	%*Store a copy of the diabetes diagnoses -- CDL: added 1.30.2025 to ensure variables correct;
	data temp.diabetes_dx_&lmpindex._&lookbackdt._&lookbackdays;
	set _diabetes_diagnoses;
	run;

	*Add on the additional information that we need for diabetes:
		(1) outcome from the diabetes reference file and
		(2) dt_lmp from the pregnancies
	Further, we need additional information from the pregnancies: dt_lookback and dt_gapreg.
	We want to be sure that we dont look at codes after the person disenrolled or had their outcome;
	proc sql;
		create table _diabetes_diagnoses2 as
		select a.*, b.outcome, c.dt_lmp, c.dt_lookback, c.dt_gapreg
		from _diabetes_diagnoses as a
		left join covref.diabetes_dx as b
		on a.code=b.code
		left join &INPUT_DATA as c
		on a.enrolid=c.enrolid and a.idxpren=c.idxpren
		;
		quit;


	%*PREGESTATIONAL
	First, create columns for pregestational diabetes -- We want to know if they have these codes prior to their LMP or <=
	90 days after their LMP. In one analyses, we will require that the 90 day occur prior to the index date and in the other
	we will look up to 30d after the LMP.
	This code collects these counts by location of claim (inpt or outpt);
	proc sql;

		%*Hard stop at index;
		create table _pregestational_index as
		select enrolid, idxpren, location, count(distinct svcdate) as numdx
		from (select * from _diabetes_diagnoses2 where dt_lookback <= svcdate <= min(dt_lmp+90, dt_index)) as a
		group by enrolid, idxpren, location
		;

		%*30-day allotment - include end of enrollment;
		create table _pregestational_index30 as
		select enrolid, idxpren, location, count(distinct svcdate) as numdx
		from (select * from _diabetes_diagnoses2 where dt_lookback <= svcdate <= min(dt_lmp+90, min(dt_index+30, dt_gapreg))) as a
		group by enrolid, idxpren, location
		;

		quit;

	%*Now add these variables onto the same dataset and transpose it;
	proc sql;
		create table _pregestational as
		select enrolid, idxpren, location, numdx, "Idx" as timing
		from _pregestational_index
		union corr
		select enrolid, idxpren, location, numdx, "Idx30" as timing
		from _pregestational_index30
		;
		quit;
	proc transpose data=_pregestational out=_pregestational_wide (drop = _NAME_) prefix=dx_DiabPreG;
		by enrolid idxpren;
		var numdx;
		id location timing;
	run;

	%*Repeat the above steps but now include the variable outcome in the group by statement. 
	Ignore the location of the diagnosis;
	proc sql;

		%*Hard stop at indx;
		create table _pregestational_type_index as
		select enrolid, idxpren, outcome, count(distinct svcdate) as numdx
		from (select * from _diabetes_diagnoses2 where dt_lookback <= svcdate <= min(dt_lmp+90, dt_index)) as a
		group by enrolid, idxpren, outcome
		;

		%*30-day allotment - include end of enrollment;
		create table _pregestational_type_index30 as
		select enrolid, idxpren, outcome, count(distinct svcdate) as numdx
		from (select * from _diabetes_diagnoses2 where dt_lookback <= svcdate <= min(dt_lmp+90, min(dt_index+30, dt_gapreg))) as a
		group by enrolid, idxpren, outcome
		;

		quit;

	%*Now add these variables onto the same dataset;
	proc sql;
		create table _pregestational_type as
		select enrolid, idxpren, outcome, numdx, "Idx" as timing
		from _pregestational_type_index
		union corr
		select enrolid, idxpren, outcome, numdx, "Idx30" as timing
		from _pregestational_type_index30
		;
		quit;

	%*Create a modified outcome variable that we can use in a transpose statement;
	data _pregestational_type;
	set _pregestational_type;
		if outcome = "Pregestational diabetes mellitus, Not otherwise specified" then outc = "PDMNOS";
			else if outcome = "Diabetes mellitus, Not otherwise specified" then outc = "DMNOS";
			else if outcome = "Gestational diabetes" then outc = "GDM";
			else if outcome = "Type-1 diabetes mellitus" then outc = "T1DM";
			else if outcome = "Type-2 diabetes mellitus" then outc = "T2DM";
	run;

	%*Now, transpose the dataset;
	proc transpose data=_pregestational_type out=_pregestational_type_wide (drop = _NAME_) prefix=dx_DiabPreG;
		by enrolid idxpren;
		var numdx;
		id outc timing;
	run;

	%*Add these variables onto the pregnancy dataset;
	proc sql;
		create table _pregnancies_dx2 as
		select a.*, b.*, c.* /*Intentionally have a warning on idxpren and enrolid*/
		from &INPUT_DATA as a
		left join _pregestational_wide as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		left join _pregestational_type_wide as c
		on a.enrolid=c.enrolid and a.idxpren=c.idxpren
		;
		quit;

	%**********GESTATIONAL
	Now we create columns for diabetes codes necessary for identifying gestational diabetes;

	proc sql;

		%*Some codes we want to know >= 141d gestation, hard stop at index;
		create table _late_gestation_index as
		select enrolid, idxpren, outcome, count(distinct svcdate) as numdx
		from (select * 
				from _diabetes_diagnoses2 
				where dt_lmp + 141 <= dt_index and dt_lmp + 141 <= svcdate <= dt_index)
		group by enrolid, idxpren, outcome
		;

		%*Some codes we want to know >= 141d gestation, with 30 d period after index;
		create table _late_gestation_index30 as
		select enrolid, idxpren, outcome, count(distinct svcdate) as numdx
		from (select * 
				from _diabetes_diagnoses2 
				where dt_lmp + 141 <= min(dt_index + 30, dt_gapreg) and 
					dt_lmp + 141 <= svcdate <= min(dt_index + 30, dt_gapreg))
		group by enrolid, idxpren, outcome
		;

		%*Some codes we want to know <= 140 d gestation, hard stop at index;
		create table _early_gestation_index as
		select enrolid, idxpren, count(distinct svcdate) as numdx
		from (select * 
				from _diabetes_diagnoses2 
				where dt_lookback <= svcdate <= min(dt_index, dt_lmp + 140))
		group by enrolid, idxpren
		;

		%*Some codes we want to know <= 140 d gestation, with 30-d grace period after index;
		create table _early_gestation_index30 as
		select enrolid, idxpren, count(distinct svcdate) as numdx
		from (select * 
				from _diabetes_diagnoses2 
				where dt_lookback <= svcdate <= min(dt_lmp + 140, min(dt_index+30, dt_GApreg)))
		group by enrolid, idxpren
		;
		quit;

	%*Union together all of the numbers that we want;
	proc sql;
		create table _gestational as
		/*Number of diabetes codes between LMP+141 and index date*/
		select enrolid, idxpren, "ge141" as timing, "idx" as idx, "DM" as outcome, sum(numdx) as numdx 
		from _late_gestation_index
		group by idxpren
		union corr
		/*Number of gestational diabetes codes between LMP+141 and index date*/
		select enrolid, idxpren, "ge141" as timing, "idx" as idx, numdx, "GDM" as outcome
		from _late_gestation_index
		where outcome = "Gestational diabetes"
		union corr
		/*Number of diabetes codes between LMP+141 and min(index date+30, end of enrollment)*/
		select enrolid, idxpren, "ge141" as timing, "idx30" as idx, "DM" as outcome, sum(numdx) as numdx 
		from _late_gestation_index30
		group by idxpren
		union corr
		/*Number of gestational diabetes codes between LMP+141 and min(index date+30, end of enrollment)*/
		select enrolid, idxpren, "ge141" as timing, "idx30" as idx, numdx, "GDM" as outcome
		from _late_gestation_index30
		where outcome = "Gestational diabetes"
		union corr
		/*Number of diabetes codes between lookback date and min(LMP+140, index date)*/
		select enrolid, idxpren, "le140" as timing, "idx" as idx, "DM" as outcome, sum(numdx) as numdx
		from _early_gestation_index
		group by idxpren
		union corr
		/*Number of diabetes codes between lookback date and min(LMP+140, indexdate + 30, end of enrollment)*/
		select enrolid, idxpren, "le140" as timing, "idx30" as idx, "DM" as outcome, sum(numdx) as numdx
		from _early_gestation_index30
		group by idxpren
		;
		quit;

	%*Now, transpose the dataset;
	proc transpose data=_gestational out=_gestational_wide (drop = _NAME_) prefix=dx_Diab;
		by enrolid idxpren;
		var numdx;
		id outcome timing idx;
	run;

	%*Add the columns to the pregnancy dataset;
	proc sql;
		create table _pregnancies_dx3 as
		select a.*, b.* /*Intentionally have a warning on idxpren and enrolid*/
		from _pregnancies_dx2 as a
		left join _gestational_wide as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		;
		quit;

	%*All variables with missing values for the dx_ need to be replaced with 0.;
	data &OUTPUT_DATA;
	set _pregnancies_dx3;
		array dx_vars{*} dx_:;
		do i=1 to dim(dx_vars);
			if missing(dx_vars[i]) then dx_vars[i]=0;
		end;
		drop i;
	run;

	%*Delete the unnecessary datasets;
	proc datasets gennum=all noprint;
		delete _diabetes: _pregestational: _pregnancies_dx2
			_pregnancies_dx3 _late: _early: _gestational: ;
	run; quit; run;

%mend;







/*
MACRO: identify_htn_dx
PURPOSE: The purpose of this macro is to identify the additional variables, based on diagnosis codes,
that are necessary to define hypertension types in the dataset.

INPUTs:
- INPUT_DATA = input pregnancy datset
- OUTPUT_DATA = output pregnancy datset
- DIAGNOSES_DATA = dataset with diagnosis codes from the claims data
*/

%macro identify_htn_dx(INPUT_DATA, OUTPUT_DATA, DIAGNOSES_DATA=_diagnoses2);

	/*TESTING:
	%let diagnoses_data = _diagnoses2;
	%let input_data = _pregnancies_diab; 
	*/

	%*Import the dataset with all claims with diagnoses from the referent file and subset to those with diabetes diagnosis codes.
		Each row is a diagnosis code.;
	data _htn_diagnoses;
	set &DIAGNOSES_DATA;
		where diagnosis in ("GHtn", "Preec", "UHtn", "HELLP", "Eclamp");
	run;


	*Add on the additional pregnancy information that we need for HTN:
	We want to be sure that we dont look at codes after the person had their outcome;
	proc sql;
		create table _htn_diagnoses2 as
		select a.*, c.dt_lmp, c.dt_lookback, c.dt_gapreg
		from _htn_diagnoses as a
		left join &INPUT_DATA as c
		on a.enrolid=c.enrolid and a.idxpren=c.idxpren
		;
		quit;

	%*We want counts between 0, 4, 15, and 20 weeks gestation and the index date (or index date + 30);

	%let lmps = 0 4 15 20;
	%let numlmp = %sysfunc(countw(&lmps));

	%do d=1 %to &numlmp;

		%let loop&d = %scan(&lmps, &d);

		%*Count the number of inpatient and outpatient encounters;
		proc sql;
			create table _htn_diagnoses_lmp&&loop&d as
			select enrolid, idxpren, diagnosis, location, timing,
				count(distinct svcdate) as lmp&&loop&d
			from _htn_diagnoses2 where dt_lmp + (7 * &&loop&d) <= svcdate
			group by enrolid, idxpren, diagnosis, location, timing
			;
			quit;

		%*Transpose the dataset - all counts;
		proc transpose data=_htn_diagnoses_lmp&&loop&d out=_htn_diagnoses_count_wide&&loop&d (drop = _NAME_) prefix=dx_lmp&&loop&d;
			by enrolid idxpren;
			var lmp&&loop&d;
			id diagnosis location timing;
		run;

		%*Count the number of inpatient admissions with dx in primary or secondary position;
		proc sql;
			create table _htn_diagnses_ps_&&loop&d as
			select enrolid, idxpren, diagnosis, location, timing,
				count(distinct svcdate) as lmp&&loop&d
			from _htn_diagnoses2 where dxloc = "InptAdm" and dxnum <= 2 and dt_lmp + (7 * &&loop&d) <= svcdate
			group by enrolid, idxpren, diagnosis, location, timing
			;
			quit;

		%*Transpose the dataset - all counts;
		proc transpose data=_htn_diagnses_ps_&&loop&d out=_htn_diagnoses_count_wide_ps&&loop&d (drop = _NAME_) prefix=dx_lmp&&loop&d.._ps;
			by enrolid idxpren;
			var lmp&&loop&d;
			id diagnosis location timing;
		run;

	%end;

	
	%*Add the columns to the pregnancy dataset;
	proc sql;
		create table _pregnancies_dxlmp as
		select distinct a.*, b.*, c.*, d.*, e.*, f.*, g.*, h.*, i.* /*Intentionally have a warning on idxpren and enrolid*/
		from &input_data as a
		left join _htn_diagnoses_count_wide0 as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		left join _htn_diagnoses_count_wide_ps0 as c
		on a.enrolid=c.enrolid and a.idxpren=c.idxpren
		left join _htn_diagnoses_count_wide4 as d
		on a.enrolid=d.enrolid and a.idxpren=d.idxpren
		left join _htn_diagnoses_count_wide_ps4 as e
		on a.enrolid=e.enrolid and a.idxpren=e.idxpren
		left join _htn_diagnoses_count_wide15 as f
		on a.enrolid=f.enrolid and a.idxpren=f.idxpren
		left join _htn_diagnoses_count_wide_ps15 as g
		on a.enrolid=g.enrolid and a.idxpren=g.idxpren
		left join _htn_diagnoses_count_wide20 as h
		on a.enrolid=h.enrolid and a.idxpren=h.idxpren
		left join _htn_diagnoses_count_wide_ps20 as i
		on a.enrolid=i.enrolid and a.idxpren=i.idxpren
		;
		quit;

	%*All variables with missing values for the dx_ need to be replaced with 0.;
	data &OUTPUT_DATA;
	set _pregnancies_dxlmp;
		array dx_vars{*} dx_:;
		do i=1 to dim(dx_vars);
			if missing(dx_vars[i]) then dx_vars[i]=0;
		end;
		drop i;
	run;

	%*Delete the unnecessary datasets;
	proc datasets gennum=all noprint;
		delete _diabetes: _pregestational: _pregnancies_dx2
			_pregnancies_dx3 _late: _early: _gestational: ;
	run; quit; run;

%mend;







/*
MACRO: get_procedures
PURPOSE: The purpose of this macro is to grab all procedure codes for a reference file.

INPUTS:
INPUT_PREGS -- Input dataset where each row is a pregnancy
INPUT_PROC -- Input reference dataset for procedure codes
BASE -- Dataset on which to append the procedure codes
*/
%macro get_procedures(input_pregs, input_proc, base);

	%let numYr = %sysfunc(countw(&years));

	%do d=1 %to &numYr;

		%let loop&d = %scan(&years, &d);

		proc sql;
			create table _proc_subset as
			select distinct a.enrolid as enrolid format=best12. length=8, a.svcdate as svcdate format=MMDDYY10. length=8,
					b.idxpren as idxpren format=best12. length=8, a.proc&&loop&d as code format=$7. length=7,
					b.dt_index as dt_index format=MMDDYY10. length=8, b.dt_lookback as dt_lookback format=MMDDYY10. length=8,
					b.dt_lmp as dt_lmp format=MMDDYY10. length=8,
					b.dt_gapreg as dt_gapreg format=MMDDYY10. length=8,
					c.diagnosis format=$10. length=10
			from der.allproc&&loop&d as a
			inner join &input_pregs as b
			on a.enrolid=b.enrolid and b.dt_lookback <= a.svcdate and a.svcdate <= min(b.dt_index + 30, b.dt_gapreg)
			inner join &input_proc as c
			on a.proc&&loop&d = c.code
			;
			quit;

		proc append base=&base data=_proc_subset; run;

	%end;

%mend;


/*
MACRO: identify_pr_covariates
PURPOSE: The purpose of this macro is to output necessary covariates using procedure codes. This includes the complex
approach for diabetes.

INPUTS:
INPUT_DATA -- Input pregnancy-level dataset
OUTPUT_DATA -- Output dataset name
PROC_REF -- Reference file for the procedure codes
*/

%macro identify_pr_covariates(INPUT_DATA, OUTPUT_DATA, PROC_REF);

	/*TESTING:
	%let input_data = _pregnancies_htn;
	%let output_data = _pregnancies_pr;
	%let proc_ref = cov_procedures;
	*/

	%*Grab all the procedure codes between the lookback date and 30 d after the index date;

	%*Create an empty file;
	data _covariate_proc;
		length enrolid 8 svcdate 8 idxpren 8 code $7 dt_index 8 diagnosis $10 dt_lookback 8 dt_lmp 8 dt_gapreg 8;
		format enrolid best12. svcdate MMDDYY10. idxpren best12. code $7. dt_index MMDDYY10. diagnosis $10. dt_lookback MMDDYY10. dt_lmp MMDDYY10. dt_gapreg MMDDYY10.;
		stop;
	run;

	%get_procedures(input_pregs=&INPUT_DATA, input_proc=cov_procedures, base=_covariate_proc);

	%*CDL: 1.28.2025 added to ensure that there are no duplicate rows, though there should not be;
	proc sort data=_covariate_proc nodup;
		by enrolid idxpren svcdate;
	run;


	%*STEP 1: Collect simple counts around the index and lookback dates;

	%*Create an indicator for pre versus post index;
	data _covariate_proc2;
	set _covariate_proc;
		if svcdate <= dt_index then timing = "Pre";
			else timing = "Pos";
	run;

	%*First, identify counts before and after the index date;
	proc sql;
		create table _proc_counts as 
		select distinct enrolid, idxpren, diagnosis, timing, count(distinct svcdate) as numPr /*CDL: 1.28.2025 added distinct */
		from _covariate_proc2 
		group by enrolid, idxpren, diagnosis, timing
		;
		quit;

	%*Transpose the dataset so that each count is a column;
	proc transpose data=_proc_counts out=_proc_counts_wide (drop = _NAME_) prefix=pr_;
		by enrolid idxpren;
		var numpr;
		id diagnosis timing;
	run;

	%*Add these variables back onto the pregnancy dataset;
	proc sql;
		create table _pregnancies_proc as
		select a.*, b.* /*This will throw a warning on idxpren and enrolid but that is expected and okay*/
		from &INPUT_DATA as a
		left join _proc_counts_wide as b
		on a.enrolid=b.enrolid and a.idxpren = b.idxpren
		;
		quit;




	%*STEP 2: Get the covariate information required to identify diabetes based upon oral glucose tolerance
	tests;


	%*STEP 2a: The date of the first OGTT between LMP+91 and LMP +140 and the number of DM diagnosis codes 
	after this date.;

	%*Get this information from the the list of procedures for Diabetes;
	proc sql;
		
		%*Only consider information up to the index date;
		create table early_gtt0 as 
		/*CDL: 1.28.2025 added distinct */
		select distinct enrolid, idxpren, "91to140" as ga, "Idx" as lkfwd_dt, min(svcdate) as dt_gtt91to140 format=MMDDYY10.,
			min(dt_lmp+140, dt_index) as dt_lookforward format = MMDDYY10.
		from _covariate_proc2 (where = (diagnosis = "Diabetes" and
										dt_index >= dt_lmp + 91 and
										dt_lmp+91 <= svcdate <= min(dt_lmp+140, dt_index)))
		group by enrolid, idxpren
		; 

		%*Consider information up to 30 days after the index date;
		create table early_gtt30 as
		/*CDL: 1.28.2025 added distinct */
		select distinct enrolid, idxpren, "91to140" as ga, "Idx30" as lkfwd_dt, min(svcdate) as dt_gtt91to140 format=MMDDYY10.,
			min(dt_lmp+140, min(dt_index+30, dt_gapreg)) as dt_lookforward format=MMDDYY10.
		from _covariate_proc2 (where = (diagnosis = "Diabetes" and
										dt_index >= dt_lmp + 91 and
										dt_lmp+91 <= svcdate <= min(dt_lmp +140, min(dt_index+30, dt_gapreg))))
		group by enrolid, idxpren
		; 

		quit;

	%*Stack the datasets;
	proc sql;
		create table early_gtt as
		select * from early_gtt0
		union /*all*/ /*CDL: 1.28.2025 removed the all operator. Ensures no duplicate rows from the two datasets*/
		select * from early_gtt30
		;
		quit;

	%*Now, we want to know the number of diagnosis codes that occurred between that GTT date and the index or index +30;
	proc sql;
		create table _gtt_dx as
		select distinct a.*, count(distinct b.svcdate) as numDx /*CDL: 1.29.2025 - removed b.timing. Not necessary and created errors*/
		from early_gtt as a
		left join _diagnoses2 (where = (diagnosis = "Diabetes")) as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren and a.dt_gtt91to140 <= b.svcdate <= a.dt_lookforward 
			/*CDL: 1.28.25 - Added a. to date variables. 1.29.2025 - Added b. to svcdate*/
		group by a.enrolid, a.idxpren, a.ga, a.lkfwd_dt, a.dt_gtt91to140, a.dt_lookforward
		;
		quit; 

	%*If the date of the OGTT is not defined, then numDx should be 0 - unnecessary step;
	data _gtt_dx;
	set _gtt_dx;
		if dt_gtt91to140 = . then numDx = 0;
			else numdx = numDx;
	run;

	%*Finally, transpose the dataset;

	%*First, transpose the dates to their wide version.;
	proc sort data=early_gtt nodup; by enrolid idxpren; run; /*CDL: 1.28.2025 added nodupkey */
	proc transpose data=early_gtt out=_dt_gtt_wide (drop = _NAME_) prefix=dt_gtt;
		by enrolid idxpren;
		var dt_gtt91to140;
		id ga lkfwd_dt;
	run;

	%*Now transpose the diagnosis codes;
	proc sort data=_gtt_dx nodup; by enrolid idxpren; run; /*CDL: 1.28.2025 added nodupkey*/
	proc transpose data=_gtt_dx out=_dx_gtt_wide (drop = _NAME_) prefix=dx_gtt;
		by enrolid idxpren;
		var numdx;
		id ga lkfwd_dt;
	run;

	%*Now make the dataset with all the information;
	proc sql;
		create table _gtt_dx_wide as 
		select a.*, b.* /*provides a warning, but that is okay*/
		from _dt_gtt_wide as a
		left join _dx_gtt_wide as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		;
		quit;




	%*STEP 2B: Calculate the number of oral glucose tolerance tests after LMP +141
	;
	proc sql;

		%*Hard stop at index date;
		create table _late_ogtt0 as
		/*CDL: 1.28.2025 - added distinct*/
		select distinct enrolid, idxpren, "ge141" as ga, "Idx" as timing, count(distinct svcdate) as numPr
		from _covariate_proc2 
		where diagnosis = "Diabetes" and dt_lmp+141 <= dt_index and
			dt_lmp+141 <= svcdate <= dt_index
		group by enrolid, idxpren
		;

		%*Allow 30 days after index date or end of continuous enrollment;
		create table _late_ogtt30 as
		/*CDL: 1.28.2025 - added distinct*/
		select distinct enrolid, idxpren, "ge141" as ga, "Idx30" as timing, count(distinct svcdate) as numPr
		from _covariate_proc2 
		where diagnosis = "Diabetes" and dt_lmp+141 <= dt_index and
			dt_lmp+141 <= svcdate <= min(dt_index+30, dt_gapreg)
		group by enrolid, idxpren
		;
		quit;

	%*Stack the datasets;
	proc sql;
		create table _late_gtt as
		select * from _late_ogtt0
		union
		select * from _late_ogtt30
		;
		quit;

	%*Transpose the dataset so that each is a column;
	proc transpose data=_late_gtt out=_late_gtt_wide (drop = _NAME_) prefix=pr_gttge141;
		by enrolid idxpren;
		var numpr;
		id timing;
	run;


	%*Put all the information onto the pregnancies dataset.;

	proc sql;
		create table _pregnancies_proc_diab as
		/*CDL: 1.28.2025 - added distinct*/
		select distinct a.*, b.dt_gtt91to140Idx, b.dt_gtt91to140Idx30, b.dx_gtt91to140Idx, b.dx_gtt91to140Idx30,
			c.pr_gttge141Idx, c.pr_gttge141Idx30
		from _pregnancies_proc as a
		left join _gtt_dx_wide as b
		on a.enrolid = b.enrolid and a.idxpren = b.idxpren
		left join _late_gtt_wide as c
		on a.enrolid = c.enrolid and a.idxpren = c.idxpren
		;
		quit;


	%*All variables with missing values for the dx_ or pr_ need to be replaced with 0.;
	data &OUTPUT_DATA;
	set _pregnancies_proc_diab;

		%*DX variables;
		array dx_vars{*} dx_:;
		do i=1 to dim(dx_vars);
			if missing(dx_vars[i]) then dx_vars[i]=0;
		end;
		
		%*PR variables;
		array pr_vars{*} pr_:;
		do j=1 to dim(pr_vars);
			if missing(pr_vars[j]) then pr_vars[j]=0;
		end;
		
		drop i j;
	run;

	%*All variables with missing values for the pr_ need to be replaces with 0.;

	%*Delete the unnecessary datasets;
	proc datasets gennum=all noprint;
		delete _covariate: _proc: early: _gtt: _late: _pregnancies_proc_diab _pregnancies_proc;
	run; quit; run;

%mend;




/*
MACRO: identify_rx_covariates
PURPOSE: The purpose of this macro is to identify covariates that are defined using prescription fills

INPUTS:
- INPUT_DATA: input dataset of the pregnancies to reference against
- OUTPUT_DATA: name for the output pregnancy dataset
- REF_DATA: Theh reference dataset for the medications
*/

%macro identify_rx_covariates(INPUT_DATA, OUTPUT_DATA, REF_DATA);

	%*Now get all the medications between dt_lookback and 30 days after index;
	%*Create an empty dataset that we are going to append the records too. We need to manually set the length
	and formats so that we avoid errors.;
	data _cov_meds;
		length enrolid 8 svcdate 8 idxpren 8 dt_index 8 dt_lookback 8 dt_lmp 8 atc_label $56 medication $10 dt_gapreg 8;
		format enrolid best12. svcdate MMDDYY10. idxpren best12. dt_index MMDDYY10. dt_lookback MMDDYY10. dt_lmp MMDDYY10. atc_label $56. medication $10. dt_gapreg MMDDYY10.;
		stop;
	run;
		
	%*Now grab the medications for asthma;
	%get_meds(input_pregs=&INPUT_DATA, input_meds= &REF_DATA, base=_cov_meds);

	%*Make sure that there are no duplicates - there should not be; /*CDL: 1.28.2025 - added this sort to be sure that there are no duplicates*/
	proc sort data=_cov_meds nodup; by enrolid svcdate; run;

	%*Output a version of this dataset that we can review later -- T2DM meds are unexpected and want to investigate more deeply;
	data temp.cov_meds_&lmpindex._&lookbackdt._&lookbackdays; set _cov_meds; run;

	%*Add an indicator for if the medication fill occurred before or after index;
	data _cov_meds2;
	set _cov_meds;
		if svcdate <= dt_index then timing = "Pre";
			else timing = "Pos";
	run;

	%*I want counts of the number of fills prior to the index date;
	proc sql;
		create table _cov_meds_count as
		/*CDL: 1.28.2025 - added distinct*/
		select distinct enrolid, idxpren, medication, timing, count(*) as numRx
		from _cov_meds2
		group by enrolid, idxpren, medication, timing
		;
		quit;

	%*Transpose the dataset so that we have columns for each of the timings;
	proc sort data=_cov_meds_count nodup; by enrolid idxpren; run; /*CDL: 1.28.2025 - added sort to ensure no duplicate observations -  there should not be*/
	proc transpose data=_cov_meds_count out=_cov_meds_count_wide (drop = _NAME_) prefix=rx_;
		by enrolid idxpren;
		var numRx;
		id medication timing;
	run;

	%*Add these columns back to the pregnancies dataset;
	proc sql;
		create table _pregnancies_rx1 as
		select a.*, b.* /*This will throw a warning for idxpren and enrolid being on both datasets but that is okay*/
		from &INPUT_DATA as a
		left join _cov_meds_count_wide as b
		on a.enrolid=b.enrolid and a.idxpren = b.idxpren
		;
		quit;


	%****Now we want information for non-metformin antidiabetics, indexing around LMP;

	%*First, grab all the fills for non-metformin antidiabetics;
	data _cov_meds_diabetes;
	set _cov_meds2;
		where medication in ("t1t2dm" "t2dm");
	run;

	%*Calculate counts indexed on LMP;
	proc sql;

		%****COUNTS DT_LOOKBACK TO LMP+90, up to index date;

		%*Hard stop at index date;
		create table _cov_meds_diab_count_early0 as
		/*CDL: 1.28.2025 - added distinct*/
		select distinct enrolid, idxpren, "le90GA" as ga, "idx" as timing, sum(case when upcase(atc_label) contains 'INSULIN' then 1 else 0 end) as sumInsulin,
			count(*) as numRx, sum(medication = "t2dm") as numT2DMRx
		from _cov_meds_diabetes
		where dt_lookback <= svcdate <= min(dt_lmp + 90, dt_index)
		group by enrolid, idxpren
		;

		%*Allow 30 days after index date;
		create table _cov_meds_diab_count_early30 as
		/*CDL: 1.28.2025 - added distinct*/
		select distinct enrolid, idxpren, "le90GA" as ga, "idx30" as timing, sum(case when upcase(atc_label) contains 'INSULIN' then 1 else 0 end) as sumInsulin,
			count(*) as numRx, sum(medication = "t2dm") as numT2DMRx
		from _cov_meds_diabetes
		where dt_lookback <= svcdate <= min(dt_lmp + 90, min(dt_index + 30, dt_gapreg))
		group by enrolid, idxpren
		;

		%*****COUNTS DT_LOOKBACK TO LMP+140, up to index date;
		%*Hard stop at index date;
		create table _cov_meds_diab_count_late0 as
		/*CDL: 1.28.2025 - added distinct*/
		select distinct enrolid, idxpren, "le140GA" as ga, "idx" as timing, sum(case when upcase(atc_label) contains 'INSULIN' then 1 else 0 end) as sumInsulin,
			count(*) as numRx, sum(medication = "t2dm") as numT2DMRx
		from _cov_meds_diabetes
		where dt_lookback <= svcdate <= min(dt_lmp + 140, dt_index)
		group by enrolid, idxpren
		;

		%*Allow 30 days after index date;
		create table _cov_meds_diab_count_late30 as
		/*CDL: 1.28.2025 - added distinct*/
		select distinct enrolid, idxpren, "le140GA" as ga, "idx30" as timing, sum(case when upcase(atc_label) contains 'INSULIN' then 1 else 0 end) as sumInsulin,
			count(*) as numRx, sum(medication = "t2dm") as numT2DMRx
		from _cov_meds_diabetes
		where dt_lookback <= svcdate <= min(dt_lmp + 140, min(dt_index + 30, dt_gapreg))
		group by enrolid, idxpren
		;
		quit;


	%*Union the datasets;
	proc sql;
		create table _meds_diabetes_all as
		select * from _cov_meds_diab_count_early0
		union corr
		select * from _cov_meds_diab_count_early30
		union corr
		select * from _cov_meds_diab_count_late0
		union corr
		select * from _cov_meds_diab_count_late30
		;
		quit;

	%*Transpose the dataset just for the insulin counts;
	proc sort data=_meds_diabetes_all nodup; by enrolid idxpren; run; /*CDL: 1.28.2025 added to ensure no duplicates though there should not be any*/
	proc transpose data=_meds_diabetes_all out=_meds_insulin (drop = _NAME_) prefix=rx_insulin;
		by enrolid idxpren;
		var sumInsulin;
		id ga timing;
	run;

	%*Transpose the dataset just for the non-metformin antidiabetics counts;
	proc transpose data=_meds_diabetes_all out=_meds_nonmet_antidiab(drop = _NAME_) prefix=rx_nonmet_antidiab;
		by enrolid idxpren;
		var numRx;
		id ga timing;
	run;

	%*Transpose the dataset just for the non-metformin T2DM antidiabetics counts;
	proc transpose data=_meds_diabetes_all out=_meds_t2dm_antidiab(drop = _NAME_) prefix=rx_t2dm_antidiab;
		by enrolid idxpren;
		var numT2DMRx;
		id ga timing;
	run;

	%*Finally, add these variables onto the pregnancies dataset.;
	proc sql;
		create table _pregnancies_rx as
		/*CDL: 1.28.2025 - added distinct*/
		select distinct a.*, b.*, c.*, d.* /*This will lead to a warning on enrolid and idxpren being present in both data sets but that is okay*/
		from _pregnancies_rx1 as a
		left join _meds_insulin as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		left join _meds_nonmet_antidiab as c
		on a.enrolid=c.enrolid and a.idxpren=c.idxpren
		left join _meds_t2dm_antidiab as d
		on a.enrolid=d.enrolid and a.idxpren=d.idxpren
		;
		quit;

	%*All variables with missing values for the rx_ need to be replaced with 0.;
	data &OUTPUT_DATA;
	set _pregnancies_rx;

		%*RX variables;
		array rx_vars{*} rx_:;
		do i=1 to dim(rx_vars);
			if missing(rx_vars[i]) then rx_vars[i]=0;
		end;
			
		drop i;
	run;

	%*Delete the unnecessary datasets;
	proc datasets gennum=all noprint;
		delete _cov_meds: _pregnancies_rx: _meds_diabetes: _meds: ;
	run; quit; run;

%mend;








/*
MACRO: get_inpt_adm
PURPOSE: To identify distinct inpatient admission dates

INPUT:
INPUT_PREGS -- Input dataset with all the pregnancies that we are identifying information for
BASE -- Dataset on which we will append inpatient admission dates.
*/

%*Grab all inpatient admission dates;
%macro get_inpt_adm(input_pregs, base);

	%let numYr = %sysfunc(countw(&years));

	%do d=1 %to &numYr;

		%let loop&d = %scan(&years, &d);

		proc sql;
			create table _adm_subset as
			select distinct a.enrolid as enrolid format=best12. length=8, a.admdate as admdate format=MMDDYY10. length=4
			from raw.inptadm&&loop&d as a
			inner join &input_pregs as b
			on a.enrolid=b.enrolid and b.dt_lookback <= a.admdate <= b.dt_index
			;
			quit;

		proc append base=&base data=_adm_subset; run;

	%end;

%mend;





/*
MACRO: get_healthcare_use
PURPOSE: To identify the number of pre-index outpatient and inpatient prenatal encounters and inpatient admissions

INPUTS:
INPUT_DATA -- Input dataset of pregnancies that want to analyze
OUTPUT_DATA -- Output dataset of pregnancies
*/

%macro get_healthcare_use(INPUT_DATA, OUTPUT_DATA);

	%*First, identify the number of distinct prenatal encounter dates prior to the index date;
	proc sql;
		create table _pregnancies_pnc1 as
		select distinct a.*, sum(b.pren_outpatient)as num_OutptPNC, sum(b.pren_inpatient) as num_InptPNC
		from &INPUT_DATA as a
		left join (select distinct * from temp.codeprenatal_meg1_dts) as b
		on a.patient_deid=b.patient_deid and a.dt_indexprenatal <= b.enc_date <= a.dt_index
		group by a.patient_deid, a.idxpren
		;
		quit;

	%*Second, identify the number of distinct encounter dates with a Specific GA code prior to the index date.
	- Identify those between the indexing prenatal encounter (dt_indexprenatal) and the index date (dt_index).;
	proc sql;
		create table _pregnancies_pnc2 as
		select distinct a.*, count(distinct c.enc_date) as num_SpecificGA_preidx
		from _pregnancies_pnc1 as a
		left join (select distinct * from temp.gestagepren (where = (code_hierarchy = "Specific gestational age"))) as c
		on a.patient_deid=c.patient_deid and a.dt_indexprenatal <= c.enc_date <= a.dt_index
		group by a.patient_deid, a.idxpren
		;
		quit;

	%*Third, identify the number of specific gestational age codes associated with the pregnancy.;
	proc sql;
		create table _pregnancies_pnc as
		select distinct a.*, count(distinct c.enc_date) as num_SpecificGA_preg
		from _pregnancies_pnc2 as a
		left join (select distinct * from temp.gestagepren (where = (code_hierarchy = "Specific gestational age"))) as c
		on a.patient_deid=c.patient_deid and a.dt_indexprenatal <= c.enc_date <= a.dt_gapreg + 7
		group by a.patient_deid, a.idxpren
		;
		quit;

	%*Create an empty dataset that we are going to append the records too. We need to manually set the length
	and formats so that we avoid errors.;
	data _inpt_adm;
		length enrolid 8 admdate 4;
		format enrolid best12. admdate MMDDYY10.;
		stop;
	run;
			
	%*Now grab the medications for asthma;
	%get_inpt_adm(input_pregs=_pregnancies_pnc, base=_inpt_adm);

	%*Add sort statement to make certain that there are no duplicate rows (though there should not be);
	proc sort data=_inpt_adm nodup; by enrolid; run; /*CDL: 1.28.2025 - added to be sure that there were no issues*/

	%*Merge this information onto the pregnancies dataset that created, which now has the distinct outpatient and inpatent PNC encounter dates;
	proc sql;
		create table &OUTPUT_DATA as
		select distinct a.*, count(distinct c.admdate) as num_InptAdm
		from _pregnancies_pnc as a
		left join _inpt_adm as c
		on a.enrolid=c.enrolid and a.dt_lookback <= c.admdate <= a.dt_index
		group by a.enrolid, a.idxpren
		;
		quit;

	%*Delete unnecessary datasets;
	proc datasets gennum=all noprint;
		delete _pregnancies_pnc _inpt: ;
	run;

%mend;








/*
MACRO: get_all_covariates
PURPOSE: To define all of the necessary covariates for analyses

INPUTS: 
LMPINDEX -- Gestational age assumed at the indexing date for pregnancies with UNK outcomes and no GA information
LOOKBACKDT -- The date from which we at looking back to define covarites: dt_index (primary) and dt_lmp (sensitivity analysis)
LOOKBACKDAYS -- The number of days that we are looking back from the lookbackdt to define covariates
*/

%macro get_allcovariates(LMPINDEX, LOOKBACKDT, LOOKBACKDAYS);

	/*For testing:
	%let lmpindex = 63;
	%let lookbackdt = dt_index;
	%let lookbackdays = 270;
	*/

	%*Identify all the years that we will need to pull sample data from;
	%let years = 2016 2017 2018 2019 2020 2021 2022;

	%*Pull in the pregnancies file that we will be working with;
	data _pregnancies;
	set temp.preg_incl_excl_&lmpindex._&lookbackdt._&lookbackdays;

		%*Make sure remove pregnancies with their indexing prenatal enocunter after their pregnancy outcome date;
		/*CDL: 1.23.2025 -- Included this statement here becuase deletion originally not included in the first data cleaning step.
		Should not affect any after a run that includes that code.*/
		if dt_indexprenatal > dt_gapreg then delete;
	run;

	%**********
		STEP 1: Identify information relevant to diagnosis codes. For most diagnoses, we just want the number of 
		diagnosis that occurred <= 270 day prior to or <= 30 days after the index date. However, antiphospholipid
		syndrome and diabetes, we implement more complex algorithms
	**********;

	%*Combine all the diagnosis codes with the relevant covariate information in a column labeled diagnosis;
	%combine_codes(cov_diagnoses, "DX10");


	%*First, identify all the relevant diagnosis code based information for diagnoses with simple algorithms that
	rely on counts only. Implement same macro as in the inclusion and exclusion criteria

	Covariates that STILL NEED TO BE HANDLED:
	- diabetes
	- antiphosphilipid syndrome;
	%identify_dx_covariates(REF=cov_diagnoses, PREG_INPT=_pregnancies, OUTPUT=_pregnancies_dx);


	%*******
		DIABETES: We need additional variables to assess this. Work off of the _diagnoses2 dataset, as
		it has all of the diagnosis codes that matched with the list above;

	%identify_diabetes_dx(INPUT_DATA=_pregnancies_dx, OUTPUT_DATA=_pregnancies_diab, DIAGNOSES_DATA=_diagnoses2);

	%*******
		HYPERTENION: We need additional variables on hypertension diagnoses that are indexed on the LMP at which they occur;

	%identify_htn_dx(_pregnancies_diab, _pregnancies_htn, DIAGNOSES_DATA=_diagnoses2);


	%***Combine all the procedure codes with the relevant covariate information in a column labeled diagnosis;

	%combine_codes(cov_procedures, "CPT" "HCPCS" "PR10" );

	%***Identify the relevant covariates based upon procedure codes;

	%identify_pr_covariates(_pregnancies_htn, _pregnancies_pr, cov_procedures);

	%***Combine all the medication order NDC codes;
	proc sql;
		create table medications as
		select distinct atc_label, ndc9, "ADHD" as medication from rxcov.adhd_meds_rx
		union corr
		/*CDL: ADDED 4.16.2025 -- Wanted to be able to identify stimulants only*/
		select distinct atc_label, ndc9, "Simulant" as medication from rxcov.adhd_stimulants_rx
		union corr
		select distinct atc_label, ndc9, "Anticonvul" as medication from rxcov.anticonvulsants_rx
		union corr
		select distinct atc_label, ndc9, "Antidep" as medication from rxcov.antidepressants_rx
		union corr
		select distinct atc_label, ndc9, "Antipsy" as medication from rxcov.antipsychotics_rx
		union corr
		select distinct atc_label, ndc9, "Benzo" as medication from rxcov.benzodiazepines_rx
		union corr
		select distinct atc_label, ndc9, "HyperThy" as medication from rxcov.hyperthyroid_rx
		union corr
		select distinct atc_label, ndc9, "Nausea" as medication from rxcov.nausea_vomiting_rx
		union corr
		select distinct atc_label, ndc9, "HypoThy" as medication from rxcov.hypothyroid_rx
		union corr
		select distinct atc_label, ndc9, "Metfor" as medication from rxcov.metformin_rx
		union corr
		select distinct atc_label, ndc9, "Migraine" as medication from rxcov.migraine_rx
		union corr
		select distinct atc_label, ndc9, "MoodStab" as medication from rxcov.mood_stabilizers_rx
		union corr
		select distinct atc_label, ndc9, "Oud" as medication from rxcov.oud_rx
		union corr
		select distinct atc_label, ndc9, "PTSD" as medication from rxcov.ptsd_meds_rx
		union corr
		select distinct atc_label, ndc9, "Smk" as medication from rxcov.smoking_rx
		union corr
		select distinct atc_label, ndc9, "Statin" as medication from rxcov.statins_rx
		union corr
		select distinct atc_label, ndc9, "t1t2dm" as medication from rxcov.t1t2dm_antidiabetics_rx
		union corr
		select distinct atc_label, ndc9, "t2dm" as medication from rxcov.t2dm_antidiabetics_rx
		union corr
		select distinct atc_label, ndc9, "terat" as medication from rxcov.teratogenic_meds_rx
		union corr
		select distinct atc_label, ndc9, "anxiety" as medication from rxcov.anxiety_rx
		union corr
		select distinct atc_label, ndc9, "glp1wgt" as medication from rxcov.wgtloss_glp1_rx
		union corr
		select distinct atc_label, ndc9, "otherwgt" as medication from rxcov.wgtloss_other_rx
		union corr
		/*CDL: ADDED 11.16.2025 all rows below*/
		select distinct atc_label, ndc9, "lda" as medication from rxcov.lda_rx
		union corr
		select distinct atc_label, ndc9, "aspirin" as medication from rxcov.nsaid_aspirin_rx
		union corr
		select distinct atc_label, ndc9, "nsaid" as medication from rxcov.nsaid_rx
		;
		quit;

	%*Now use the macro to identify the RX covariates;

	%identify_rx_covariates(_pregnancies_pr, _pregnancies_cov_rx, medications);

	%*Identify healthcare utilization variables. These include number of prior outpatient prenatal encounters
	and number of prior inpatient admissions;

	%get_healthcare_use(_pregnancies_cov_rx, _pregnancies_health_use);


	%*Output the final dataset;
	data temp.preg_del_covar_&lmpindex._&lookbackdt._&lookbackdays;
	set _pregnancies_health_use;
	run;

	%*Delete the unnecessary dataset;
	proc datasets gennum=all noprint;
		delete _:;
	run;

%mend;



















/********************************************************************************************************************************************

										04 - IDENTIFY POST-INDEX ANTIHYPERTENSIVES, PREECLAMPSIA, and PTB

Derive post-index variables:
- First fill for an antihypertensive other than the exposure
- First fill for the other study medication
- First date where we see a gap in usage: days supply + 0, 7, 14, and 30 d gap
- The number of fills for non-exposure antihypertensive
- The number of distinct ATC labels.
- First preeclampsia diagnosis after the index date, on or before the outcome date.
- The number of indicated preterm birth procedure codes 7 days before or after the outcome date
- The number of spontaneous preterm birth diagnosis codes 7 days before or after the outcome date

********************************************************************************************************************************************/


%macro get_dx_outcomes(input_pregs, input_ref, base);

	%let numYr = %sysfunc(countw(&years));

	%do d=1 %to &numYr;

		%let loop&d = %scan(&years, &d);

		%*Create dataset with all relevant diagnosis codes for identified pregnancy for the year defined by &loop&d;
		proc sql;
			create table _diagnoses_subset as
			select distinct a.enrolid format=best12. length=8, a.svcdate format=MMDDYY10. length=8, 
				c.idxpren format=best12. length=8, c.dt_index format=MMDDYY10. length=8, 
				a.dxLoc format=$9. length=9, a.dxNum format=best12. length=8, b.code format=$7. length=7, 
				b.diagnosis format=$8. length=8
			from %if &&loop&d = 2016 %then %do; der.alldx102016 %end; %else %do; der.alldx&&loop&d. %end; as a
			inner join &input_pregs as c
			on a.enrolid=c.enrolid and c.dt_index < a.svcdate <= c.dt_gapreg+14
			inner join &input_ref as b
			on a.dx&&loop&d. = b.code
			;
			quit;

		%*Append the dataset onto the base dataset;
		proc append base=&base data=_diagnoses_subset; run;

	%end;

%mend;





/*
MACRO: identify_htn_outcomes
PURPOSE: The purpose of this macro is to identify the occurrence of preeclampsia after the index
date. We are interested in the first diagnosis.

INPUTs:
- INPUT_DATA = input pregnancy datset
- OUTPUT_DATA = output pregnancy datset
- DIAGNOSES_DATA = dataset with diagnosis codes from the claims data
*/

%macro identify_htn_outcomes(INPUT_DATA, OUTPUT_DATA);

	%*Create the reference codelist for preeclampsia;
	proc sql;
		create table preeclampsia as
		select distinct code, "Preec" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.preeclampsia_dx where codetype = "DX10"
		;
		quit;

	%*Now get the relevant diagnosis codes for preeclampsia;
	data _diagnoses_outcomes;
		length enrolid 8 svcdate 8 idxpren 8 dt_index 8 dxLoc $9 dxNum 8 code $7 diagnosis $8;
		format enrolid best12. svcdate MMDDYY10. idxpren best12. dt_index MMDDYY10. dxLoc $9. dxnum best12. code $7. diagnosis $8.;
		stop;
	run;

	%get_dx_outcomes(input_pregs=&INPUT_DATA, input_ref = preeclampsia, base = _diagnoses_outcomes);

	%*Do some data cleaning on the codes;
	data _diagnoses2;
	set _diagnoses_outcomes;
		%*Determine if the diagnosis code was on an outpatient record (outpatient service claim) or inpatient record (inpatient service or
			inpatient admission claim);
		if dxLoc = "OutptServ" then location = "Outpt";
			else if dxLoc in ("InptServ" "InptAdm") then location = "Inpt";
			else location = "";
		if location = "" then delete; /*Should not be a problem*/
	run;

	%*Now, we want the first instance after the index date, on or up to 2 weeks after the outcome date;
	proc sql;
		create table _diagnoses_outcomes2 as
		select distinct enrolid, idxpren, location, min(svcdate) as dt_preeclampsia format=MMDDYY10.
		from _diagnoses2
		group by enrolid, idxpren, location
		;
		quit;

	%*Transpose the dataset;
	proc transpose data=_diagnoses_outcomes2 out=_diagnoses2_tranpose(drop = _NAME_) prefix=dt_preec_outc;
		by enrolid idxpren;
		var dt_preeclampsia;
		id location;
	run;

	%*Now add those dates onto the output pregnancies dataset;
	proc sql;
		create table &OUTPUT_DATA as
		select a.*, b.dt_preec_outcOutpt, b.dt_preec_outcInpt
		from &INPUT_DATA as a
		left join _diagnoses2_tranpose as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		;
		quit;

	%*Delete the unnecessary datasets;
	proc datasets gennum=all noprint;
		delete _diagnoses: ;
	run; quit; run;

%mend;



/*
MACRO: get_proc_ptb_outcomes
PURPOSE: To identify all procedure codes used to identify indicated preterm birth from the claims data
*/


%macro get_proc_ptb_outcomes(input_pregs, input_ref, base);

	%let numYr = %sysfunc(countw(&years));

	%do d=1 %to &numYr;

		%let loop&d = %scan(&years, &d);

		%*Create dataset with all relevant diagnosis codes for identified pregnancy for the year defined by &loop&d;
		proc sql;
			create table _proc_subset as
			select distinct a.enrolid format=best12. length=8, a.svcdate format=MMDDYY10. length=8, 
				c.idxpren format=best12. length=8, c.dt_index format=MMDDYY10. length=8, 
				a.procLoc format=$9. length=9, b.code format=$7. length=7, 
				b.diagnosis format=$8. length=8
			from der.allproc&&loop&d. as a
			inner join &input_pregs as c
			on a.enrolid=c.enrolid and c.dt_gapreg - 7 <= a.svcdate <= c.dt_gapreg + 7
			inner join &input_ref as b
			on a.proc&&loop&d. = b.code
			;
			quit;

		%*Append the dataset onto the base dataset;
		proc append base=&base data=_proc_subset; run;

	%end;

%mend;




/*
MACRO: identify_ptb_outcomes
PURPOSE: The purpose of this macro is to identify diagnosis codes related to indicated
versus spontaneous preterm birth within 7 days of the index date.

INPUTs:
- INPUT_DATA = input pregnancy datset
- OUTPUT_DATA = output pregnancy datset
- DIAGNOSES_DATA = dataset with diagnosis codes from the claims data
*/

%macro identify_ptb_outcomes(INPUT_DATA, OUTPUT_DATA);


	%**********************
		Spontaneous PTB;


	%*Create the reference codelist for spontaneous preterm birth;
	proc sql;
		create table spontaneous_ptb as
		select distinct code, "SponPTB" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.ptb_spontaneous_dx 
		;
		quit;

	%*Now get the relevant diagnosis codes for spontaneous preterm birth;
	data _diagnoses_outcomes;
		length enrolid 8 svcdate 8 idxpren 8 dt_index 8 dxLoc $9 dxNum 8 code $7 diagnosis $8;
		format enrolid best12. svcdate MMDDYY10. idxpren best12. dt_index MMDDYY10. dxLoc $9. dxnum best12. code $7. diagnosis $8.;
		stop;
	run;

	%get_dx_outcomes(input_pregs=&INPUT_DATA, input_ref = spontaneous_ptb, base = _diagnoses_outcomes);

	%*Do some data cleaning on the codes;
	data _diagnoses2;
	set _diagnoses_outcomes;
		%*Determine if the diagnosis code was on an outpatient record (outpatient service claim) or inpatient record (inpatient service or
			inpatient admission claim);
		if dxLoc = "OutptServ" then location = "Outpt";
			else if dxLoc in ("InptServ" "InptAdm") then location = "Inpt";
			else location = "";
		if location = "" then delete; /*Should not be a problem*/
	run;

	%*Add the outcome date information back onto the datset;
	proc sql;
		create table _diagnoses2_outc as
		select a.*, b.dt_gapreg
		from _diagnoses2 as a
		left join &INPUT_DATA as b
		on a.idxpren = b.idxpren
		;
		quit;

	%*Now, only retain those rows where the service date occurs within 7 days of the outcome date;
	data _diagnoses3;
	set _diagnoses2_outc;

		%*Create a variable for if the code was pre (including outcome date) or post outcome;
		if svcdate <= dt_gapreg then timing = "pre";
			else timing = "post";

		if dt_gapreg - 7 <= svcdate and svcdate <= dt_gapreg + 7 then output;
	run;

	%*Now count all instances up to and including the outcome date and after the outcome date.;
	proc sql;
		create table _spon_ptb_counts as
		select distinct enrolid, idxpren, location, timing, count(distinct svcdate) as numdx
		from _diagnoses3
		group by enrolid, idxpren, location, timing
		;
		quit;


	%*Transpose the dataset so that each diagnosis, location, and timing  (pre v post index) are their own column.
	These columns contain counts of the diagnoses that satisfy the requirement.
	NOTE: If no one has a code that meets a criteria, then that column will not be created.;
	proc transpose data=_spon_ptb_counts out=_spon_ptb_counts_trans (drop = _NAME_) prefix=sponptb_;
		by enrolid idxpren;
		id location timing;
		var numdx;
	run;




	%**********************
		Indicated PTB;


	%*Create the reference codelist for spontaneous preterm birth;
	proc sql;
		create table indicated_ptb as
		select distinct code, "IndPTB" as diagnosis, 0 as exclusion_criteria, "" as outcome format=$57. length=57, "" as algorithm_step format=$9. length=9
		from covref.ptb_indicated_dx 
		;
		quit;

	%*Now get the relevant diagnosis codes for spontaneous preterm birth;
	data _proc_outcomes;
		length enrolid 8 svcdate 8 idxpren 8 dt_index 8 procLoc $9 code $7 diagnosis $8;
		format enrolid best12. svcdate MMDDYY10. idxpren best12. dt_index MMDDYY10. procLoc $9. code $7. diagnosis $8.;
		stop;
	run;

	%get_proc_ptb_outcomes(input_pregs=&INPUT_DATA, input_ref = indicated_ptb, base = _proc_outcomes);

	%*Do some data cleaning on the codes;
	data _procedures2;
	set _proc_outcomes;
		%*Determine if the procedure code was on an outpatient record (outpatient service claim) or inpatient record (inpatient service or
			inpatient admission claim);
		if procLoc = "OutptServ" then location = "Outpt";
			else if procLoc in ("InptServ" "InptAdm") then location = "Inpt";
			else location = "";
		if location = "" then delete; /*Should not be a problem*/
	run;

	%*Add the outcome date information back onto the datset;
	proc sql;
		create table _procedures2_outc as
		select a.*, b.dt_gapreg
		from _procedures2 as a
		left join &INPUT_DATA as b
		on a.idxpren = b.idxpren
		;
		quit;

	%*Now, only retain those rows where the service date occurs within 7 days of the outcome date;
	data _procedures3;
	set _procedures2_outc;

		%*Create a variable for if the code was pre (including outcome date) or post outcome;
		if svcdate <= dt_gapreg then timing = "pre";
			else timing = "post";

		if dt_gapreg - 7 <= svcdate <= dt_gapreg + 7 then output;
	run;

	%*Now count all instances up to and including the outcome date and after the outcome date.;
	proc sql;
		create table _ind_ptb_counts as
		select distinct enrolid, idxpren, location, timing, count(distinct svcdate) as numpr
		from _procedures3
		group by enrolid, idxpren, location, timing
		;
		quit;


	%*Transpose the dataset so that each diagnosis, location, and timing  (pre v post index) are their own column.
	These columns contain counts of the diagnoses that satisfy the requirement.
	NOTE: If no one has a code that meets a criteria, then that column will not be created.;
	proc transpose data=_ind_ptb_counts out=_ind_ptb_counts_trans (drop = _NAME_) prefix=indptb_;
		by enrolid idxpren;
		id location timing;
		var numpr;
	run;


	%*Now add those dates onto the output pregnancies dataset -- Return just the minimum information;
	proc sql;
		create table &OUTPUT_DATA as
		select a.enrolid, a.idxpren, b.*, c.* /*This will cause warnings becuase some duplicate variables across but are matching on those*/
		from &INPUT_DATA as a
		left join _spon_ptb_counts_trans as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		left join _ind_ptb_counts_trans as c
		on a.enrolid=c.enrolid and a.idxpren=c.idxpren
		;
		quit;

	%*Delete the unnecessary datasets;
	proc datasets gennum=all noprint;
		delete _diagnoses: _procedures: _spon: _ind:;
	run; quit; run;

%mend;








/*
MACRO: get_postdx_antihypertensives
PURPOSE: Go through all outpatient fills to identify those fills for antihypertensives
that are either nifedipine, labetalol, or neither. This loops through all the separate
year files from the raw data. In particular, we want to identify all those fills that occurred
between the index date and the minimum .
*/
%macro get_postidx_antihypertensives(INPUT_DATA, BASE);

	%let numYr = %sysfunc(countw(&years));

	%do d=1 %to &numYr;

		%let loop&d = %scan(&years, &d);

		%*This is a subset of all the antihypertensive fills where a dispensing quantity was greater than 0. This is
		created for each year and then appended (i.e., proc append) onto the _prior_antihypertensives;
		proc sql;
			create table _antihtn_subset as
			select a.enrolid as enrolid format=best12. length=8, a.svcdate as svcdate format=MMDDYY10. length=8,
				a.ndcnum as ndc11 format=$11. length=11, a.metqty as metqty format=best12. length=4,
				a.daysupp as daysupp format=best12. length=3,
				b.idxpren as idxpren format=BEST12. length=8, b.dt_index as dt_index format=MMDDYY10. length=8, 
				b.exposure as exposure format=$56. length=56, c.atc_label as atc_label format=$56. length=56,
				c.ndc9 as ndc9 format=$9. length=9
			from (select distinct * from raw.outptdrug&&loop&d where metqty > 0) as a
			inner join &INPUT_DATA as b /*Only grab those fills for pregnancies in the current dataset*/
			on a.enrolid = b.enrolid and b.dt_index <= a.svcdate <= b.dt_gapreg /*Fills must occur between the lookback date and the index date*/
			inner join (select distinct atc_label, ndc9 from rxcov.antihypertensives_rx) as c /*Link based upon NDC-9 code, as there may be missed packagings with NDC-11*/
			on substr(a.ndcnum, 1, 9) = c.ndc9
			;
			quit;

		proc append base=&BASE data=_antihtn_subset; run;

		/*Output the number of pregnancies in that year*/
		proc sql noprint;
			select count(*) as num_pregnancies_loop_year
			into :num_pregnancies_loop_year
			from _antihtn_subset;
			quit;

		%put Number of pregnancies with qualifying antihypertensives in &&loop&d: &num_pregnancies_loop_year;

	%end;

%mend;



/*
MACRO: identify_postid_antihypertensives
PURPOSE: To identify key variables that describe post-index antihypertensive fills

INPUTS:
LMPINDEX -- Assumed gestational age at index date for UNK pregnancies with no GA information
LOOKBACKDT -- Date from which we apply lookback period: dt_index (primary) and dt_lmp (sensitivity)
LOOKBACKDAYS -- Lookback period from LOOKBACKDT in days
*/

%macro identfy_postid_antihypertensives(LMPINDEX, LOOKBACKDT, LOOKBACKDAYS);

	/*For testing:
	%let lmpindex = 63;
	%let lookbackdt = dt_index;
	%let lookbackdays = 270;
	*/

	%*Identify all the years that we will need to pull sample data from;
	%let years = 2016 2017 2018 2019 2020 2021 2022;

	%*Upload the pregnancies dataset;
	data _pregnancies;
	set temp.preg_covar_&lmpindex._&lookbackdt._&lookbackdays;
	run;

	%*Get the post-index preeclampsia outcomes -- Added 3.27.2025;
	%identify_htn_outcomes(INPUT_DATA=_pregnancies, OUTPUT_DATA=_postidx_preec);

	%*Get the preterm birth information round the outcome date -- Added 4.16.2025;
	%identify_ptb_outcomes(INPUT_DATA = _pregnancies, OUTPUT_DATA=_postidx_ptb)

	%*Get all antihypertensive fills that occurred between the index date and dt_gapreg;

	%*Create an empty dataset that we are going to append the records too. We need to manually set the length
	and formats so that we avoid errors.;
	data _postidx_antihypertensives;
	    length enrolid 8 svcdate 8 ndc11 $11 metqty 4 daysupp 3 idxpren 8 dt_index 8 exposure $56 atc_label $56 ndc9 $9;
	    format enrolid best12. svcdate MMDDYY10. ndc11 $11. metqty best12. daysupp best12. idxpren best12. dt_index MMDDYY10. exposure $56. atc_label $56. ndc9 $9.;
	    stop;
	run;

	%get_postidx_antihypertensives(INPUT_DATA=_pregnancies, BASE=_postidx_antihypertensives);

	%*Explicit statement to remove all duplicates, though there should not be any. CDL: ADDED 1.28.2025;
	proc sort data=_postidx_antihypertensives nodup; by enrolid idxpren; run;
	
	%*We only want to save this dataset for the index date.;

	%if &lookbackdt = dt_index %then %do;

		%*Save this dataset in case needed;
		data temp.postidx_antihypertensives_&lmpindex; set _postidx_antihypertensives; run;

	%end;

	%*****Identify key variables;

	%*Identify the date of the first fill for an antihypertensive other than the exposure;
	proc sql;
		create table first_fill_any as
		/*CDL: 1.28.2025 - added distinct*/
		select distinct enrolid, idxpren, min(svcdate) as dt_first_nonexp_antihtn format=MMDDYY10.,
			count(*) as num_nonexp_antihtn, count(distinct atc_label) as num_distinct_nonexp_antihtn
		from _postidx_antihypertensives (where = (exposure ne atc_label))
		group by enrolid, idxpren
		;
		quit;

	%*First fill for the other study medication;
	data _postidx_antihypertensives2;
	set _postidx_antihypertensives;
		if exposure = 'NIFEDIPINE' and ATC_LABEL = "LABETALOL" then other_exposure = 1;
			else if exposure = "LABETALOL" and ATC_LABEL = "NIFEDIPINE" then other_exposure = 1;
			else other_exposure = 0;
	run;

	proc sql;
		create table _other_exposure as
		/*CDL: 1.28.2025 - added distinct*/
		select distinct enrolid, idxpren, min(svcdate) as dt_first_otherexp_antihtn format=MMDDYY10.,
			count(*) as num_otherexp_antihtn
		from _postidx_antihypertensives2 (where = (other_exposure = 1))
		group by enrolid, idxpren
		;
		quit;


	%*Identify the first date where we see a gap in usage of &medgap;

	%*Subset to only those fills for the exposure;
	data _exposure;
	set _postidx_antihypertensives;
		where exposure = atc_label;
	run;

	proc sort data=_exposure;
		by enrolid idxpren svcdate;
	run;
	data _exposure2;
	set _exposure;
		format next_fill MMDDYY10. last_dt MMDDYY10.;
		by enrolid idxpren svcdate;
		retain count0 count7 count30 count45;

		last_id = lag1(idxpren);
		next_fill = svcdate + daysupp;
		last_dt = lag1(next_fill);

		if first.idxpren then count0 = 1;
			else if idxpren = last_id and svcdate - last_dt > 0 then count0 = count0+1;
			else count0 = count0;	

		if first.idxpren then count7 = 1;
			else if idxpren = last_id and svcdate - last_dt > 7 then count7 = count7+1;
			else count7 = count7;

		if first.idxpren then count30 = 1;
			else if idxpren = last_id and svcdate - last_dt > 30 then count30 = count30+1;
			else count30 = count30;

		if first.idxpren then count45 = 1;
			else if idxpren = last_id and svcdate - last_dt > 45 then count45 = count45+1;
			else count45 = count45;
	run;

	%*Now the date of the first gap will be the next_fill date for the last row where count = 1;
	data _exposure_0;
	set _exposure2 (where = (count0 = 1));
		by enrolid idxpren svcdate;
		if last.idxpren then output;
	run;

	data _exposure_7;
	set _exposure2 (where = (count7 = 1));
		by enrolid idxpren svcdate;
		if last.idxpren then output;
	run;

	data _exposure_30;
	set _exposure2 (where = (count30 = 1));
		by enrolid idxpren svcdate;
		if last.idxpren then output;
	run;

	data _exposure_45;
	set _exposure2 (where = (count45 = 1));
		by enrolid idxpren svcdate;
		if last.idxpren then output;
	run;


	%*Now put all of these variables onto the pregnancy dataset;
	proc sql;
		create table out.preg_del_cohort_&lmpindex._&lookbackdt._&lookbackdays as
		/*CDL: 1.28.2025 - added distinct*/
		select distinct a.*, b.dt_first_nonexp_antihtn, b.num_nonexp_antihtn, b.num_distinct_nonexp_antihtn,
			c.dt_first_otherexp_antihtn, c.num_otherexp_antihtn,
			e.svcdate as dt_lastfill_gap0, e.next_fill as dt_lastfill_daysupp_gap0 format=MMDDYY10.,
			f.svcdate as dt_lastfill_gap7, f.next_fill as dt_lastfill_daysupp_gap7 format=MMDDYY10.,
			g.svcdate as dt_lastfill_gap30, g.next_fill as dt_lastfill_daysupp_gap30 format=MMDDYY10.,
			h.svcdate as dt_lastfill_gap45, h.next_fill as dt_lastfill_daysupp_gap45 format=MMDDYY10.,
			i.dt_preec_outcOutpt, i.dt_preec_outcInpt, j.* /*This will cause some warnings on enrolid and idxpren.*/
		from _pregnancies as a
		left join first_fill_any as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		left join _other_exposure as c
		on a.enrolid=c.enrolid and a.idxpren=c.idxpren
		left join _exposure_0 as e
		on a.enrolid=e.enrolid and a.idxpren=e.idxpren
		left join _exposure_7 as f
		on a.enrolid=f.enrolid and a.idxpren=f.idxpren
		left join _exposure_30 as g
		on a.enrolid=g.enrolid and a.idxpren=g.idxpren
		left join _exposure_45 as h
		on a.enrolid=h.enrolid and a.idxpren=h.idxpren
		left join _postidx_preec as i
		on a.enrolid=i.enrolid and a.idxpren=i.idxpren
		left join _postidx_ptb as j
		on a.enrolid=j.enrolid and a.idxpren=j.idxpren
		;
		quit;

%mend;





/************************************************************************************************************

											05 - PULL THE COHORT DATA

************************************************************************************************************/

%*******gestational age at index is 63 days and restrict to deliveries


%********gestational age at index is 63 days;
%derive_initial_cohort(lmpindex = 63); *maxfill = 36*7+6 or 36w6d;
proc datasets gennum=all noprint; delete _:; run;

%*Index date analysis;
%identify_new_users(63, lookbackdt = dt_lmp, lookbackdays = 180);
%define_incl_excl(lmpindex=63, gap=31, lookbackdt = dt_lmp, lookbackdays=180);
%get_allcovariates(lmpindex=63, lookbackdt = dt_lmp, lookbackdays=180);
%identfy_postid_antihypertensives(lmpindex=63, lookbackdt = dt_lmp, lookbackdays=180);

