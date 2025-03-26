/********************************************************************************************************************************************
PROGRAM: 01b_clean_preg_cohort.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to clean the pregnancy cohort that we derived from the MarketScan claims data..
	
Goal: 
Output data: 

Date: 12.19.2024
********************************************************************************************************************************************/











/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 02 - MACROS FOR TOO LONG PREGNANCIES
	- 03 - MACROS FOR INDEX PRENATAL TOO EARLY
	- 04 - MACROS FOR OVERLAP
	- 05 - CLEAN UP PREGNANCIES

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
/*change "saveLog=" to "Y" when program is closer to complete*/
%setup(sample=full, programname=01b_clean_preg_cohort, savelog=N);


******************************************************************************************************************************************;
/*Create local mirros of the libraries from the set up macro - Run locally*/
/*libname lraw slibref=raw server = server;*/
/*libname lder slibref=der server = server;*/
/*libname lout slibref=out server=server; */
/*libname lana slibref=ana server=server;*/
/*libname lpreg slibref=preg server=server;*/
/*libname lwork slibref=work server=server;*/
/*libname lcovref slibref=covref server=server;*/
/*libname ltemp slibref=temp server=server;*/

*Run this format file locally (Ctrl+A on the file with local run) if you want to view datasets;
/*%inc "/local/projects/marketscan_preg/Latour_23_2322/programs/FormatStatements_CDWH.sas";*/










/********************************************************************************************************************************************

													02 - MACROS FOR TOO LONG PREGNANCIES

Some pregnancies were identified as having too long of a gestational length. These have to be revised so that they are not too long.
This only occurred to pregnancies with UNK outcomes.

********************************************************************************************************************************************/

/*
MACRO: revise_too_long
PURPOSE: The goal of this macro is to revise date-based information for those pregnancies
that are identified as being too long.

The process should look as follow:
1. Subset to those pregnancies with UNK outcomes with a gestational length gt 307d
2. For those without any GA codes, reassign their LMP by assuming that they were 4w0d
	gestation at their indexing prenatal encounter.
3. Reapply the gestational age algorithm but assume the minimum GA for the
	gestational age in days. Ignore prenatal care codes. For any min GAs of 14, reset
	those to 28
	- Those that no longer have GA information  (i.e., only GA codes were prenatal care),
	reassign their LMP by assuming they were 4w0d gestation at their indexing prenatal enc.
4. For those pregnancies that continue to be too long, assign their LMP as the latest LMP 
	calculated by the algorithm.
5. For those pregnancies that continue to be too long:
	- If they do not have a Z3A cde, assign their LMP by assumign the GA at their indexing
	prenatal encounter was 4w0d
	- If they have a Z3A code, revise the index and LTFU dates. Grab all prenatal encounters that 
	occur between the maximum LMP (from the algorithm using the minimum gestational age) and 307 
	days gestation. Reassign the index date and LTFU accordingly.

INPUTS:
- INPUT_DATA = INput pregnancy dataset

%let INPUT_DATA=_pregnancies2;
%let days = 307;
%let lmpindex = 63;
*/



%macro revise_too_long(INPUT_DATA, output_data, LMPINDEX, DAYS=307, round=0);

	%*Save a proc freq output of the original dataset as an analytic dataset;
	proc freq data=&input_data noprint;
		table preg_outcome_clean / missing out=ana.preg_outc_pre_cleaning_&lmpindex._&round;
	run;

	%*First, subset to the pregnancies of interest;
	data _pregnancies_long;
	set &INPUT_DATA;
		ga_length = dt_gapreg - dt_lmp;
		if preg_outcome_clean = "UNK" and ga_length > &DAYS then output;
	run;

	%*Save a proc freq output as an analytic dataset;
	proc freq data=_pregnancies_long noprint;
		table preg_outcome_clean / missing out=ana.preg_outc_toolong_&lmpindex._&round;
	run;
	/*Output the counts of pregnancies with too long of gestational lengths*/
	proc sql noprint;
		select count(distinct idxpren) into :num_preg_toolong from _pregnancies_long;
		quit;
	%put Number of pregnancies with too long gestational length: &num_preg_toolong;

	%*Identify those pregnancies with no GA information and revise their LMP by
	assuming their gestational age was 4wk at their indexing prenatal encounter;
	data _pregnancies_long_noga;
	set _pregnancies_long;
		dt_lmp = dt_indexprenatal - 28;
		if LMP_AlgBased = 0 then output;
	run;

	%*For the rest of the pregnancies, revise their LMP estimates using the minimum
	gestational age indicated by a code. Prenatal care codes do not have a minimum gestational age, so they are not used.;
	%*This macro is stored in the macros folder within programs;
	options nomprint;
	%getGA(INPUT_DATA=_pregnancies_long (where = (lmp_algbased = 1)), output_data= ga_rev, GA=min_gest_age);
	options mprint;

	%*Reassign the LMPs based on the revisions.;
	data _pregnancies_long_rev2;
	set ga_rev;
		%*Assign the LMP so that it is not missing;
		if lmp_algbased = 1 then dt_lmp = dt_lmp_alg;
			else if lmp_algbased = 0 and preg_outcome_clean ne "UNK" then dt_lmp = dt_lmp_table; %*Should not be in this dataset;
			else if lmp_algbased = 0 then dt_lmp = dt_indexprenatal - &lmpindex;

		ga_length = dt_gapreg - dt_lmp; %*Calculate new GA length variable;
	run;

	%*Output a dataset of pregnancies with revised LMPs that no longer have too long of gestational lengths;
	data _pregnancies_long_rev2_out;
	set _pregnancies_long_rev2;
		where ga_length <= &DAYS;
	run;

	%*Refine the LMPs for those pregnancies with GA length that continues to be too long.
	Specifically, assign the latest (max) LMP provided by the updated algorithm;
	data _pregnancies_long2;
	set _pregnancies_long_rev2 (where = (ga_length > &DAYS));
		if dt_LMP_Max ne . then dt_lmp = dt_LMP_Max;
			else dt_lmp = dt_lmp;

		%*Revised galength;
		ga_length2 = dt_gapreg - dt_lmp;
	run;

	%*Output those pregnancies where the max LMP returns a reasonable gestational length
	once we use the maximum LMP output by the algorithm;
	data _pregnancies_long2_max;
	set _pregnancies_long2;
		where ga_length2 <= &DAYS;
	run;

	%*For those pregnancies with GA information that are still too long, identify those with Z3A codes;

	%*Combine back with the temp.gestagepren datset to get Z3A codes and determine if the person has 
	at least 1 Z3A code;
	proc sql;
		create table _pregnancies_long3 as
		select distinct a.*, sum(b.code_hierarchy ne "") as numZ3A
		from _pregnancies_long2 (where = (ga_length2 > &DAYS)) as a
		left join temp.gestagepren (where = (code_hierarchy = "Specific gestational age")) as b
		on a.patient_deid = b.patient_deid and a.dt_indexprenatal <= b.enc_date <= a.dt_gapreg + 7
		group by idxpren
		;
		quit;

	%*For those pregnancies without Z3A codes, assign their LMP as 4 weeks prior to the indexing date;
	data _pregnancies_long4;
	set _pregnancies_long3 (where = (numZ3A = 0));
		dt_lmp = dt_indexprenatal - 28;
	run;

	%*For those pregnancies with an LMP based on Z3A codes, revise the indexprenatal and outcome date based on
	all prenatal encounter dates between their LMP and &DAYS gestation;
	proc sql;
		create table _pregnancies_long_prenatal as
		select distinct a.*, min(b.enc_date) as dt_indexprenatal_rev, max(b.enc_date) as dt_gapreg_rev
		from _pregnancies_long3 (where = (numZ3A >= 1)) as a
		left join temp.codeprenatal_meg1_dts as b
		on a.patient_deid=b.patient_deid and a.dt_lmp <= b.enc_date <=  a.dt_lmp + &days
		group by a.patient_deid, a.idxpren
		;
		quit;

	%*Now revise the variables;
	data _pregnancies_long_prenatal3;
	set _pregnancies_long_prenatal;
		dt_indexprenatal = dt_indexprenatal_rev;
		dt_gapreg = dt_gapreg_rev;
		drop dt_indexprenatal_rev dt_gapreg_rev;
	run;

	%*Final dataset;
	proc sql;
		create table pregnancies_long_rev as
		select distinct idxpren, dt_lmp, dt_indexprenatal, dt_gapreg, 
				1 as LMP_long_revision, "No GA codes" as LMP_long_revision_reason
		from _pregnancies_long_noga
		union corr
		select distinct idxpren, dt_lmp, dt_indexprenatal, dt_gapreg, 
				1 as LMP_long_revision, "Reasonable with min_gest_age in algorithm" as LMP_long_revision_reason
		from _pregnancies_long_rev2_out
		union corr
		select distinct idxpren, dt_lmp, dt_indexprenatal, dt_gapreg, 
				1 as LMP_long_revision, "Reasonable if use max_dt_lmp in algorithm" as LMP_long_revision_reason
		from _pregnancies_long2_max
		union corr
		select distinct idxpren, dt_lmp, dt_indexprenatal, dt_gapreg, 
				1 as LMP_long_revision, "No Z3A codes and algorithm returned too long GA" as LMP_long_revision_reason
		from _pregnancies_long4
		union corr
		select distinct idxpren, dt_lmp, dt_indexprenatal, dt_gapreg,
				1 as LMP_long_revision, "Revised index or gapreg based on Z3A codes" as LMP_long_revision_reason
		from _pregnancies_long_prenatal3
		;
		quit;

	%*Join these back onto the original pregnancies;

	%*First, join together the LMP, indexprenatal, and gapreg information. This allows us to overwrite the dt_lmp, etc.
	values and more easily use them later with the getGA macro;
	proc sql;
		create table _pregnancies_long_int as
		select distinct a.idxpren, 
			case when missing(b.idxpren) then a.dt_lmp else b.dt_lmp end as dt_lmp format=MMDDYY10., 
			case when missing(b.idxpren) then a.dt_indexprenatal else b.dt_indexprenatal end as dt_indexprenatal format=MMDDYY10.,
			case when missing(b.idxpren) then a.dt_gapreg else b.dt_gapreg end as dt_gapreg format=MMDDYY10., 
			case when missing(b.idxpren) then 0 else 1 end as LMP_long_revision,
			case when missing(b.idxpren) then "Not too long" else LMP_long_revision_reason end as LMP_long_revision_reason
		from &INPUT_DATA as a
		left join pregnancies_long_rev as b
		on a.idxpren = b.idxpren
		;
		quit;

	%*Now, join that information onto all other information;
	proc sql;
		create table &output_data as
		select distinct a.*, b.dt_lmp, b.dt_indexprenatal, b.dt_gapreg, b.LMP_long_revision, b.LMP_long_revision_reason
		from &INPUT_DATA (drop = dt_lmp dt_indexprenatal dt_gapreg) as a
		left join _pregnancies_long_int as b
		on a.idxpren = b.idxpren
		;
		quit;

	%*Save a proc freq output of the revised dataset as an analytic dataset;
	proc freq data=&output_data noprint;
		table preg_outcome_clean / missing out=ana.preg_outc_post_long_&lmpindex._&round;
	run;

	%*Delete unnecessary datasets;
	proc datasets gennum = all;
		delete _pregnancies_long:;
	run;

%mend;













/********************************************************************************************************************************************

											03 - MACROS FOR INDEX PRENATAL TOO EARLY

Some pregnancies had their indexing prenatal encounter occur prior to their LMP OR too early in their pregnancy such that the pregnancy 
would not have been known at the time (e.g., 28 days). We revise the dates for these pregnancies.

********************************************************************************************************************************************/



/*
MACRO: revise_index_before_lmp
PURPOSE: The purpose of this macro is to revise the LMPs and index dates for those pregnancies that 
originally have their index date occurring prior to their estimated LMP.

Steps are detailed in the macro below and documentation.
%let INPUT_DATA=_toolong;
%let days=28;
*/

%macro revise_index_before_lmp(INPUT_DATA, output_data, LMPINDEX, DAYS=28, round=0);


	%*Subset to those pregnancies where the indexing prenatal occurs prior to the lmp AND
		those where the indexingprenatal encounter occurs <4 weeks after the LMP
	Used the revised values from the data cleaning above.
	Make sure to ignore those pregnancies WITHOUT indexing prenatal encounters;
	data _pregnancies_pre;
	set &INPUT_DATA;
		where . < dt_indexprenatal < dt_lmp + &DAYS; 
	run;

	
	%*Output counts of the pregnancy outcomes for later review;
	proc freq data=_pregnancies_pre noprint;
		table preg_outcome_clean / missing out=ana.preg_outc_early_idx_&lmpindex._&round;
	run;

	%*Output:
		(1) The number of pregnancies input into this macro
		(2) The number of pregnancies with an indexing prenatal PRIOR to the LMP
		(3) The number of pregnancies with an indexing prenatal after LMP but prior to 28 days gestation;
	proc sql noprint;
		select count(distinct idxpren) into :num_preg_input from &INPUT_DATA;
		select count(distinct idxpren) into :num_preg_revised from _pregnancies_pre;
		select count(distinct idxpren) into :num_preg_prior from _pregnancies_pre where . < dt_indexprenatal < dt_lmp;
		select count(distinct idxpren) into :num_preg_pre28 from _pregnancies_pre where dt_lmp <= dt_indexprenatal < dt_lmp + &DAYS;
		quit;
	%put Number of pregnancies overall: &num_preg_input;
	%put Number of pregnancies being revised in this macro (total): &num_preg_revised;
	%put Number of pregnancies with dt_indexprenatal < dt_lmp: &num_preg_prior;
	%put Number of pregnancies with dt_lmp <= dt_indexprenatal < dt_lmp + &DAYS: &num_preg_pre28;

	%***********
		STEP 1: Deal with those pregnancies with observed outcomes but no GA information
	***********;

	data _pregnancies_pre_outc;
	set _pregnancies_pre;

		%*Identify those pregnancies with observed outcomes but no GA information.;
		where lmp_algbased = 0 and preg_outcome_clean ne "UNK"; /*There should be no UNK pregnancies where lmp_algased = 0 but included for explicitness*/
		
		%*Revise the LMPs based upon maximums for each outcome.;
		if preg_outcome_clean in ("EM", "IAB", "UAB", "SAB") then dt_lmp_rev = dt_gapreg - 90;
			else if preg_outcome_clean in ("LBM","SB","LBS","UDL") then dt_lmp_rev = dt_gapreg - 301;

		%*Create an indicator as to whether the indexing prenatal encounter still occurs prior to the LMP;
		stillshort = (. < dt_indexprenatal < dt_lmp_rev + &DAYS);
	run;

	%*Output counts;
	proc sql noprint; 
		select count(distinct idxpren), sum(preg_outcome_clean = "LBS"), sum(preg_outcome_clean = "LBM"),
			sum(preg_outcome_clean = "SB"), sum(preg_outcome_clean = "MLS"), sum(preg_outcome_clean = "EM"),
			sum(preg_outcome_clean = "IAB"), sum(preg_outcome_clean = "SAB"), sum(preg_outcome_clean = "UAB"),
			sum(preg_outcome_clean = "UDL")
		into :num_preg_outc, :num_LBS, :num_LBM, :num_SB, :num_MLS, :num_EM, :num_IAB, :num_SAB,
				:num_UAB, :num_UDL
		from _pregnancies_pre_outc; 
		quit;

	%put Number of pregnancies with LMP assigned based upon outcome only: &num_preg_outc;
	%put Number of LBS: &num_LBS;
	%put Number of LBM: &num_LBM;
	%put Number of SB: &num_SB;
	%put Number of MLS: &num_MLS;
	%put Number of EM: &num_EM;
	%put Number of IAB: &num_IAB;
	%put Number of SAB: &num_SAB;
	%put Number of UAB: &num_UAB;
	%put Number of UDL: &num_UDL;

	%*Output those that are not too short;
	data _pregnancies_pre_step1;
	set _pregnancies_pre_outc;
		where stillshort = 0;
		dt_lmp = dt_lmp_rev;
		LMP_late_rev_reason = "LMP before index with revised GA for Outcomes";
		keep idxpren dt_indexprenatal dt_gapreg preg_outcome_clean dt_lmp LMP_late_rev_reason;
	run;
	
	%*Output counts;
	proc sql noprint; 
		select count(distinct idxpren), sum(preg_outcome_clean = "LBS"), sum(preg_outcome_clean = "LBM"),
			sum(preg_outcome_clean = "SB"), sum(preg_outcome_clean = "MLS"), sum(preg_outcome_clean = "EM"),
			sum(preg_outcome_clean = "IAB"), sum(preg_outcome_clean = "SAB"), sum(preg_outcome_clean = "UAB"),
			sum(preg_outcome_clean = "UDL")
		into :num_preg_outc, :num_LBS, :num_LBM, :num_SB, :num_MLS, :num_EM, :num_IAB, :num_SAB,
				:num_UAB, :num_UDL
		from _pregnancies_pre_step1; 
		quit;
	%put Number of pregnancies of sufficient GA length after max GA based on outcome: &num_preg_outc;
	%put Number of LBS: &num_LBS;
	%put Number of LBM: &num_LBM;
	%put Number of SB: &num_SB;
	%put Number of MLS: &num_MLS;
	%put Number of EM: &num_EM;
	%put Number of IAB: &num_IAB;
	%put Number of SAB: &num_SAB;
	%put Number of UAB: &num_UAB;
	%put Number of UDL: &num_UDL;

	%*For those pregnancies where the indexing prenatal encounter still occurs prior to their index date, we revise
	their indexing prenatal encounter date as the first prenatal encounter after their original LMP and prior to 
	the pregnancy outcome date.;

	%*Grab all of the prenatal encounter dates for these pregnancies;
	proc sql;
		create table _pregnancies_pre_pncidx as
		select distinct idxpren, preg_outcome_clean, dt_lmp, dt_gapreg, min(b.enc_date) as dt_indexprenatal,
			"LMP after index w revised GA by Outcome, Revised Index date" as LMP_late_rev_reason 
		from _pregnancies_pre_outc (where = (stillshort = 1)) as a
		left join temp.codeprenatal_meg1_dts as b
		on a.patient_deid = b.patient_deid and a.dt_lmp + &DAYS <= b.enc_date <= a.dt_gapreg
		group by idxpren, preg_outcome_clean, dt_lmp, dt_gapreg
		;
		quit;
	%*Output count;
	proc sql noprint; 
		select count(distinct idxpren), sum(preg_outcome_clean = "LBS"), sum(preg_outcome_clean = "LBM"),
			sum(preg_outcome_clean = "SB"), sum(preg_outcome_clean = "MLS"), sum(preg_outcome_clean = "EM"),
			sum(preg_outcome_clean = "IAB"), sum(preg_outcome_clean = "SAB"), sum(preg_outcome_clean = "UAB"),
			sum(preg_outcome_clean = "UDL")
		into :num_preg_outc, :num_LBS, :num_LBM, :num_SB, :num_MLS, :num_EM, :num_IAB, :num_SAB,
				:num_UAB, :num_UDL
		from _pregnancies_pre_pncidx; 
		quit;
	%put Number of pregnancies with defined outcome and revised indexing prenatal: &num_preg_outc;
	%put Number of LBS: &num_LBS;
	%put Number of LBM: &num_LBM;
	%put Number of SB: &num_SB;
	%put Number of MLS: &num_MLS;
	%put Number of EM: &num_EM;
	%put Number of IAB: &num_IAB;
	%put Number of SAB: &num_SAB;
	%put Number of UAB: &num_UAB;


	%************
		STEP 2: DEAL WITH THOSE WITH LMPS ASSIGNED BY GA CODES
	************;

	%*Identify those with the LMP assigned via the gestational age algorithm;
	data _pregnancies_pre_alg;
	set _pregnancies_pre;
		where lmp_algbased = 1;
	run;

	%*Output count;
	proc sql noprint; 
		select count(distinct idxpren), sum(preg_outcome_clean = "LBS"), sum(preg_outcome_clean = "LBM"),
			sum(preg_outcome_clean = "SB"), sum(preg_outcome_clean = "MLS"), sum(preg_outcome_clean = "EM"),
			sum(preg_outcome_clean = "IAB"), sum(preg_outcome_clean = "SAB"), sum(preg_outcome_clean = "UAB"),
			sum(preg_outcome_clean = "UDL"), sum(preg_outcome_clean = "UNK")
		into :num_preg_outc, :num_LBS, :num_LBM, :num_SB, :num_MLS, :num_EM, :num_IAB, :num_SAB,
				:num_UAB, :num_UDL, :num_UNK
		from _pregnancies_pre_alg; 
		quit;
	%put Number of pregnancies with LMP assigned via GA algorithm: &num_preg_outc;
	%put Number of LBS: &num_LBS;
	%put Number of LBM: &num_LBM;
	%put Number of SB: &num_SB;
	%put Number of MLS: &num_MLS;
	%put Number of EM: &num_EM;
	%put Number of IAB: &num_IAB;
	%put Number of SAB: &num_SAB;
	%put Number of UAB: &num_UAB;
	%put Number of UNK: &num_UNK;

	%*Implement the gestational age algorithm;
	options nomprint nomlogic nosymbolgen;
	%getGA(_pregnancies_pre_alg, ga_rev, GA=max_gest_age);
	options mprint mlogic symbolgen;

	%*Reassign the LMPs based on the revisions. If not revised, we will assume that they 
	only had Prenatal Care codes for GA assignment. This is part of the revised algorithm in MarketScan;
	data _pregnancies_pre_rev;
	set ga_rev;

		%*Assign the LMP so that it is not missing;
		if lmp_algbased = 1 then dt_lmp_rev = dt_lmp_alg;
			else if lmp_algbased = 0 and preg_outcome_clean ne "UNK" then dt_lmp_rev = dt_lmp_table; 
			%*Those assigned GAs based on PNC codes now have their GAs based on no codes, reassign their LMP
			by assuming they were &lmpindex days gestation at their indexing prenatal encounter, as originally done in the algorithm;
			else if lmp_algbased = 0 then dt_lmp_rev = dt_indexprenatal - &lmpindex;

		ga_length = dt_gapreg - dt_lmp_rev;
		toolong = ga_length > 307;
	run;

	%*Those pregnancies that are now too long with unknown outcomes, assign them a LMP of 28 at their indexing PNC.
	Of note, this only occurs to those with unknown outcomes.;
	proc sql;
		create table _pregnancies_pre_rev_long as
		select distinct idxpren, dt_indexprenatal, dt_indexprenatal - 28 as dt_lmp, dt_gapreg,
			"GA only by Prenatal Care Codes, Assumed 4w GA at index" as LMP_late_rev_reason
		from _pregnancies_pre_alg
		where idxpren in (select distinct idxpren from _pregnancies_pre_rev where toolong = 1)
		;
		quit;
	%*Output counts;
	proc sql noprint; select count(*) into :num_preg_toolong from _pregnancies_pre_rev_long; quit;
	%put Number of UNK pregnancies too long after revised GA algorithm: &num_preg_toolong;

	%*Review those pregnancies now with revised LMPs. Now need to deal with those where index prenatal is still prior to LMP;
	data _pregnancies_pre_rev2;
	set _pregnancies_pre_rev;
		where toolong = 0;
		indexearly = (. < dt_indexprenatal < dt_lmp_rev + &DAYS);
	run;

	%*Output those pregnancies where the index is not prior to the LMP anymore;
	proc sql;
		create table _pregnancies_pre_gaAlg as
		select distinct idxpren, dt_indexprenatal, dt_lmp_rev as dt_lmp, dt_gapreg, preg_outcome_clean,
			"Index after LMP with max_gest_age as gestational_age_days" as LMP_late_rev_reason
		from _pregnancies_pre_rev2
		where indexearly = 0;
		quit;
	%*Output count;
	proc sql noprint; 
		select count(distinct idxpren), sum(preg_outcome_clean = "LBS"), sum(preg_outcome_clean = "LBM"),
			sum(preg_outcome_clean = "SB"), sum(preg_outcome_clean = "MLS"), sum(preg_outcome_clean = "EM"),
			sum(preg_outcome_clean = "IAB"), sum(preg_outcome_clean = "SAB"), sum(preg_outcome_clean = "UAB"),
			sum(preg_outcome_clean = "UDL"), sum(preg_outcome_clean = "UNK")
		into :num_preg_outc, :num_LBS, :num_LBM, :num_SB, :num_MLS, :num_EM, :num_IAB, :num_SAB,
				:num_UAB, :num_UDL, :num_UNK
		from _pregnancies_pre_gaAlg; 
		quit;
	%put Number of pregnancies sufficiently early index prenatal after revised GA algorithm: &num_preg_outc;
	%put Number of LBS: &num_LBS;
	%put Number of LBM: &num_LBM;
	%put Number of SB: &num_SB;
	%put Number of MLS: &num_MLS;
	%put Number of EM: &num_EM;
	%put Number of IAB: &num_IAB;
	%put Number of SAB: &num_SAB;
	%put Number of UAB: &num_UAB;
	%put Number of UNK: &num_UNK;

	%*Assign those remaining pregnancies their minimum LMP according to the algorithm and see if any can be 
	considered finalized;
	data _pregnancies_pre_rev3;
	set _pregnancies_pre_rev2;
		where indexearly = 1;
		if dt_lmp_min ne . then dt_lmp_rev = dt_lmp_min;
			else dt_lmp_rev = dt_lmp_rev;

		indexearly2 = (. < dt_indexprenatal < dt_lmp_rev + &DAYS);
	run;

	%*Those with their LMPs prior to their indexing prenatal encounters are now finalized;
	proc sql;
		create table _pregnancies_pre_minLMP as
		select distinct idxpren, dt_indexprenatal, dt_lmp_rev as dt_lmp, dt_gapreg, preg_outcome_clean,
			"Index after LMP with max_gest_age as gestational_age_days and minimum LMP" as LMP_late_rev_reason
		from _pregnancies_pre_rev3
		where indexearly2 = 0
		;
		quit;
	%*Output count;
	proc sql noprint; 
		select count(distinct idxpren), sum(preg_outcome_clean = "LBS"), sum(preg_outcome_clean = "LBM"),
			sum(preg_outcome_clean = "SB"), sum(preg_outcome_clean = "MLS"), sum(preg_outcome_clean = "EM"),
			sum(preg_outcome_clean = "IAB"), sum(preg_outcome_clean = "SAB"), sum(preg_outcome_clean = "UAB"),
			sum(preg_outcome_clean = "UDL"), sum(preg_outcome_clean = "UNK")
		into :num_preg_outc, :num_LBS, :num_LBM, :num_SB, :num_MLS, :num_EM, :num_IAB, :num_SAB,
				:num_UAB, :num_UDL, :num_UNK
		from _pregnancies_pre_minLMP; 
		quit;
	%put Number of pregnancies sufficiently early index prenatal after using min LMP from revised GA algorithm: &num_preg_outc;
	%put Number of LBS: &num_LBS;
	%put Number of LBM: &num_LBM;
	%put Number of SB: &num_SB;
	%put Number of MLS: &num_MLS;
	%put Number of EM: &num_EM;
	%put Number of IAB: &num_IAB;
	%put Number of SAB: &num_SAB;
	%put Number of UAB: &num_UAB;
	%put Number of UNK: &num_UNK;


	%*Identify those that still have an indexing prenatal prior to their LMP. These will be linked to their 
	prenatal care encounters to determinen if they have recorded Z3A codes;
	data _pregnancies_pre_z3a;
	set _pregnancies_pre_rev3;
		where indexearly2 = 1;
	run;
	%*Output count;
	proc sql noprint; 
		select count(distinct idxpren), sum(preg_outcome_clean = "LBS"), sum(preg_outcome_clean = "LBM"),
			sum(preg_outcome_clean = "SB"), sum(preg_outcome_clean = "MLS"), sum(preg_outcome_clean = "EM"),
			sum(preg_outcome_clean = "IAB"), sum(preg_outcome_clean = "SAB"), sum(preg_outcome_clean = "UAB"),
			sum(preg_outcome_clean = "UDL"), sum(preg_outcome_clean = "UNK")
		into :num_preg_outc, :num_LBS, :num_LBM, :num_SB, :num_MLS, :num_EM, :num_IAB, :num_SAB,
				:num_UAB, :num_UDL, :num_UNK
		from _pregnancies_pre_z3a; 
		quit;
	%put Number of pregnancies evaluated for Specific gestational age codes: &num_preg_outc;
	%put Number of LBS: &num_LBS;
	%put Number of LBM: &num_LBM;
	%put Number of SB: &num_SB;
	%put Number of MLS: &num_MLS;
	%put Number of EM: &num_EM;
	%put Number of IAB: &num_IAB;
	%put Number of SAB: &num_SAB;
	%put Number of UAB: &num_UAB;
	%put Number of UNK: &num_UNK;


	%*Link them back with their gestational age encounters to determine which ones have Z3A codes;
	proc sql;
		create table _pregnancies_pre_pnc2 as
		select distinct a.*, count(distinct b.enc_date) as numZ3A
		from _pregnancies_pre_z3a as a
		left join (select * from temp.gestagepren where code_hierarchy = "Specific gestational age") as b
		on a.patient_deid = b.patient_deid and a.dt_indexprenatal <= b.enc_date <= a.dt_gapreg + 7
		group by idxpren
		;
		quit;

	%*Those with 0 Z3A codes will be assigned a gestational age by assuming that they were 4w0d
	gestation at their indexing prenatal encounter;
	proc sql;
		create table _pregnancies_pre_noZ3A as
		select distinct idxpren, dt_indexprenatal, dt_indexprenatal - 28 as dt_lmp, dt_gapreg, preg_outcome_clean,
			"No Specific GA codes, so assumed GA of 4W at indexing PNC" as LMP_late_rev_reason
		from _pregnancies_pre_pnc2
		where numZ3A = 0;
		quit;
	%*Output count;
	proc sql noprint; 
		select count(distinct idxpren), sum(preg_outcome_clean = "LBS"), sum(preg_outcome_clean = "LBM"),
			sum(preg_outcome_clean = "SB"), sum(preg_outcome_clean = "MLS"), sum(preg_outcome_clean = "EM"),
			sum(preg_outcome_clean = "IAB"), sum(preg_outcome_clean = "SAB"), sum(preg_outcome_clean = "UAB"),
			sum(preg_outcome_clean = "UDL"), sum(preg_outcome_clean = "UNK")
		into :num_preg_outc, :num_LBS, :num_LBM, :num_SB, :num_MLS, :num_EM, :num_IAB, :num_SAB,
				:num_UAB, :num_UDL, :num_UNK
		from _pregnancies_pre_noZ3A; 
		quit;
	%put Number of pregnancies evaluated for Specific gestational age codes but had none: &num_preg_outc;
	%put Number of LBS: &num_LBS;
	%put Number of LBM: &num_LBM;
	%put Number of SB: &num_SB;
	%put Number of MLS: &num_MLS;
	%put Number of EM: &num_EM;
	%put Number of IAB: &num_IAB;
	%put Number of SAB: &num_SAB;
	%put Number of UAB: &num_UAB;
	%put Number of UNK: &num_UNK;

	%*Those with at least 1 Z3A code. Grab all the prenatal encounter dates between their max LMP and
	the outcome date.;
	proc sql;
		create table _pregnancies_pre_Z3A as
		select distinct a.idxpren,  min(b.enc_date) as dt_indexprenatal, dt_lmp_alg as dt_lmp, a.dt_gapreg, preg_outcome_clean,
			"At least 1 Specific GA code, revised the date of the indexing PNC" as LMP_late_rev_reason
		from _pregnancies_pre_pnc2 (where = (numZ3A > 0)) as a
		left join temp.codeprenatal_meg1_dts as b
		on a.patient_deid=b.patient_deid and a.dt_lmp_alg+&days <= b.enc_date <= a.dt_gapreg /*CDL: ADDED +&days to min 1.24.2025 based on error identified by MGP*/
		group by idxpren, dt_indexprenatal, dt_lmp, dt_gapreg, preg_outcome_clean
		;
		quit;
	%*Output counts;
	proc sql noprint; 
		select count(distinct idxpren), sum(preg_outcome_clean = "LBS"), sum(preg_outcome_clean = "LBM"),
			sum(preg_outcome_clean = "SB"), sum(preg_outcome_clean = "MLS"), sum(preg_outcome_clean = "EM"),
			sum(preg_outcome_clean = "IAB"), sum(preg_outcome_clean = "SAB"), sum(preg_outcome_clean = "UAB"),
			sum(preg_outcome_clean = "UDL"), sum(preg_outcome_clean = "UNK")
		into :num_preg_outc, :num_LBS, :num_LBM, :num_SB, :num_MLS, :num_EM, :num_IAB, :num_SAB,
				:num_UAB, :num_UDL, :num_UNK
		from _pregnancies_pre_Z3A; 
		quit;
	%put Number of pregnancies w revised index dates because they had specific GA codes: &num_preg_outc;
	%put Number of LBS: &num_LBS;
	%put Number of LBM: &num_LBM;
	%put Number of SB: &num_SB;
	%put Number of MLS: &num_MLS;
	%put Number of EM: &num_EM;
	%put Number of IAB: &num_IAB;
	%put Number of SAB: &num_SAB;
	%put Number of UAB: &num_UAB;
	%put Number of UNK: &num_UNK;

	%*Create the final dataset of all the pregnancies with revised information;
	proc sql;
		create table _pregnancies_pre_final as
		select distinct * from _pregnancies_pre_step1 (drop = preg_outcome_clean)
		union corr
		select distinct * from _pregnancies_pre_pncidx (drop = preg_outcome_clean)
		union corr
		select distinct * from _pregnancies_pre_rev_long 
		union corr 
		select distinct * from _pregnancies_pre_gaAlg (drop = preg_outcome_clean)
		union corr
		select distinct * from _pregnancies_pre_minLMP (drop = preg_outcome_clean)
		union corr
		select distinct * from _pregnancies_pre_noZ3A (drop = preg_outcome_clean)
		union corr
		select distinct * from _pregnancies_pre_Z3A (drop = preg_outcome_clean)
		;
		quit;

	%*Add that information back to the LMPs;
	proc sql;
		create table _pregnancies_pre_all as
		select distinct a.idxpren,
			case when missing(b.idxpren) then a.dt_lmp else b.dt_lmp end as dt_lmp format=MMDDYY10.,
			case when missing(b.idxpren) then a.dt_indexprenatal else b.dt_indexprenatal end as dt_indexprenatal format=MMDDYY10.,
			case when missing(b.idxpren) then a.dt_gapreg else b.dt_gapreg end as dt_gapreg format=MMDDYY10.,
			case when missing(b.idxpren) then 0 else 1 end as LMP_late_rev,
			case when missing(b.idxpren) then "Not required" else b.LMP_late_rev_reason end as LMP_late_rev_reason
		from &INPUT_DATA (keep = idxpren dt_lmp dt_indexprenatal dt_gapreg) as a
		left join _pregnancies_pre_final as b
		on a.idxpren = b.idxpren
		;
		quit;

	%*Now, finally add it back to the original dataset;
	proc sql;
		create table &output_data as
		select distinct a.*, b.dt_lmp, b.dt_indexprenatal, b.dt_gapreg, b.LMP_late_rev, b.LMP_late_rev_reason
		from &INPUT_DATA (drop = dt_lmp dt_indexprenatal dt_gapreg) as a
		left join _pregnancies_pre_all as b
		on a.idxpren = b.idxpren
		;
		quit;

	%*Output counts of the final pregnancy outcomes after revising based on early indices;
	proc freq data=&output_data noprint;
		table preg_outcome_clean / missing out=ana.preg_outc_post_early_idx_&lmpindex._&round;
	run;

	%*Delete the unnecessary datasets;
	proc datasets gennum=all;
		delete _pregnancies_pre:;
	run;


%mend;

















/********************************************************************************************************************************************

															04 - MACROS FOR OVERLAP

Some pregnancies were identified as overlapping with other pregnancies based upon their LMP and outcome date. This has to be resolved.
The macros below deal with this.

********************************************************************************************************************************************/




/*
MACRO: implement_hierarchy
PURPOSE: To implement a hierarchy that is regularly applied within the algorithm:
(1) If there is only one remaining pregnancy of the pair, then output that pregnancy
(2) If there are two remaining pregnancies in the pair, then output the one with the highest hierarchy
(2) If same hierarchy, then output based on presence of Specific gestational age codes
(3) If no specific gestational age codes, then output the first.

These are ordered according to dt_LMP
*/

%macro implement_hierarchy(input, output);

	%*Determine what is the minimum value and subset to those. Further, count the revised number of pregnancies within the variable COUNT;
	proc sql;
		create table _hier1 as
		select distinct *, count(*) as num_preg_rev
		from  (
			select *, min(hier) as min_hier
			from &input
			group by patient_deid, enrolid, count
			having hier = min_hier
			)
		group by patient_deid, enrolid, count
		;
		quit;

	%*Sort through the hierarchies.;
	proc sort data=_hier1;
		by patient_deid enrolid count dt_lmp;
	run;
	data &output;
	set _hier1;
		by patient_deid enrolid count dt_lmp;

		%*If only one pregnancy, then retain that one.;
		if first.count and num_preg_rev = 1 then output;
		
		%*If NEITHER (num_preg_w_Z3A=0) or BOTH (num_pref_w_Z3A=2) of the pair has a Z3A code, then output the first;
		else if first.count and num_preg_w_Z3A ne 1 then output;

		%*If the first of the pair where only 1 has a Z3A code and this has a Z3A code then output;
		else if first.count and num_preg_w_Z3A = 1 and numZ3A > 0 then output;
		%*Otherwise, output the other pregnancy in the pair that has a Z3A code;
		else if num_preg_w_Z3A = 1 and numZ3A > 0 then output;

		drop hier min_hier num_preg_rev;
	run;

%mend;







/*
MACRO: clean_overlap
PURPOSE: to address those pregnancies that overlap with one another within a person

INPUTS:
INPUT_DATA -- Input pregnancy dataset
OUTPUT_DATA -- Name for the output pregnancy dataset
*/

%macro clean_overlap(INPUT_DATA, OUTPUT_DATA, LMPINDEX);

	/*Testing:
	%let INPUT_DATA = _before_index;
	*/
	
	%let num_preg_overlap = 1;
	%let num = 0;
	
	%do %until(&num_preg_overlap = 0); %*Repeat this macro until there are no overlapping pregnancies.;
	
		%let num = %eval(&num +1); 
		%put Round: &num;
		
		%if &num=1 %then %do;
			data _pregnancies;
			set &input_data;
			run;
		%end;
		%else %do;
			data _pregnancies;
			set _pregnancies_final_%eval(&num-1);
			run;
		%end;

		%*************
			STEP 1: Identify those pregnancies that have any amount of overlap.;
	
		%*Create a variable COUNT that identifies those pregnancies that overlap within a person;
		proc sort data=_pregnancies;
			by enrolid dt_lmp;
		run;
		data _pregnancies2;
		set _pregnancies;
			format last_end MMDDYY10.;
			by enrolid dt_lmp;
			retain count;
	
			last_id = lag1(enrolid);
			last_end = lag1(dt_gapreg);
			last_lmp = lag1(dt_lmp);
	
			if first.enrolid then do;
				count = 1;
				overlap = .;
				lmp_overlap = .;
			end;
				else if enrolid = last_id and dt_lmp <= last_end then do;
					count = count;
					overlap = abs(last_end - dt_lmp);
					lmp_overlap = abs(last_lmp - dt_lmp);
				end;
				else do;
					count = count + 1;
					overlap = .;
					lmp_overlap = .;
				end;
		  run;

		%*Create the variable num_preg which counts the number of pregnancies within each COUNT value. This will be used
		to identify the number of pregnancies that overlap;
		proc sql;
			create table _pregnancies3 as
			select distinct patient_deid, enrolid, idxpren, dt_lmp, dt_gapreg, count, count(*) as num_preg
			from _pregnancies2
			group by enrolid, count
			;
			quit;
	
		%*Output counts;
		proc sql noprint;
			select count(*) into :num_preg_orig from _pregnancies;
			select count(*) into :num_preg_nooverlap from _pregnancies3 where num_preg = 1;
			select count(*) into :num_preg_overlap from _pregnancies3 where num_preg > 1;
			quit;
		%put Number of pregnancies in the original cohort: &num_preg_orig;
		%put Number of pregnancies with no overlap: &num_preg_nooverlap;
		%put Number of pregnancies with overlap: &num_preg_overlap;		
		
		%*Subset to overlapping pregnancies.;
		proc sql;
			create table _preg_overlap0 as
			select distinct * 
			from _pregnancies2 
			where idxpren in (select distinct idxpren from _pregnancies3 where num_preg > 1)
			;
			quit;
			
		%*Count the number of pregnancies within each series of overlaps;
		proc sql;
			create table _preg_overlap0_2 as
			select distinct a.*, b.num_preg
			from _preg_overlap0 as a
			left join _pregnancies3 as b
			on a.enrolid=b.enrolid and a.idxpren=b.idxpren
			;
			quit;
			
		%***Subset to the FIRST two overlapping pregnancies.;

		%*Create a variable - ROW - which indicates the order of the pregnancy in time BASED UPON LMP.;
		proc sort data=_preg_overlap0_2;
			by enrolid count dt_lmp;
		run;
		data _preg_loop2;
		set _preg_overlap0_2;
			by enrolid count dt_lmp;
			retain row;

			if first.count then row = 1;
				else row = row + 1;
		run;
		%*Subset to the first two pregnancies in each pair. Those excluded will be restacked later.;
		data _preg_overlap;
		set _preg_loop2;
			where row <= 2;
		run;
		%*Output a dataset of the overlapping pregnancy outcomes;
		proc freq data=_preg_overlap noprint;
			table preg_outcome_clean / missing out=ana.preg_outc_overlap_&lmpindex._&num;
		run;


		%*Grab the pregnancies that are not going to be evaluated in this round. They will be stacked back later.;
		proc sql;
			create table _no_pair_eval as
			select * 
			from _pregnancies
			where idxpren not in (select distinct idxpren from _preg_overlap)
			;
			quit;
		
		%*Create an indicator for the quality of evidence for the overlapping pregnancies outcome;
		data _preg_overlap2;
		set _preg_overlap;
	
			if Preg_CodeType = 'PR' then outc_quality = 1;
				else if preg_outcome_clean in ("IAB", "SAB", "UAB") and Preg_CodeType = "RX" then outc_quality = 2;
				else if preg_outcome_clean ne "UNK" then outc_quality = 3;
				else outc_quality = 4;
		run;		

		%*Create an indicator as to whether there is a Z3A code.;
		proc sql;
			create table _preg_overlap3 as
			select distinct aa.*, bb.numZ3A
			from _preg_overlap2 as aa
			left join (select distinct a.patient_deid, a.idxpren, count(distinct enc_date) as numZ3A
						from _preg_overlap as a
						left join temp.gestagepren (where = (code_hierarchy = "Specific gestational age")) as b
						on a.patient_deid=b.patient_deid and a.dt_indexprenatal <= b.enc_date <= a.dt_GApreg
						group by a.patient_deid, a.idxpren) as bb
			on aa.patient_deid=bb.patient_deid and aa.idxpren=bb.idxpren
			;
			quit;

			

		%*Collect only the relevant information for pairwise evaluation;
		proc sql;
			create table _pair_eval_summary as
			select distinct patient_deid, enrolid, idxpren, dt_indexprenatal, dt_gapreg, dt_lmp, count, preg_outcome_clean, 
					num_preg, outc_quality, numZ3A,
					sum(numZ3A > 0) as num_preg_w_Z3A, max(lmp_overlap) as max_lmp_overlap
			from _preg_overlap3 
			group by enrolid, count
			;
			quit;

		%******
			STEP 2: Implement hierarchies for all the overlapping pregnancies.;

		%*Output data to assess what pregnancies were compared in this group;
		%*Get counts of the pairwise comparisons
		We want the outcomes compared, the outcome quality, and the numZ3A codes;
		proc sort data=_pair_eval_summary out=pairwise;
			by enrolid count dt_lmp; /*NOTE: now sorted by LMP*/
		run;
		data ana.pairwise_outc_round&num._&lmpindex ;
		    set pairwise;
		    by enrolid count dt_lmp;
		    retain preg_outcome_rolled outc_info_rolled num_specificGA_rolled;
		    
		    /* Initialize the combined variable at the first occurrence of each group */
		    if first.count then do;
				preg_outcome_rolled = preg_outcome_clean;
				outc_info_rolled = strip(put(outc_quality, best.));
				num_specificGA_rolled = strip(put(numZ3A, best.));;
			end;
		    else do;
				preg_outcome_rolled = catx('-', preg_outcome_rolled, preg_outcome_clean);
				outc_info_rolled = catx('-', outc_info_rolled, outc_quality);
				num_specificGA_rolled = catx('-', num_specificGA_rolled, numZ3A);
			end;
		    
		    /* Output only at the last occurrence of each group */
		    if last.count then output;
		    
		    /* Drop intermediate variables */
		    keep enrolid count preg_outcome_rolled outc_info_rolled num_specificGA_rolled; 
		run;


		%*First, define the necessary variables to implement the hierarchies.;
		%*We need to know (in hierarchical order):
			(1) If at least 1 pregnancy in the pair has a procedure code
			(2) If only has has a Specific gestational age code - already derived
			(3) If at least 1 pregnancy is a delivery based on diagnosis codes
			(4) If at least 1 pregnancy is an abortion based on RX code
			(5) If at least 1 pregnancy is defined via diagnosis codes
			(6) If both have UNK outcomes - rederive
			;
		proc sql;
			create table _pair_eval_nolmp2 as
			select distinct *, sum(outc_quality = 1) as num_preg_w_proc, sum(preg_outcome_clean in ('LBS','LBM','SB','UDL') and outc_quality = 3) as num_preg_w_deldx,
				sum(preg_outcome_clean in ('SAB','IAB','UAB') and outc_quality = 2) as num_preg_w_aborrx, sum(outc_quality = 3) as num_preg_w_dx,
				sum(preg_outcome_clean = 'UNK') as num_preg_w_unk, sum(preg_outcome_clean in ('LBS','LBM','SB','UDL')) as num_del
			from _pair_eval_summary
			group by patient_deid, enrolid, count
			;
			quit;

		%*First, evaluate if those with two deliveries are sufficiently far apart;
		proc sort data=_pair_eval_nolmp2;
			by patient_deid enrolid count dt_gapreg;
		run;
		data _pair_eval_nolmp2_2;
		set _pair_eval_nolmp2;
			by patient_deid enrolid count dt_gapreg;

			ga_preg_diff = dt_gapreg - lag1(dt_gapreg);

			if first.count then ga_preg_diff = .;
		run;

		%*Create an indicator of how far apart on each pregnancy;
		proc sql;
			create table _pair_eval_nolmp2_3 as
			select distinct *, max(ga_preg_diff) as overall_ga_preg_diff
			from _pair_eval_nolmp2_2
			group by enrolid, count;
			quit;


		%*Now split the pairs into groups accordingly;
		data _pair_eval_nolmp3;
		set _pair_eval_nolmp2_3;
			length group $40;

			if num_del = 2 and overall_ga_preg_diff >= 168 then group = "Two deliveries";
				else if num_preg_w_proc = 1 then group = "One procedure code";
				else if num_preg_w_proc = 2 then group ="Two procedure codes";
				else if num_preg_w_Z3A = 1 then group = "One specific GA code";
				else if num_preg_w_deldx > 0 then group = "Delivery diagnosis code";
				else if num_preg_w_aborrx > 0 then group = "Abortion prescription fill";
				else if num_preg_w_dx > 0 then group = "Diagnosis code";
				else if num_preg_w_unk = 2 then group = "Unknown outcome";
		run;

		

		%*Delete unnecessary datasets for working memory;
		proc datasets gennum=all;
			delete _pair_eval_nolmp _pair_eval_nolmp2 _pair_eval_nolmp2_3;
		run;

		%**********
			Deal with those where there are two deliveries and both are biologically plausible;

		proc sort data=_pair_eval_nolmp3; 
			by enrolid count dt_gapreg;
		run;
		data pair_eval_2delp;
		set _pair_eval_nolmp3;
			where group = "Two deliveries";
			by enrolid count dt_gapreg;

			last_outc = lag1(dt_gapreg);

			if last.count then dt_lmp = last_outc + 1;
				else dt_lmp = dt_lmp;

			drop last_outc;
		run;
		%*Output counts;
		proc sql noprint;
			select count(*) into :num_preg from _pair_eval_nolmp3 where group = "Two deliveries";
			quit;
		%put Number of pregnancies with 2 deliveries that are plausibly far apart: &num_preg;


		%**********
			Deal with those where only one is defined via procedure codes.;

		%*Retain only the pregnancy where the outcome is defined with a procedure code.;
		data pair_eval_proc;
		set _pair_eval_nolmp3;
			where group = "One procedure code";
			if outc_quality = 1 then output;
		run;
		%*Output counts;
		proc sql noprint;
			select count(*) into :num_preg from _pair_eval_nolmp3 where group = "One procedure code";
			select count(*) into :num_preg_proc from pair_eval_proc;
			quit;
		%put Number of pregnancies where only 1 in the pair had an outcome defined via procedure code in round &num : &num_preg;
		%put Number of pregnancies after retaining the only preg with a procedure code in round &num : &num_preg_proc;

		%*************
			Deal with those where both pregnancy outcomes are defined via procedure codes;
		data _pair_eval_2proc1;
		set _pair_eval_nolmp3;
			where group = "Two procedure codes";
		run;

		%*Implement an outcome hierarchy to determine which to keep: (1) LBS LBM SB MLS, (2) UDL, (3) SAB (4) IAB (5) UAB, (6) EM;
		data _pair_eval_2proc2;
		set _pair_eval_2proc1;
			if preg_outcome_clean in ('LBS','LBM','SB','MLS') then hier = 1;
				else if preg_outcome_clean = 'UDL' then hier = 2;
				else if preg_outcome_clean = 'SAB' then hier = 3;
				else if preg_outcome_clean = 'IAB' then hier = 4;
				else if preg_outcome_clean = 'UAB' then hier = 5;
				else if preg_outcome_clean = 'EM' then hier = 6;
				else hier = 7; /*None should fall in here -- CDL: Changed from 6 to 7 on 1.28.2025*/
		run;

		%implement_hierarchy(input=_pair_eval_2proc2, output=pair_eval_2proc);

		%*Output counts;
		proc sql noprint;
			select count(*) into :num_preg from _pair_eval_2proc1;
			select count(*) into :num_preg_2proc from pair_eval_2proc;
			quit;
		%put Number of pregs in pairs with both had a procedure code in round &num : &num_preg;
		%put Number of pregs after evaluating procedure code pairs in round &num : &num_preg_2proc;

		%*Delete unnecessary datasets;
		proc datasets gennum = all;
			delete _pair_eval_2proc:;
		run;


		%*Deal with those where only one pregnancy has a specific GA code;

		%*Retain only the pregnancy with the specific GA code;
		data pair_eval_1ga;
		set _pair_eval_nolmp3;
			where group = "One specific GA code";
			if numZ3A > 0 then output;
		run;

		%*Output counts;
		proc sql noprint;
			select count(*) into :num_preg from _pair_eval_nolmp3 where group = "One specific GA code";
			select count(*) into :num_preg_1ga from pair_eval_1ga;
			quit;
		%put Number of pregnancies in pairs where one has a specific GA code in round &num : &num_preg;
		%put Number of pregnancies after evaluating pairs with one specific GA code in round &num : &num_preg_1ga;


		%************
			Deal with those where at least 1 pregnancy is a delivery defined via dx codes;

		%*Subset to the relevant population;
		data _pair_eval_deldx1a;
		set _pair_eval_nolmp3;
			where group = "Delivery diagnosis code";
		run;

		%*If only one of the pregnancies is a delivery defined via dx codes, then retain that one;
		data _pair_eval_deldx1;
		set _pair_eval_deldx1a;
			if num_preg_w_deldx = 1 and outc_quality = 3 and preg_outcome_clean in ('LBM','LBS','SB','MLS', 'UDL') then output;
				else if num_preg_w_deldx > 1 then output;
		run;

		%*Implement an outcome hierarchy to determine which to keep: (1) LBS LBM, (2) SB MLS, (3) UDL, (4) Other;
		data _pair_eval_deldx2;
		set _pair_eval_deldx1;
			if preg_outcome_clean in ('LBS','LBM') then hier = 1;
				else if preg_outcome_clean in ('SB','MLS') then hier = 2;
				else if preg_outcome_clean = 'UDL' then hier = 3;
				else hier = 4;
		run;

		%implement_hierarchy(input=_pair_eval_deldx2, output=pair_eval_deldx);

		%*Output counts;
		proc sql noprint;
			select count(*) into :num_preg from _pair_eval_deldx1a;
			select count(*) into :num_preg_deldx from pair_eval_deldx;
			quit;
		%put Number of pregnancies in pairs where at least one was a delivery with a dx code in round &num : &num_preg;
		%put Number of pregnancies after evaluating pairs where at least one was a delivery with dx code in round &num : &num_preg_deldx;

		%*Delete unnecessary datasets;
		proc datasets gennum=all;
			delete _pair_eval_deldx:;
		run;


		%************
			Deal with those where at least 1 pregnancy is an abortion defined via rx;

		data _pair_eval_aborrx1;
		set _pair_eval_nolmp3;
			where group = "Abortion prescription fill";
		run;

		%*If only one of the pregnancies is an abortion defined via rx fill, then output it.;
		data _pair_eval_aborrx2;
		set _pair_eval_aborrx1;
			if num_preg_w_aborrx = 1 and outc_quality = 2 and preg_outcome_clean in ('SAB','IAB','UAB') then output;

			else if num_preg_w_aborrx ne 1 then output;
		run;

		%*Create a hierarchy variable;
		data _pair_eval_aborrx3;
		set _pair_eval_aborrx2;
			if outc_quality = 2 and preg_outcome_clean = 'SAB' then hier = 1;
				else if outc_quality = 2 and preg_outcome_clean = 'IAB' then hier = 2;
				else if outc_quality = 2 and preg_outcome_clean = 'UAB' then hier = 3;
				else hier = 4;
		run;

		%implement_hierarchy(input=_pair_eval_aborrx3, output=pair_eval_aborrx);

		%*Output counts;
		proc sql noprint;
			select count(*) into :num_preg from _pair_eval_nolmp3 where group = "Abortion prescription fill";
			select count(*) into :num_preg_aborrx from pair_eval_aborrx;
			quit;
		%put Number of pregnancies in pairs where at least 1 is an abortion with rx fill in round &num : &num_preg;
		%put Number of pregnancies after evaluating pairs where at least 1 is an abortion with rx fill in round &num : &num_preg_aborrx;

		%*Delete unnecessary datasets;
		proc datasets gennum=all;
			delete _pair_eval_aborrx:;
		run;


		%***************
			Deal with those where at least 1 pregnancy has an outcome defined via a diagnosis code;
		data _pair_eval_dx1;
		set _pair_eval_nolmp3;
			where group = "Diagnosis code";
		run;

		%*If only one pregnancy is an outcome defined via DX codes then retain that one;
		data _pair_eval_dx2;
		set _pair_eval_dx1;
			
			if num_preg_w_dx = 1 and outc_quality = 3 then output;
			
			else if num_preg_w_dx ne 1 then output;
		run;

		%*Now create a hierarchy variable to reflect the decisions in our data cleaning;
		data _pair_eval_dx3;
		set _pair_eval_dx2;
			if outc_quality = 3 and preg_outcome_clean = 'SAB' then hier=1;
				else if outc_quality = 3 and preg_outcome_clean = 'IAB' then hier=2;
				else if outc_quality = 3 and preg_outcome_clean = 'UAB' then hier = 3;
				else if outc_quality = 3 and preg_outcome_clean = 'EM' then hier = 4;
				else hier = 6;
		run;

		%implement_hierarchy(input=_pair_eval_dx3, output=pair_eval_dx);

		%*Output counts;
		proc sql noprint;
			select count(*) into :num_preg from _pair_eval_nolmp3 where group = "Diagnosis code";
			select count(*) into :num_preg_aborrx from pair_eval_dx;
			quit;
		%put Number of pregnancies in pairs where at least 1 is an outcome with dx code in round &num : &num_preg;
		%put Number of pregnancies after evaluating pairs where at least 1 is an outcome with dx code in round &num : &num_preg_aborrx;

		%*Delete unnecessary datasets;
		proc datasets gennum=all;
			delete _pair_eval_dx:;
		run;


		%*****************
			Deal with those where both pregnancies are unknown outcomes;
		data _pair_eval_unk1;
		set _pair_eval_nolmp3;
			where group = "Unknown outcome";
		run;

		%*All have the same hierarchy, but this allows us to implement the same hierarchy;
		data _pair_eval_unk2;
		set _pair_eval_unk1;
			hier = 1;
		run;

		%implement_hierarchy(input=_pair_eval_unk2, output=pair_eval_unk);

		%*Output counts;
		proc sql noprint;
			select count(*) into :num_preg from _pair_eval_nolmp3 where group = "Unknown outcome";
			select count(*) into :num_preg_aborrx from pair_eval_unk;
			quit;
		%put Number of pregnancies in pairs where both are UNK outcomes in round &num : &num_preg;
		%put Number of pregnancies after evaluating pairs where both are UNK outcomes in round &num : &num_preg_aborrx;

		%*Delete unnecessary datasets;
		proc datasets gennum=all noprint;
			delete _pair_eval_unk:;
		run;




		%**************
			Do the final evaluation of those that were evaluated;

		%******Now stack all of the pregnancies;
		proc sql;
			create table _stack1 as
			select distinct *, "Two deliveries" as clean from pair_eval_2delp
			union corr
			select distinct *, "One procedure code" as clean from pair_eval_proc
			union corr
			select distinct *, "Two procedure codes" as clean from pair_eval_2proc
			union corr
			select distinct *, "One specific GA code" as clean from pair_eval_1ga
			union corr
			select distinct *, "Delivery via diagnosis code" as clean from pair_eval_deldx
			union corr
			select distinct *, "Abortion via prescription fill" as clean from pair_eval_aborrx
			union corr
			select distinct *, "At least 1 outcome via diagnosis code" as clean from pair_eval_dx
			union corr
			select distinct *, "Both outcomes unknown" as clean from pair_eval_unk
			;
			quit;

		%*Get a new index prenatal encounter date based upon the LMP of the retained pregnanacy. This indexing prenatal encounter MUST occur
		on or after 28 days of gestation.;
		proc sql;
			create table _stack1_revidx as
			select distinct a.patient_deid, a.enrolid, a.idxpren, a.group, a.dt_lmp as dt_lmp_rev, min(b.enc_date) as dt_indexprenatal_rev format=MMDDYY10.
			from _stack1 as a
			left join temp.codeprenatal_meg1_dts as b
			on a.patient_deid=b.patient_deid and /*a.dt_indexprenatal*/ a.dt_lmp + 28 <= b.enc_date <= a.dt_gapreg
			group by a.enrolid, a.patient_deid, a.idxpren, a.group, a.dt_lmp
			;
			quit; 

		%*Now output the pregnancies within the overlapping ones and revise the dates of the indexing prenatal encounters;
		proc sql;
			create table _preg_stack1 as 
			select distinct a.*, b.dt_lmp_rev format=MMDDYY10. length=8, 
					b.dt_indexprenatal_rev format=MMDDYY10. length=8,
					b.group
			from &input_data as a
			inner join _stack1_revidx as b
			on a.enrolid=b.enrolid and a.idxpren=b.idxpren
			;
			quit;

		%*Now output the revised indexing prenatal date estimate;
		data _preg_stack1_2;
		set _preg_stack1;
			dt_indexprenatal = dt_indexprenatal_rev;
			dt_lmp = dt_lmp_rev; %*dt_lmp_rev will only be different for deliveries;

			drop dt_indexprenatal_rev dt_lmp_rev;
		run;

		%*Identify only those pregnancies that did NOT have two deliveries a plausible distance apart;
		proc sql;
			create table _preg_stack1_3 as 
			select * 
			from _preg_stack1_2
			where idxpren not in (select distinct idxpren from _stack1 where clean = "Two deliveries")
			;
			quit;

			
		%*Re-run the gestational age algorithm to ensure that the GA is assigned based upon all encounters around it.
		This will only affect pregnancies with unknown outcomes, so we only retain those.;
		options nomprint nomlogic nosymbolgen;
		%getga(input_data = _preg_stack1_3, output_data = _ga_rev_int, GA = gestational_age_days);
		options mprint mlogic symbolgen;

		data _pregnancies_ga_rev;
		set _ga_rev_int;
			%*Assign the LMP so that it is not missing;
			if lmp_algbased = 1 then dt_lmp = dt_lmp_alg;
				else if lmp_algbased = 0 and preg_outcome_clean ne "UNK" then dt_lmp = Dt_LMP_table;
				%*Those assigned GAs based on PNC codes now have their GAs based on no codes, reassign their LMP
				by assuming they were lmpindex days gestation at their indexing prenatal encounter;
				else if lmp_algbased = 0 then dt_lmp = dt_indexprenatal - &lmpindex;
		run;

		%*Then re-do the other algorithms on them, as they should stand on their own.;
		%revise_too_long(_pregnancies_ga_rev, _long, &lmpindex, days=307, round=&num);
		/*Output; _pregnancies2_rev*/ 

		/*Some pregnancies have their indexing prenatal encounter prior to their LMP or too early in gestation (we define as 28 days). 
		We will revise these pregnancies.*/
		%revise_index_before_lmp(_long, _before, &lmpindex, days=28, round=&num);

		%*Retain all the information from the pregnancies except for the information on LMP and index prenatal encounter date;
		proc sql;
			create table _evaluated_&num as
			select distinct a.patient_deid, a.enrolid, a.idxpren, 
					case when missing(c.idxpren) then b.dt_indexprenatal else c.dt_indexprenatal end as dt_indexprenatal format=MMDDYY10. length=8,
					/*CDL: MODIFIED becuase it was not carrying over the updated outcome date - 2.21.2025*/
/*					a.dt_gapreg, */
					case when missing(c.idxpren) then b.dt_gapreg else c.dt_gapreg end as dt_gapreg format=MMDDYY10. length=8,
					case when missing(c.idxpren) then b.dt_lmp else c.dt_lmp end as dt_lmp format=MMDDYY10. length=8,
					a.count, a.preg_outcome_clean,
					a.num_preg, a.outc_quality, a.numZ3A, a.num_preg_w_Z3A, a.max_lmp_overlap, a.num_preg_w_proc, a.num_preg_w_deldx,
					a.num_preg_w_aborrx, a.num_preg_w_dx, a.num_preg_w_unk, a.group, a.clean
			from _stack1 as a
			left join _before as b
			on a.patient_deid=b.patient_deid and a.idxpren=b.idxpren
			left join (select distinct * from _preg_stack1_2 where group = "Two deliveries") as c
/*			left join (select distinct * from _preg_stack1_2 where idxpren in (select distinct idxpren from _stack1 where clean = "Two deliveries")) as c*/
			on a.patient_deid=c.patient_deid and a.idxpren=c.idxpren
			;
			quit;

		%*Now union the rest of the pregnancies onto the dataset and create the dataset that we want for the next
		round of evaluation;
		proc sql;
			create table _pregnancies_final_&num as
			select distinct * from _evaluated_&num (drop = clean)
			union corr
			select distinct * from _no_pair_eval
			;
			quit;

	%end;


	%*The final set of pregnancies comes from _pregnancies_final_%eval(&num-1);
	data _final_pregnancies;
	set _pregnancies_final_%eval(&num-1);
	run;

	%*Output a dataset of the revised overlapping pregnancy outcomes;
	proc freq data=_final_pregnancies noprint;
		table preg_outcome_clean / missing out=ana.preg_outc_post_overlap_&lmpindex;
	run;

	%*Create a reference set for how each was determined;
	proc sql;
		create table _reference1 as
		%do i=1 %to %eval(&num-1);
			%if &i = 1 %then %do;
				select distinct *, &i as round from _evaluated_&i
			%end;
			%else %do;
				union corr
				select distinct *, &i as round from _evaluated_&i
			%end;
		%end;
		;
		quit;


	%*********
		STEP 3: Clean evaluated pregnancies;

	%*Sort it by pregnancy and round;
	proc sort data=_reference1;
		by enrolid idxpren round;
	run;
	%*Output the last round-s value;
	data _reference2;
	set _reference1;
		by enrolid idxpren round;
		if last.idxpren then output;
	run;

	%*Append the reason for why the pregnancy was retained when compared to pairs around it;
	proc sql;
		create table _final2 as
		select distinct a.*, b.clean as overlap_reason_retained
		from _final_pregnancies as a
		left join _reference2 as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		;
		quit;


	%*********
		STEP 4: Output the final pregnancy dataset;

	%*First, inner join the overlapping pregnancies so that we are subset appropriately. Retain the revised values of dt_indexprenatal,
		dt_lmp, dt_gapreg, and preg_outcome_clean. Further, retain the reason why the pregnanacy was retained.;
	proc sql;
		create table &OUTPUT_DATA as
		select distinct a.*, b.dt_indexprenatal as dt_indexprenatal format=MMDDYY10. length=8, b.dt_lmp as dt_lmp format=MMDDYY10. length=8, 
			b.dt_gapreg as dt_gapreg format=MMDDYY10. length=8,
			b.preg_outcome_clean as preg_outcome_clean format=$10. length=10, 
			b.overlap_reason_retained
		from &INPUT_DATA (drop = dt_indexprenatal dt_lmp dt_gapreg preg_outcome_clean) as a
		inner join _final2 as b
		on a.enrolid=b.enrolid and a.idxpren=b.idxpren
		;
		quit;


%mend;













/********************************************************************************************************************************************

															05 - CLEAN UP PREGNANCIES

The final step is to clean up the pregnancies dataset. This is done by applying a macro that uses all of the above macros.

********************************************************************************************************************************************/






*We are interested in identifying cohorts of pregnancies where we make different assumptions about the gestational
age at which someone enters the dataset. Here, we consider two assumptions: 6*7 or 42 and 9*7 or 63
%let lmpindex = 63;

%macro reassign_lmps(lmpindex);

	/**We want to identify those pregnancies with UNK outcomes that were assigned 140d for their gestational length */
	*For these pregnancies, assign their LMP by subtracting 6*7 or 42 from the date of their indexing prenatal encounter;
	data _pregnancies2;
	set temp.pregnancies;

		%*For those pregnancies with unknown outcomes and no gestational age codes, reassign their LMP as occurring 
		at the prespecified gestational age at the index date (&lmpindex);
		if preg_outcome_clean = 'UNK' and lmp_algbased = 0 then dt_LMP = DT_INDEXPRENATAL - &lmpindex;
			else dt_LMP = dt_LMP;

		%*Create indicators for the original values of dt_lmp, indexprenatal and gapreg. These will be overwritten for 
		pregnancies that we clean.;
		dt_lmp_orig = dt_lmp;
		dt_indexprenatal_orig = dt_indexprenatal;
		dt_gapreg_orig = dt_gapreg;
	run;

		
	/*Some pregnancies with unknown outcomes have too long of gestational lengths, which we define as
	> 43w6d or 307 days.*/
	%revise_too_long(_pregnancies2, _toolong, &lmpindex, days=307, round=0);
	/*Output; _toolong*/ 

	data test;
	set _toolong;
		where . < dt_indexprenatal < dt_lmp + 28;
	run;

	proc freq data=test;
	table preg_outcome_clean/missing;
	run;

	data test2;
	set test;
		where LMP_AlgBased = 0 and preg_outcome_clean ne "UNK";
	run;

	proc freq data=test2;
		table preg_outcome_clean / missing;
	run;

	/*Some pregnancies have their indexing prenatal encounter prior to their LMP or too early in gestation (we define as 28 days). 
	We will revise these pregnancies.*/
	%revise_index_before_lmp(_toolong, _before_index, &lmpindex, days=28, round=0);
	

	%*Some pregnancies overlap with others. That is dealt with below;
	%clean_overlap(INPUT_DATA=_before_index, OUTPUT_DATA=_pregnancies_clean, LMPINDEX=&lmpindex);

	*Output the final pregnancy cohort. Subset to those pregnancies with
	LMPs within the relevant window;
	data ana.preg_cohort_&lmpindex;
	set _pregnancies_clean; 
		format dt_LMP date9.;
		if '01JUL2016'd <= dt_LMP <= '01FEB2022'd;
	run;

	%*Output counts;
	proc freq data=ana.preg_cohort_&lmpindex noprint;
		table preg_outcome_clean / missing out=ana.preg_outc_final_&lmpindex;
	run;

	%*Finally, delete the macro datasets;
	proc datasets gennum=all;
		delete _: pair_eval:;
	run;


%mend;



*Now run this for two different LMPs at index;
%reassign_lmps(63);
%reassign_lmps(42);
%reassign_lmps(84);








/***************************************************************************************************************************

												05 - DERIVE LMP RANGES PREGNANCIES

Create datasets for those sensitivity analyses where we assume that we have over and under estimated the LMP for 
non-delivery outcomes by 3 weeks.

I decided to do this by following this logic:
(1) Revise the LMPs of all the non-live birth pregnancies in the originally derived cohort.
(2) Run the cleaning algorithm on this cohort.

Determined: This was the most logically consistent way to do it with our original approach. GA at index is important not only
for identifying outcomes but also for identifying the target population.

***************************************************************************************************************************/

/*
MACRO: reassign_lmps_range
PURPOSE: To revise the LMPs for those pregnancies where the LMP was revised based upon the range provided.
*/


%macro reassign_lmps_range(lmpindex=63,revision=,suffix=);

	/**We want to identify those pregnancies with UNK outcomes that were assigned 140d for their gestational length */
	*For these pregnancies, assign their LMP by subtracting 6*7 or 42 from the date of their indexing prenatal encounter;
	data _pregnancies2;
	set temp.pregnancies;

		%*For those pregnancies with unknown outcomes and no gestational age codes, reassign their LMP as occurring 
		at the prespecified gestational age at the index date (&lmpindex);
		if preg_outcome_clean = 'UNK' and lmp_algbased = 0 then dt_LMP = DT_INDEXPRENATAL - &lmpindex;
			else dt_LMP = dt_LMP;

		%*Create indicators for the original values of dt_lmp, indexprenatal and gapreg. These will be overwritten for 
		pregnancies that we clean.;
		dt_lmp_orig = dt_lmp;
		dt_indexprenatal_orig = dt_indexprenatal;
		dt_gapreg_orig = dt_gapreg;

		%*CDL: 2.14.25 -- FIXED SO THAT ONLY MODIFY LMP FOR NON-LIVE BIRTH OUTCOMES;
		%*NEW STEP -- For those pregnancies with non-live birth outcomes, reassign their LMP based on the provided range.;
		if preg_outcome_clean in ('LBS' 'LBM' 'UDL') then dt_lmp = dt_lmp;
			else dt_lmp = dt_lmp &revision;

	run;
		
	/*Some pregnancies with unknown outcomes have too long of gestational lengths, which we define as
	> 43w6d or 307 days.*/
	%revise_too_long(_pregnancies2, _toolong, &lmpindex, days=307, round=0);
	/*Output; _toolong*/

	/*Some pregnancies have their indexing prenatal encounter prior to their LMP or too early in gestation (we define as 28 days). 
	We will revise these pregnancies.*/
	%revise_index_before_lmp(_toolong, _before_index, &lmpindex, days=28, round=0);
	

	%*Some pregnancies overlap with others. That is dealt with below;
	%clean_overlap(INPUT_DATA=_before_index, OUTPUT_DATA=_pregnancies_clean, LMPINDEX=&lmpindex);

	*Output the final pregnancy cohort. Subset to those pregnancies with
	LMPs within the relevant window;
	data ana.preg_cohort_&lmpindex._&suffix;
	set _pregnancies_clean; 
		format dt_LMP date9.;
		if '01JUL2016'd <= dt_LMP <= '01FEB2022'd;
	run;

	%*Output counts;
	proc freq data=ana.preg_cohort_&lmpindex noprint;
		table preg_outcome_clean / missing out=ana.preg_outc_final_&lmpindex;
	run;

	%*Finally, delete the macro datasets;
	proc datasets gennum=all;
		delete _: pair_eval:;
	run;


%mend;


/*%reassign_lmps_range(lmpindex=63, revision=- 21,suffix=m21);*/
/*%reassign_lmps_range(lmpindex=63, revision=+ 21, suffix=p21);*/

