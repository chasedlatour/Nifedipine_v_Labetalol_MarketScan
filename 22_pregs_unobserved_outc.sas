/********************************************************************************************************************************************
PROGRAM: 22_pregs_unobserved_outc.sas
PROGRAMMER: Chase Latour
PURPOSE: To provide descriptive statistics regarding pregnancies with unobserved outcomes. These analyses focus on those pregnancies
with unobserved outcomes that are in the primary cohort.
	
Goal: 
Output data: 

Date: 5.26.2025
********************************************************************************************************************************************/






/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - DESCRIBE GESTATIONAL AGE CODES

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
%setup(sample=full, programname=22_pregs_unobserved_outc, savelog=Y);

******************************************************************************************************************************************;
/*Create local mirros of the libraries from the set up macro - Run locally*/
/*libname lout slibref=out server=server; */
/*libname lana slibref=ana server=server;*/
/*libname lwork slibref=work server=server;*/
/*libname ltemp slibref=temp server=server;*/
/*libname lpreg slibref=preg server=server;*/

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

													01 - DESCRIBE GESTATIONAL AGE CODES

We are interested in understanding first what types of gestational age codes overlap with the pregnancies with unobserved outcomes.
In particular, we are interested in codes from claims between the indexing prenatal claim and the outcome. All of the pregnancies
in the primary cohort should have a non-missing indexing prenatal claim.

********************************************************************************************************************************************/


%*First, create a local copy of the primary cohort;
data primary;
set ana.primary_cohort;
run;


%*Link on the gestational-age claims;
proc sql;
	create table ga_codes as
	select a.*, b.*
	from primary as a
	left join temp.gestage as b
	on a.enrolid=b.enrolid and a.dt_indexprenatal <= b.enc_date <= a.dt_gapreg + 7
	;
	quit;


%*Provide a descriptive of the code_hierarchy variable;
proc freq data=ga_codes;
	table code_hierarchy / missing;
run;


%*Create a hierarchical assignment according to the pregnancy outcome;
data ga_codes2;
set ga_codes;

	if code_hierarchy = 'Weight only - preterm' then code_hierarchy = '';
		else code_hierarchy = code_hierarchy;

	if preg_outcome_clean in ('SB' 'UDL' 'LBS' 'LBM' 'MLS') then do;
		if code_hierarchy = 'Specific gestational age' then code_hier_num = 1; 
			else if code_hierarchy = 'Extreme prematurity' then code_hier_num = 2; 
			else if code_hierarchy in ('Other preterm' 'Post-term') then code_hier_num = 3;
			else code_hier_num = 9; 
	end;
		else if preg_outcome_clean in ('SAB' 'UAB' 'IAB') then do;
			if code_hierarchy = 'Specific gestational age' then code_hier_num = 1;
				else if code_hierarchy = 'Abortion' then code_hier_num = 2; 
				else code_hier_num = 9; 
		end;
		else if preg_outcome_clean = 'EM' then do;
			if code_hierarchy = 'Specific gestational age' then code_hier_num = 1; 
			else if code_hierarchy = 'Extreme prematurity' then code_hier_num = 2; 
			else if code_hierarchy = 'Other preterm' then code_hier_num = 3;
			else code_hier_num = 9; 
		end;
		else if preg_outcome_clean = 'UNK' then do;
			if code_hierarchy = 'Specific gestational age' then code_hier_num = 1; 
			else if code_hierarchy = 'Extreme prematurity' then code_hier_num = 2; 
			else if code_hierarchy in ('Other preterm' 'Post-term') then code_hier_num = 3;
			else if code_hierarchy = 'Missing' then code_hier_num = 4;
			else if code_hierarchy = 'Prenatal Care' then code_hier_num = 5;
			else code_hier_num = 9; 
		end;

run;

proc freq data= ga_codes2;
	table code_hierarchy * code_hier_num / missing;
run;

%*Now output the values of the variables with the minimum hierarchy;
proc sql;
	create table ga_codes3 as
	select distinct enrolid, idxpren, code_hier_num, preg_outcome_clean
	from ga_codes2
	group by enrolid, idxpren
	having code_hier_num = min(code_hier_num)
	;
	quit;

%*Now do a proc-freq output;
proc freq data=ga_codes3;
	table preg_outcome_clean*code_hier_num / missing;
run;


%*Finding: Almost 75% of pregnancies with unobserved outcomes were assigned a GA via a specific GA code. The rest
were majority assigned becuase of a Missing code (e.g., 1st, second trimester);








%********************************************
	Look more deeply at those encounters for pregnancies with unobserved outcomes that did not have a specific GA code;

%*Subset to those pregnancies;
data unobserved;
set ga_codes3;
	where preg_outcome_clean = 'UNK' and code_hier_num in (4,9);
run;

%*182 pregnancies;

%*Look at all the prenatal encounters for these pregnancies;
proc sql;
	create table unobserved_pnc as
	select a.*,b.idxpren, b.preg_outcome_clean, b.code_hier_num, b.dt_indexprenatal, b.dt_gapreg
	from preg.codeprenatal as a
	right join (select aa.*, bb.dt_indexprenatal, bb.dt_gapreg
				from unobserved as aa
				left join primary as bb
				on aa.idxpren=bb.idxpren and aa.enrolid=bb.enrolid) as b
	on a.enrolid=b.enrolid and b.dt_indexprenatal <= a.enc_date <= b.dt_gapreg
	;
	quit;

proc sql;
	select count(distinct idxpren) as num_preg from unobserved_pnc
	;
	quit;

%*Sort the data so that it can be easily viewed;
proc sort data=unobserved_pnc;
	by enrolid idxpren enc_date;
run;

%*Select the variables in an easily viewable order;
proc sql;
	create table unobserved_pnc2 as
	select enrolid, idxpren, enc_date, code, code_type, description, dt_indexprenatal, dt_gapreg
	from unobserved_pnc
	;
	quit;

%*Look at the number of prenatal codes and distinct claim dates;
proc sql;
	create table unobserved_pnc3 as
	select distinct enrolid, idxpren, enc_date, code, code_type, description,
		dt_indexprenatal, dt_gapreg, count(distinct enc_date) as num_pnc_dates, count(idxpren) as num_distinct_rows
	from unobserved_pnc2
	group by enrolid, idxpren
	;
	quit;

%*Look at the distribution;
proc means data=unobserved_pnc3 min p10 p25 p50 p75 p90 max nmiss;
	var num_pnc_dates num_distinct_rows;
run;

%*Look at those with less than 7 distinct PNC dates;
data unobserved_pnc_lt7;
set unobserved_pnc3;
	where num_pnc_dates < 7;
run;
