/********************************************************************************************************************************************
PROGRAM: 01a_import_preg_data.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to import the relevant pregnancy information into this project directory from the raw_data
directory. This should NOT be re-run once completed.
	

Date: 1.15.2025
********************************************************************************************************************************************/











/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - IMPORT NECESSARY DATASETS

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
%setup(sample=full, programname=01a_import_preg_data, savelog=Y);

*Point to the temp files in the pregnancy derivation folder.;
/*libname ogtemp "/local/projects/marketscan_preg/raw_data/data/random1pct";*/
libname ogtemp "/local/projects/marketscan_preg/raw_data/data/full";
/*libname logtemp slibref=ogtemp server=server;*/

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
%inc "/local/projects/marketscan_preg/Latour_23_2322/programs/FormatStatements_CDWH.sas";









/********************************************************************************************************************************************

												01 - IMPORT NECESSARY DATASETS

********************************************************************************************************************************************/

*We are only interested in those pregnancies from algorithm 30-3 with LMPs that occurred on or after January 1, 2015;
proc sql;
	create table temp.pregnancies as
	select distinct *
	from preg.pregnancy_lmp_simp_all
	where algorithm = "30-3" and year(dt_lmp) >= 2015;
	quit;

*Get all the gestational age encounters and save in the temporary folder to ensure no revisions in subsequent runs
	of the pregnancy cohort. Subset to only those encounters from people in our pregnancy cohort.;
proc sql;
	create table temp.gestage as
	select *
	from ogtemp.gestage
	where patient_deid in (select distinct patient_deid from temp.pregnancies)
	;
	quit;


*Output the getational age encounters;
proc sql;
	create table temp.gestagepren as
	select distinct patient_deid, enc_date, parent_code, code, max(prenatal_enc) as Pren_GA_enc,
             preg_outcome, code_hierarchy, gestational_age_days, gest_age_wks, 
				case when min_gest_age = 14 then 28 else min_gest_age end as min_gest_age, max_gest_age,
             zhu_test, zhu_hierarchy
    from (select * from temp.gestage where code_hierarchy not in ("Prenatal care","Prenatal Care","Weight only - preterm"))
    group by patient_deid, parent_code, enc_date, code, preg_outcome, code_hierarchy, 
               gestational_age_days, min_gest_age, max_gest_age 
	;
    quit;

*Output the prenatal care encounter dates. We only need those encounter dates that occurred on or after 2015;
data temp.codeprenatal_meg1_dts;
set ogtemp.codeprenatal_meg1_dts;
	where year(enc_date) >= 2015;
run;

