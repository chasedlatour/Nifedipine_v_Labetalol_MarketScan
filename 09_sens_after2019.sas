/********************************************************************************************************************************************
PROGRAM: 09_sens_after2019.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to run the sensitivity analysis where we restrict to individuals with index dates after January 1, 2019.
	
Goal: 
Output data: 

Date: 12.19.2024
********************************************************************************************************************************************/


	



/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - RUN ANALYSIS

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
%setup(sample=full, programname=09_sens_after2019, savelog=Y);

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

														01 - RUN ANALYSIS

********************************************************************************************************************************************/

%repeat_in_subset(
	inds=ana.preg_covar_63_dt_index_270,
	sens=p2019,
	where= dt_index >= 21915, /*January 1, 2020*/
	numboot=1000
);




%count_missing_zero(inds1=ana.p2019_boot_1, inds2=ana.p2019_boot_2);
%count_missing_zero(inds1=ana.p2019_boot_ptb_1, inds2=ana.p2019_boot_ptb_2);
