/********************************************************************************************************************************************
PROGRAM: 15_sens_censor_disenroll.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to run the sensitivity analysis where people are censored at their last prenatal encounter prior
to their first disenrollment after the index date, regardless of their later outcome.
	
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
%setup(sample=full, programname=15_sens_censor_disenroll, savelog=Y);

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
	/*CDL: ADDED region 1.5.2026*/
	value $region
		"1" = "Northeast"
		"2" = "North Central"
		"3" = "South"
		"4" = "West"
		"5" = "Unknown";
run;



/********************************************************************************************************************************************

														01 - RUN ANALYSIS

********************************************************************************************************************************************/

*Revise the outcome variable that we are going to use;

data disenroll;
set ana.preg_covar_63_dt_index_270;

	*Determine if the pregnancy outcome occurs after the first lapse in continuous enrollment;
	disenroll = cont_enrl_end_any < dt_gapreg;

	*Revise the pregnancy outcome date and the pregnancy outcomes for those pregnancies that disenroll prior to their outcome;
	if disenroll = 1 then do;
		preg_outcome_clean = 'UNK';
		dt_gapreg = cont_enrl_end_any; *CDL: CHANGED. Previously cont_enrl_end_pnc_any;
	end;
	else do;
		preg_outcome_clean = preg_outcome_clean;
		dt_gapreg = dt_gapreg;
	end;

run;

/**If wanted to get a count on the number of pregnancies affected in the primary dataset;*/
/*data count;*/
/*set ana.primary_cohort;*/
/*	disenroll = cont_enrl_end_any < dt_gapreg;*/
/*run;*/
/*proc freq data=count;*/
/*	table preg_outcome_clean*disenroll;*/
/*run;*/


*Now conduct the analysis;
%repeat_in_subset(
	inds=disenroll,
	sens=denrl,
	where=NA,
	numboot=1000
);


%count_missing_zero(inds1=ana.denrl_boot_1, inds2=ana.denrl_boot_2);
%count_missing_zero(inds1=ana.denrl_boot_ptb_1, inds2=ana.denrl_boot_ptb_2);
