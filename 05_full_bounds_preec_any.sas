/********************************************************************************************************************************************
PROGRAM: 05_full_bounds_preec_any.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to conduct a sensitivity analysis where we map out the full set of bounds for the treatment effect
estimat for preeclampsia (any service location), under different assumptions about pregnancies with unobserved outcomes.
	
Goal: 
Output data: 

Date: 
********************************************************************************************************************************************/






/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - IDENTIFY THE BOUND OUTCOMES
	- 02 - NONE EXPERIENCED THE STUDY OUTCOME
	- 03 - BOTH EXPERIENCED THE STUDY OUTCOME
	- 04 - NIFEDIPINE EXPERIENCED THE STUDY OUTCOME
	- 05 - LABETALOL EXPERIENCED TEH STUDY OUTCOME

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
%setup(sample=full, programname=05_full_bounds_preec_any, savelog=Y);


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

													01 - IDENTIFY THE BOUND OUTCOMES

********************************************************************************************************************************************/


*Now create the outcomes for the bounds;
data full_bounds;
set ana.preeclampsia_outc;
	
	*Assume that none experienced the study outcome;
	if preec_any_outc = 'UNK' then none_outc = 0;
		else if preec_any_outc = 'PRE' then none_outc = 1;
		else none_outc = 0;

	*Assume that both experienced the study outcome;
	if preec_any_outc = 'UNK' then both_outc = 1;
		else if preec_any_outc = 'PRE' then both_outc = 1;
		else both_outc = 0;

	*Assume that nifedipine only experienced the study outcome;
	if preec_any_outc = 'UNK' and trt=1 then nif_outc = 1;
		else if preec_any_outc = 'UNK' and trt = 0 then nif_outc = 0;
		else if preec_any_outc = 'PRE' then nif_outc = 1;
		else nif_outc = 0;

	*Assume that labetalol only experienced the study outcome;
	if preec_any_outc = 'UNK' and trt=1 then lab_outc = 0;
		else if preec_any_outc = 'UNK' and trt = 0 then lab_outc = 1;
		else if preec_any_outc = 'PRE' then lab_outc = 1;
		else lab_outc = 0;

run;

/**Check that everything looks write;*/
/*proc freq data=full_bounds;*/
/*	table none_outc*trt both_outc*trt nif_outc*trt lab_outc*trt / missing;*/
/*run;*/
/**/
/*proc freq data=full_bounds (where=(preec_any_outc = 'UNK'));*/
/*	table none_outc*trt both_outc*trt nif_outc*trt lab_outc*trt / missing;*/
/*run;*/





		

/********************************************************************************************************************************************

												02 - NONE EXPERIENCED THE STUDY OUTCOME

Assume that neither nifedipine nor labetalol initiators without pregnnacy outcomes experienced pregnancy loss.

********************************************************************************************************************************************/

%bound_analyses(
		inds=full_bounds,
		outc=none, 
		studyoutc=pany, 
		bootnum=1000,
		psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
		psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre,
		numga=2,
		trtvar=trt,
		gacatvar=ga_index_cat
		);


/********************************************************************************************************************************************

												03 - BOTH EXPERIENCED THE STUDY OUTCOME

Assume that both nifedipine nor labetalol initiators without pregnnacy outcomes experienced pregnancy loss.

********************************************************************************************************************************************/

%bound_analyses(
		inds=full_bounds,
		outc=both, 
		studyoutc=pany, 
		bootnum=1000,
		psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
		psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre,
		numga=2,
		trtvar=trt,
		gacatvar=ga_index_cat
		);


/********************************************************************************************************************************************

											04 - NIFEDIPINE EXPERIENCED THE STUDY OUTCOME

********************************************************************************************************************************************/

%bound_analyses(
		inds=full_bounds,
		outc=nif, 
		studyoutc=pany, 
		bootnum=1000,
		psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
		psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre,
		numga=2,
		trtvar=trt,
		gacatvar=ga_index_cat
		);




/********************************************************************************************************************************************

											05 - LABETALOL EXPERIENCED THE STUDY OUTCOME

********************************************************************************************************************************************/

%bound_analyses(
		inds=full_bounds,
		outc=lab, 
		studyoutc=pany, 
		bootnum=1000,
		psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre
			depressi_post anxiety_post antideprx_post benzorx_post 
			teratrx_pre num_outptpnc num_outptpnc_2,
		psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre 
			depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre,
		numga=2,
		trtvar=trt,
		gacatvar=ga_index_cat
		);


%count_missing_zero(inds1=ana.full_bounds_both_loss_boot_1, inds2=ana.full_bounds_both_loss_boot_2);
%count_missing_zero(inds1=ana.full_bounds_none_loss_boot_1, inds2=ana.full_bounds_none_loss_boot_2);
%count_missing_zero(inds1=ana.full_bounds_nif_loss_boot_1, inds2=ana.full_bounds_nif_loss_boot_2);
%count_missing_zero(inds1=ana.full_bounds_lab_loss_boot_1, inds2=ana.full_bounds_lab_loss_boot_2);


