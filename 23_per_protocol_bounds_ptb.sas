/********************************************************************************************************************************************
PROGRAM: 23_per_protocol_bounds_ptb.sas
PROGRAMMER: Chase Latour
PURPOSE: To calculate bounds on the per-protocol effect for preterm live birth. 
	
Goal: 
Output data: 

Date: 6.16.2025
********************************************************************************************************************************************/






/********************************************************************************************************************************************

TABLE OF CONTENTS:
	- 00 - SET UP LIBRARIES
	- 01 - IDENTIFY DISCONTINUATION OF THE STUDY MEDICATION
	- 02 - CALCULATE THE BOUNDS FOR PRETERM BIRTH

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
%setup(sample=full, programname=23_per_protocol_bounds_ptb, savelog=Y);

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

											01 - IDENTIFY DISCONTINUATION OF THE STUDY MEDICATION

We assume that gaps in medication fills of >30 days represent
deviation from the assigned protocol. This inherently ignores potential legitimate reasons for discontinuation (e.g., adverse side effects),
but we are not confident that these can be reasonably identified in claims data, thus necessitating this simplifying assumption. Further,
we ignore post-index fills for other antihypertensives, assuming that these instances represent appropriate up-titration of the initial
treatment, in line with the specified protocol.

********************************************************************************************************************************************/



*The file temp.postidx_antihypertensives_63 contains information for all the antihypertensive fills for a pregnancy
between their index date and the outcome date. The 63 (9 weeks) represents the gestational age that we assumed at index for pregnancies
with UNK outcomes and no gestational age information. This should index lookback based on the first fill date (i.e., initiation), not
the LMP.;


*First, merge the treatment variable onto the temp.postidx_antihypertensives_63 dataset;
proc sql;
	create table pregnancies as
	select a.trt, a.dt_gapreg, b.*
	from ana.primary_cohort as a
	left join temp.postidx_antihypertensives_63 as b
	on  a.idxpren=b.idxpren
	;
	quit;
	


*Now, we want to know information about when the person discontinued their initial antihypertensive, if they did.

For this analysis, we define discontinuation as experiencing a 30-day gap in medications at home, based upon days supply
for the medication, allowing for a 30-day gap in fills.;

*Sort by svcdate for the antihypertensive fill. Subset to those rows where the antihypertensive filled
matches the initial antihypertensive;
proc sort data=pregnancies out=cohort ( where = (exposure = atc_label));
	by enrolid idxpren svcdate;
run;

/**Confirm that everyone is still represented - Everyone should at least have thier first indexing fill;*/
/*proc sql;*/
/*	select count(distinct enrolid) as num_people, count(distinct idxpren) as num_preg*/
/*	from cohort;*/
/*	quit;*/
/**Count: 5557;*/

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
				next_fill = svcdate + daysupp + abs(days_between); 
				last_date = next_fill;
			end;
				else do;
					next_fill = svcdate + daysupp;	
					last_date = next_fill;
				end; 
		end;

	drop last_date;

	if first.idxpren then count30 = 1;
		else if days_between > 30 then count30 = count30+1;
		else count30 = count30;

run;


%*Now the END date of the first gap will be the next_fill date for the last row where count = 1;
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

*******Join the gap information and the date of the end of the gap onto the primary analysis dataset;

proc sql;
	create table per_protocol_analysis as
	select a.*, b.next_fill, b.gap
	from ana.primary_cohort as a
	left join exposure2 as b
	on a.enrolid=b.enrolid and a.idxpren=b.idxpren
	;
	quit;














/********************************************************************************************************************************************

											02 - CALCULATE THE BOUNDS FOR PRETERM BIRTH

We defined the date of protocol deviation according to the end of the 30-day grace period. We then considered outcome-specific assumptions 
about pregnancies that deviated from the protocol to construct the full bounds on the per-protocol effect. 

For preterm birth, we assumed: 
1. Nifedipine initiators experienced preterm birth at 23 weeks’ gestation or immediately upon deviation from the protocol if after 
	23 weeks’ gestation, while labetalol initiators experienced a term delivery, and 
2. Vice-versa. 

We construct the bounds first for 1 then for 2.
********************************************************************************************************************************************/



*Create the dataset to work from;

data per_protocol_bounds;
set per_protocol_analysis;

	****Assumption: Nifedipine initiators experienced preterm birth at 23 weeks’ gestation or immediately upon deviation from the protocol if after 
	23 weeks’ gestation, while labetalol initiators experienced a term delivery;

	*Event indicator;
	if trt=1 and gap = 1 then preg_outcome_nif = 'PTB';
		else if trt=0 and gap = 1 then preg_outcome_nif = 'TB';
		else preg_outcome_nif = preg_outcome_ptb;

	*New event date variable;
	if trt = 1 and gap = 1 then dt_gapreg_nif = max(dt_lmp + 161, dt_gapreg_ptb);
		else if trt = 0 and gap = 1 then dt_gapreg_nif = dt_lmp + 258;
		else dt_gapreg_nif = dt_gapreg_ptb;

	****Assumption: Labetalol initiators experienced preterm birth at 23 weeks’ gestation or immediately upon deviation from the protocol if after 
	23 weeks’ gestation, while nifedipine initiators experienced a term delivery;

	*Event indicator;
	if trt = 0 and gap = 1 then preg_outcome_lab = 'PTB';
		else if trt=1 and gap = 1 then preg_outcome_lab = 'TB';
		else preg_outcome_lab = preg_outcome_ptb;

	*New event date variable;
	if trt = 0 and gap = 1 then dt_gapreg_lab = max(dt_lmp + 161, dt_gapreg_ptb);
		else if trt = 1 and gap = 1 then dt_gapreg_lab = dt_lmp + 258;
		else dt_gapreg_lab = dt_gapreg_ptb;

run;



***Create a macro that allows you to calculate each of the bounds;


%macro perprotocol_bounds_ptb(
		allassumptions= nif lab,
		numboot = 2,
		psvars=ga_quartile age_at_index age_at_index_2
					year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
					nausea_pre recurlos_pre obesity_post chronichypertension_pre
					depressi_post anxiety_post antideprx_post benzorx_post 
					teratrx_pre num_outptpnc num_outptpnc_2,
		dovars=ga_quartile age_at_index year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
					nausea_pre recurlos_pre obesity_post chronichypertension_pre depressi_post anxiety_post antideprx_post 
					benzorx_post teratrx_pre num_outptpnc num_outptpnc_2 rural2,
		dovarsmodel=ga_quartile age_at_index year_index4
					t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp 
					nausea_pre recurlos_pre obesity_post chronichypertension_pre 
					depressi_post anxiety_post antideprx_post benzorx_post
					teratrx_pre num_outptpnc num_outptpnc_2 rural2,
		psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
					nausea_pre recurlos_pre obesity_post chronichypertension_pre 
					depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre,
		doclassvars=ga_quartile year_index4 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
					nausea_pre recurlos_pre obesity_post chronichypertension_pre
					depressi_post anxiety_post antideprx_post benzorx_post teratrx_pre rural2);

	%*Count the number of assumptions;
	%let numassumptions = %sysfunc(countw(&allassumptions));

	%*Repeat the analysis over each assumption;

	%do n=1 %to &numassumptions;

		%let assumption = %scan(&allassumptions, &n);

		%*Gather the point estimate;
		%competing3risk_weights(
			boot=0, 
			inds=per_protocol_bounds, 
			gacatvar = ga_index_cat,
			startDT=dt_index, 
			outcomevar = preg_outcome_&assumption,
			eventDT=dt_GApreg_&assumption, 
			event = 'PTB',
			cr1='IAB', 
			cr2= 'SAB' 'UAB' 'SB' 'MLS',
			cr3= 'TB',
			psvars=&psvars,
			dovars=&dovars,
			dovarsmodel = &dovarsmodel,
			psclassvars=&psclassvars, 
			doclassvars=&doclassvars,
			trtvar=trt,
			numiterations=1,
			initialseed=23244, 
			outds=ana.perpro_point_ptb_&assumption,
			outds_dist=gadist_perpro_point_ptb_&assumption,
			outds_ps = NA,
			outds_surv = NA
			);

		%*Combine the stratified estimates according to the distribution of GA at index among
		those in the treated group;
		%combine_point_estimates(
			gadist = ana.gadist_perpro_point_ptb_&assumption,
			numgastrat = 2,
			est = ana.perpro_point_ptb_&assumption._,
			outds = ana.perpro_point_ptb_overall_&assumption
			);

		%*********Finally, conduct the bootstrap;
		%competing3risk_weights(
			boot=1, 
			inds=per_protocol_bounds, 
			gacatvar = ga_index_cat,
			startDT=dt_index, 
			outcomevar = preg_outcome_&assumption,
			eventDT=dt_GApreg_&assumption, 
			event = 'PTB',
			cr1='IAB', 
			cr2='SAB' 'UAB' 'SB' 'MLS',
			cr3='TB',
			psvars=&psvars,
			dovars=&dovars,
			dovarsmodel = &dovarsmodel,
			psclassvars= &psclassvars, 
			doclassvars= &doclassvars,
			trtvar=trt, 
			numiterations=&numboot,
			initialseed=23244, 
			outds=ana.perpro_boot_ptb_&assumption,
			outds_dist=gadist_perpro_boot_ptb_&assumption,
			outds_ps = NA,
			outds_surv = NA
			);
			
			
		options mlogic mprint symbolgen notes;
		*Combine all the bootstrapped estimates;
		%combine_boot_estimates(
					inputEst= ana.perpro_boot_ptb_&assumption, 
					inputDist= ana.gadist_perpro_boot_ptb_&assumption,
					numStrata= 2, 
					output_stratified= ana.perpro_boot_ptb_strat_&assumption,
					output_overall= ana.perpro_boot_ptb_overall_&assumption);


		*Stratified estimates with confidence intervals.;
		%strat_estimates_w_CI(bootdsn=ana.perpro_boot_ptb_&assumption._1, pointdsn=ana.perpro_point_ptb_&assumption._1, output=ana.perpro_boot_ptb_&assumption._1_ci);
		%strat_estimates_w_CI(bootdsn=ana.perpro_boot_ptb_&assumption._2, pointdsn=ana.perpro_point_ptb_&assumption._2, output=ana.perpro_boot_ptb_&assumption._2_ci);

		*OVerall estimates with confidence interval;
		%overall_estimates_w_CI(stderrdsn=ana.perpro_boot_ptb_overall_&assumption, pointdsn=ana.perpro_point_ptb_overall_&assumption, output=ana.perpro_all_boot_ptb_&assumption._ci);

	%end;

%mend;
	
*options mprint;
%perprotocol_bounds_ptb(numboot=1000);

/****Assumption: Nifedipine had an outcome;*/
%count_missing_zero(inds1=ana.perpro_boot_ptb_nif_1, inds2=ana.perpro_boot_ptb_nif_2);

/****Assumption: Labetalol had an outcome;*/
%count_missing_zero(inds1=ana.perpro_boot_ptb_lab_1, inds2=ana.perpro_boot_ptb_lab_2);




