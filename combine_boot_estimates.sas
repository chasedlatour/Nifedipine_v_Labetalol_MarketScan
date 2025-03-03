/*
MACRO: combine_boot_estimates
PURPOSE: The overall purpose of this macro is to output the stratified and pooled estimates for an analysis.

The steps of this macro are as follows:
(1) Create a dataset with all the bootstrapped estimates from an analysis, stratified by gestational age at index and joined
	with infomration on the distribution of pregnancies across gestational age strata at baseline.
(2) Calculate stratified estimates based upon the combined data.
(3) Calculate overall estimates based upon the combined data, after applying standardization.

INPUTS:
- inputEst -- This is the starting part of the names of datasets with the stratified estimates, excluding the _# indicating the strata.
- inputDist -- Name of the dataset with the distribution of individuals across the gestational age categories at index 
- numStrata -- The number of gestational age strata
- output_combined -- Name of the output dataset with the stacked estimates
- output_stratified -- Name of the output dataset with the stratified estimates
- output_overall -- Name of the output dataset with the overall estimates after standardization
*/


%macro combine_boot_estimates(
			inputEst= , 
			inputDist= ,
			numStrata= , 
			output_stratified= ,
			output_overall= );
			
	/*
	Testing:
	%let inputEst= ana.primary_boot;
	%let inputDist= ana.ga_dist_primary_boot;
	%let numStrata= 2;
	%let output_stratified= ana.primary_boot_strat;
	%let output_overall= ana.primary_boot_overall;
	*/
	
	%************
		STEP 0: Get information on the gestational age distribution in each round.;

	proc sql;
		create table _strata as
		select *, sum(count) as total, count / sum(count) as prop
		from &inputDist
		where trt = 1
		group by round;
		quit;
	

	%************
		STEP 1: Create a dataset with all the bootstrapped estimates from an analysis, stratified by gestational age at index 
		and joined with information on the distribution of pregnancies across gestational age strata at baseline.;
	%*Create the combined datasets with the stratified risks and associated proportions;	
	%DO i=1 %TO &numStrata;

		/*
		%let i=1;
		*/

		proc sql;
			create table _combined_&i as 
			select a.round, &i as strata, a.e0_risk as risk0, a.e1_risk as risk1, b.prop
			from &inputEst._&i as a
			left join _strata (where = (ga_strat = &i)) as b
			on a.round=b.round;
			quit;

		%*Stack all of the datasets together;
		%if &i=1 %then %do;
			data _output_combined;
			set _combined_&i;
			run;
		%end;
		%else %do;
			proc append base=_output_combined data=_combined_&i;
			run;
		%end;

	%END;
	
	%*Determine if any rounds have missing values due to non-convergence, etc. ;
	data _combined_miss;
	set _output_combined;
		if risk0 = . then risk0_missing = 1;
			else risk0_missing = 0;
		if risk1 = . then risk1_missing = 1;
			else risk1_missing = 0;
		risk_missing = sum(risk0_missing, risk1_missing);
	run;
	
	%*Remove those rounds that returned missing risk values.;
	proc sql;
		create table _combined_nomiss as
		select distinct *, sum(risk_missing) as missing_risks
		from _combined_miss
		group by round;
		quit;
	data _combined_nomiss; set _combined_nomiss; where missing_risks = 0; run;
	%*Output counts for evaluating log;
	proc sql noprint;
		select count(*) into :num_all_combined from _combined_miss;
		select count(*) into :num_post_missing from _combined_nomiss;
		quit;
	%PUT Number of rounds due to missing risk information: %eval(&num_all_combined - &num_post_missing);
			




	%**********
		STEP 2: Calculate stratified estimates based upon the combined data.;

	%*Calculate a risk difference and risk ratio;
	data _combined_strat;
	set _combined_nomiss;
		rd = risk1 - risk0;
		rr = risk1 / risk0;
	run;
	proc sort data= _combined_strat;
		by strata round;
	run;

	%let estimates = risk0 risk1 rd rr;

	%*Output distributional statistics of each estimate as well as the standard error;

	%DO i=1 %TO 4;

		%let loop&i = %scan(&estimates, &i);

		proc univariate data=_combined_strat noprint;
			by strata;
			var &&loop&i;
			output out=&&loop&i std=stderr pctlpts=2.5 50 97.5 pctlpre=p;
		run;

	%END;

	%*Stack all of the datasets together;
	proc sql;
		create table &output_stratified as
		select "risk0" as estimate, * from risk0
		union corr
		select "risk1" as estimate, * from risk1
		union corr
		select "rd" as estimate, * from rd
		union corr
		select "rr" as estimate, * from rr
		;
		quit;

	%*Delete the datasets;
	proc datasets gennum=all noprint;
		delete risk0 risk1 rd rr;
	run;


	%***********
		STEP 3: Calculate overall estimates based upon the combined data, after applying standardization.;

	%*Now calculate the summarized risks;
	proc sort data=_output_combined out=primary_combined;
		by round strata;
	run;
	%*Calculate the standardized risks;
	data primary_combined;
	set primary_combined;
		risk0_std = risk0*prop;
		risk1_std = risk1*prop;
	run;
	proc sql;
		create table primary_combined2 as
		select distinct sum(risk0_std) as risk0, sum(risk1_std) as risk1,
				calculated risk1 - calculated risk0 as rd,
				calculated risk1 / calculated risk0 as rr
		from primary_combined
		group by round
		;
		quit;
	


	%*Output distributional stastistics for each;

	%DO i=1 %TO 4;

		%let loop&i = %scan(&estimates, &i);

		proc univariate data=primary_combined2 noprint;
			var &&loop&i;
			output out=&&loop&i std=stderr pctlpts=2.5 50 97.5 pctlpre=p;
		run;

	%END;

	%*Stack all of the datasets together;
	proc sql;
		create table &output_overall as
		select "risk0" as estimate, * from risk0
		union corr
		select "risk1" as estimate, * from risk1
		union corr
		select "rd" as estimate, * from rd
		union corr
		select "rr" as estimate, * from rr
		;
		quit;

	%*Delete the datasets;
	proc datasets gennum=all noprint;
		delete risk0 risk1 rd rr;
	run;

%mend;


