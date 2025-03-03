/*
MACRO: combine_point_estimates
PURPOSE: To combine risk estimates across gestational age strata at baseline.
PROGRAMMER: Chase Latour

INPUTS:
- gadist = Name of the gestational age distribution dataset
- numgastrat = Number of gesetational age strata
- est = Name of the dataset with the stratified treatment effect estimates, excluding GA strat number.
- outds = Name of the output dataset
*/

%macro combine_point_estimates(
		gadist = ana.ga_dist_primary_point,
		numgastrat = 2,
		est = primary_point_,
		outds = ana.primary_point_overall);
		
	%*Combine the stratified estimates according to the distribution of GA at index among
	those in the treated group;
	data _loss_strata;
	set &gadist;
		*Counts of people retained within the treated group;
		where trt = 1;
	run;
	*Get the total number of individuals;
	proc sql noprint;
		select sum(count)
		into :total
		from _loss_strata;
		quit;
	data _loss_strata;
	set _loss_strata;
		*Calculate distribution according to the total treated individuals across the two GA strata;
		prop = count / &total;
	run;
	
	%*Combine the two stratified point estimate datasets then left join on the proportion variable
	that will be used to standardize the risks.;
	proc sql;
		create table _stratified_w_prop as
		select a.*, b.prop /*Join the proportion variable onto the datasets*/
		from (
				/*Union the GA stratified datasets with point estimates together*/
				%do i=1 %to &numgastrat;
				
					%if &i=1 %then %do;
						select *, &i as ga_strat
						from &est.&i
					%end;
					%else %do;
						union corr
						select *, &i as ga_strat
						from &est.&i
					%end;

				%end;
				
				) as a 
		left join _loss_strata as b
		on a.ga_strat = b.ga_strat
		;
		quit;
		
	%*Calculate the standardized risks by multiplying the stratified risk with the proportional distribution;
	data _stratified_w_prop;
	set _stratified_w_prop;
		e0_riskstd = e0_risk * prop;
		e1_riskstd = e1_risk * prop;
	run;
	
	%*Calculate the overall standardized estimates
		For the risks, sum the standardized risks over the GA strata.
		For the RD and RR, contrast the standardized risks accordingly.;
	proc sql;
		create table &outds as
		select sum(e0_riskstd) as risk0, sum(e1_riskstd) as risk1,
			calculated risk1 - calculated risk0 as rd,
			calculated risk1 / calculated risk0 as rr
		from _stratified_w_prop
		;
		quit;
		
	%*Delete unnecessary datasets to keep the working directory clear.;
	proc datasets gennum=all noprint;
		delete _loss_strata _stratified_w_prop;
	run;

%mend;