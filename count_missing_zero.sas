/*
MACRO: count_missing_zero
PURPOSE: To count the number of GA-stratified bootstraps where at least 1 strata had a zero risk of missing effect estimate
*/
*Count the number at risk or missing an estimate;


%macro count_missing_zero(inds1=, inds2=);

	%**First stack the bootstrapped estimates across strata;
	proc sql;
		create table stacked as
		select 1 as strata, *
		from &inds1
		%IF &inds2 ne NA %THEN %DO;
			union corr
			select 2 as strata, *
			from &inds2
			%END;
		;
		quit;

	%*Output counts;
	proc sql;
		select count(*) as n_row, sum(e0_risk = 0) as risk0_0, sum(e1_risk = 0) as risk1_0,
			sum(e0_risk = .) as missing_risk0, sum(e1_risk = .) as missing_risk1,
			sum(riskDiff = .) as missing_rd, sum(lnriskRatio = .) as missing_lnrr
		from stacked
		;
		quit;

	%*Output counts;
	proc sql;
		select count(*) as n_row, sum(e0_risk = 0) as risk0_0, sum(e1_risk = 0) as risk1_0,
			sum(e0_risk = .) as missing_risk0, sum(e1_risk = .) as missing_risk1,
			sum(riskDiff = .) as missing_rd, sum(lnriskRatio = .) as missing_lnrr
		from stacked
		group by strata
		;
		quit;

%mend;

