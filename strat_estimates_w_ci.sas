/*
MACRO: strat_estimates_w_ci
PURPOSE: To calculate the final stratified estimates with the confidence interval.

INPUTS:
- bootdsn -- The dataset with the estimates from the bootstrap for the strata.
- pointdsn -- The dataset with the point estimates for the strata
- output -- The name of the output dataset
*/


%macro strat_estimates_w_ci(bootdsn=, pointdsn=, output=);

	%*Calculate the standard errors for the confidence intervals;;
	proc means data=&bootdsn noprint;
		var e0_risk e1_risk riskDiff lnriskRatio;
		OUTPUT OUT=_strata_se
	        STD=se_e0_risk se_e1_risk se_riskDiff se_riskRatio;
	run;

	%*Merge together the point estimates with the standard errors from the bootstrapped estimates
	and output the final dataset;
	data &output;
	merge &pointdsn _strata_se;

		z = probit(0.975);

		risk0 = round(100*e0_risk,0.1);
		risk0_lcl = round(100*(e0_risk - (z*se_e0_risk)),0.1);
		risk0_ucl = round(100*(e0_risk + (z*se_e0_risk)),0.1);
		risk0_est = cat(risk0, " (", risk0_lcl, ", ", risk0_ucl, ")");

		risk1 = round(100*e1_risk,0.1);
		risk1_lcl = round(100*(e1_risk - (z*se_e1_risk)),0.1);
		risk1_ucl = round(100*(e1_risk + (z*se_e1_risk)),0.1);
		risk1_est = cat(risk1, " (", risk1_lcl, ", ", risk1_ucl, ")");

		rd = round(100*riskDiff,0.1);
		rd_lcl = round(100*(riskDiff - (z*se_riskDiff)),0.1);
		rd_ucl = round(100*(riskDiff + (z*se_riskDiff)),0.1);
		rd_est = cat(rd, " (", rd_lcl, ", ", rd_ucl, ")");
			
		rr = round(riskRatio,0.01);
		rr_lcl = round(exp(lnriskRatio - (z*se_riskRatio)), 0.01);
		rr_ucl = round(exp(lnriskRatio + (z*se_riskRatio)), 0.01);
		rr_est = cat(rr, " (", rr_lcl, ", ", rr_ucl, ")");

		keep risk0_est risk1_est rd_est rr_est;
	run;	

	proc datasets gennum=all noprint;
		delete _strata_se;
	run;

%mend;

