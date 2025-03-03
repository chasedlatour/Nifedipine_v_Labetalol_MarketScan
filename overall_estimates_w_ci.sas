/*
MACRO: overall_estimates_w_ci
PURPOSE: To calculate the final stratified estimates with the confidence interval.

INPUTS:
- stderrdsn -- The dataset with the standard error estimates.
- pointdsn -- The dataset with the point estimates for the strata
- output -- The name of the output dataset
*/


%macro overall_estimates_w_ci(stderrdsn=, pointdsn=, output=);

	%*Transpose the dataset with distributional statistics for each measure;
	proc transpose data=&stderrdsn out=_se prefix=se_;
		id estimate;
		var stderr;
	run;

	%*Merge together the point estimates with the standard errors from the bootstrapped estimates
	and output the final dataset;
	data &output;
	merge &pointdsn _se;

		z = probit(0.975);

		risk0_lcl = round(100*(risk0 - (z*se_risk0)),0.1);
		risk0_ucl = round(100*(risk0 + (z*se_risk0)),0.1);
		risk0 = round(100*risk0,0.1);
		risk0_est = cat(risk0, " (", risk0_lcl, ", ", risk0_ucl, ")");

		risk1_lcl = round(100*(risk1 - (z*se_risk1)),0.1);
		risk1_ucl = round(100*(risk1 + (z*se_risk1)),0.1);
		risk1 = round(100*risk1,0.1);
		risk1_est = cat(risk1, " (", risk1_lcl, ", ", risk1_ucl, ")");

		rd_lcl = round(100*(rd - (z*se_rd)),0.1);
		rd_ucl = round(100*(rd + (z*se_rd)),0.1);
		rd = round(100*rd,0.1);
		rd_est = cat(rd, " (", rd_lcl, ", ", rd_ucl, ")");
			
		rr_lcl = round(exp(log(rr) - (z*se_rr)), 0.01);
		rr_ucl = round(exp(log(rr) + (z*se_rr)), 0.01);
		rr = round(rr,0.01);
		rr_est = cat(rr, " (", rr_lcl, ", ", rr_ucl, ")");

		keep risk0_est risk1_est rd_est rr_est;
	run;	

	proc datasets gennum=all noprint;
		delete _se;
	run;

%mend;

