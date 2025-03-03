/*
MACRO: bound_analyses
PURPOSE: To implement the sensitivity analyses where we estimate the full bounds on the treatment effect
estimate.

INPUTS:
- outc_
- INDS = input dataset name with a variable- outcome var - that contains the values both_outc, none_outc, nif_outc, lab_outc
- outc = prefix of the variable where the outcome is stored (both, none, nif, or lab)
- studyoutc = indicator of the study outcome, used for saving datasets (loss, ptb)
- PSVARS = list of variables to be included in the propensity score model (separated by spaces)
- PSCLASSVARS = list of variables to be included in the class statement for the PS model
- NUMGA = number of gestational age strata
- TRTVAR = treatment variable name  (values 0,1)
- bootnum = number of iterations (bootstraps)
- GACATVAR = variable with the gestational age strata
*/


%macro bound_analyses(
		inds=,
		outc= , 
		studyoutc=loss, 
		bootnum=,
		psvars=,
		psclassvars=,
		numga=2,
		trtvar=trt,
		gacatvar=ga_index_cat
		);

	*First, run the macro for the point estimate;
	%full_fup_weights(boot=0,
			inds=&inds, 
			gacatvar = &gacatvar,
			outcomevar = &outc._outc,
			psvars=&psvars, 
			psclassvars=&psclassvars, 
			trtvar=&trtvar, 
	      	numiterations=1, 
			initialseed=23244, 
			outds=ana.full_bounds_&outc._&studyoutc._point,
			outds_dist=ga_fullbounds_&outc._&studyoutc._point
		   );

	*Second, run the macro for the bootstrap for stderr;
	%full_fup_weights(boot=1,
			inds=&inds, 
			gacatvar = &gacatvar,
			outcomevar = &outc._outc,
			psvars=&psvars, 
			psclassvars=&psclassvars, 
			trtvar=&trtvar, 
	      	numiterations=&bootnum, 
			initialseed=23244, 
			outds=ana.full_bounds_&outc._&studyoutc._boot,
			outds_dist=ga_fullbounds_&outc._&studyoutc._boot
		   );

	options mlogic mprint symbolgen notes;

	*Combine the point estimates;
	%combine_point_estimates(
		gadist = ana.ga_fullbounds_&outc._&studyoutc._point,
		numgastrat = 2,
		est = ana.full_bounds_&outc._&studyoutc._point_,
		outds = ana.fullbounds_&outc._&studyoutc._point_all
		);

	*Combine all the bootstrapped estimates;
	%combine_boot_estimates(
				inputEst= ana.full_bounds_&outc._&studyoutc._boot, 
				inputDist= ana.ga_fullbounds_&outc._&studyoutc._boot,
				numStrata= 2, 
				output_stratified= ana.fullbounds_&outc._&studyoutc._boot_strat,
				output_overall= ana.fullbounds_&outc._&studyoutc._boot_all);


	*Stratified estimates with confidence intervals.;
	%do i=1 %to &numga;
		%strat_estimates_w_CI(bootdsn=ana.full_bounds_&outc._&studyoutc._boot_&i, pointdsn=ana.full_bounds_&outc._&studyoutc._point_&i, output=ana.fullbounds_&outc._&studyoutc._&i._ci);
	%end;

	*OVerall estimates with confidence interval;
	%overall_estimates_w_CI(stderrdsn=ana.fullbounds_&outc._&studyoutc._boot_all, pointdsn=ana.fullbounds_&outc._&studyoutc._point_all, output=ana.fullbounds_&outc._&studyoutc._all_ci);

%mend;
