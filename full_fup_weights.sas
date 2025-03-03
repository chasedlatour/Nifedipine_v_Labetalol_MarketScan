/*
MACRO: full_fup_weights
PURPOSE: To implement an analysis where we calculate risks under a scenario with full follow-up, incorporating SMR weights.

INPUTS:
- BOOT = indicator (1=yes, 0=no) as to whether we want to conduct a bootstrap or not.
- INDS = input dataset name
- OUTCOMEVAR = variable where the outcome is stored
- CENSORDT = variable name for date of censoring event
- PSVARS = list of variables to be included in the propensity score model (separated by spaces)
- PSCLASSVARS = list of variables to be included in the class statement for the PS model
- TRTVAR = treatment variable name  (values 0,1)
- NUMITERATIONS = number of iterations (bootstraps)
- INITIALSEED = optional parameter for initial seed (so results will be reproducible
- OUTDS = name of output dataset, containing one record per iteration. This will be stratified by GA at index
- GACATVAR = The variable that stratifies the pregnancies according to their gestational age at index
- GACAT = Values of GACATVAR that define the gestational age categories at index
- OUTDS_DIST = Name of the output dataset for the distribution of GA at index
*/



/*TESTING:
%let boot=0;
%let inds=full_bounds;
%let gacatvar = ga_index_cat;
%let startDT=dt_index;
%let outcomevar = none_outc;
%let psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			year_index2017*t2dmrx_post year_index2017*t1t2dmrx_post year_index2017*metforrx_post
			diabetes_simp*t2dmrx_post diabetes_simp*t1t2dmrx_post diabetes_simp*metforrx_post
			nausea_pre recurlos_pre obesity_post obesity_post*year_index2017
			chronichypertension_pre substance_use_pre recurlos_pre ckd_post
			thyroid_disorder_post thyroidrx_post thyroid_disorder_post*thyroidrx_post
			depressi_post anxiety_post antideprx_post benzorx_post 
			adhdrx_post teratrx_pre num_outptpnc num_outptpnc_2;
%let psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre substance_use_pre
			recurlos_pre ckd_post thyroid_disorder_post thyroidrx_post
			depressi_post anxiety_post antideprx_post benzorx_post adhdrx_post
			teratrx_pre ;
%let trtvar=trt ;
%let numiterations=1;
%let initialseed=23244;
%let outds=ana.primary_point_noipcw;
%let outds_dist=ga_dist_primary_point_noipcw;
*/

%macro full_fup_weights(
		boot=0,
		inds=pregnancies, 
		gacatvar = ga_index_cat,
		outcomevar = preg_outcome_clean,
		psvars=age_at_index chronichypertension_pre, 
		psclassvars=chronichypertension_pre, 
		trtvar=trt, 
      	numiterations=1, 
		initialseed=23244, 
		outds=ana.comprisk_primary,
		outds_dist=ga_dist_primary
	   );


	%PUT STEP1: If bootstrapping, then resample with replacement.;
	%if &boot = 1 %then %do;
		%*Implement the bootstrap of the samples prior to the loop. This ensures that one stream of random samples are used.;
		%PUT STEP 2: CREATE BOOTSTRAP SAMPLES;
	    proc surveyselect noprint data=&inds out=_anacohorta(rename=(replicate=b) drop=numberhits) 
	    	seed=&initialseed method=urs samprate=1 outhits rep=&numIterations; 
	%end;
	%else %do;
		data _anacohorta;
		set &inds;
		run;
	%end;

	

	%IF &boot = 1 %THEN %DO;
		%PUT ANALYZE THE BOOTSTRAPPED SAMPLES;
	%END;
	%ELSE %DO;
		%PUT ANALYZE THE SAMPLE DATA.;
	%END;

	%*Iterate the analysis over each of the bootstrapped samples. For the primary analysis, the numiterations 
	value should be 1, as there is only 1 sample to analyze;
   	%DO i=1 %TO &numIterations;

		/*Testing
		%let i=1;
		*/ 

		%PUT Bootstrap: &i ;
		%IF &i>1 %then options nomlogic nomprint nosymbolgen nonotes;;

		%*Create an identifier - _id - that is used to uniquely reference each row in the sample. This is important
		for bootstrapping since one person may appear more than once. If bootstrapping, this will also
		subset to the relevant bootstrap sample.;
		%if &boot = 1 %then %do;
			%*Subset to the replicate that we are interested in so that we can conduct analyses within this sample;
			data _anacohorta&i;
			set _anacohorta;
				where b = &i;
				_id = _N_; %*Create an ID variable for each row;
			run;
		%end;
		%else %do;
			data _anacohorta&i;
			set _anacohorta; %*Do not need to subset to boostrapped sample;
				_id = _N_; %*Create an ID variable for each row;
			run;
		%end;


		
		%PUT RUN ANALYSIS FOR EACH STRATA OF GACAT;
		
		%*Gather a vector of the gestational age categories. These must be numerical values.;
		proc sql noprint;
			select distinct &gacatvar
			into :gacat separated by " "
			from _anacohorta&i
			;
			quit;

		%put &gacat;
		%let numGA = %sysfunc(countw(&gacat, %str( )));
		%put Number of gestational age categories in round &i : &numGA;

		%*Now, we will run the analysis within each strata of gestational age for the selected sample.;
		%DO g = 1 %TO &numGA;

			/*Testing:
			%let g=1;
			*/ 
			
			%*Subset to rows within the gestational age strata;
			data _anacohort&g;
			set _anacohorta&i;
				where &gacatvar = %scan(&gacat, &g);
			run;

			%PUT STEP 3: CALCULATE STABILIZED INVERSE PROBABILITY OF TREATMENT WEIGHTS;

		  	%PUT STEP 3a: CALCULATE FIRST PS MODEL PRIOR TO TRIMMING;
	 	  	proc logistic data=_anacohort&g noprint;
				class &psclassVars;
				model &trtvar (reference = '0' ) = &psvars;
				output out=_ps p=ps;
			run;

		  	%PUT STEP 3b: IDENTIFY NON-OVERLAP AND TRIM POPULATION;
		  	*Calculate min & max PSs;
			proc sort data=_ps;
				by &trtvar;
			run;
			proc means data=_ps min max noprint;
		  		by &trtvar;
				var ps;
				output out=_sumstat min=minPS max=maxPS;
		  	run;
		  	*Apply min & max to opposite treatments;
		  	data _sumstat2;
		  	set _sumstat;
		  		if &trtvar=0 then jointo=1;
					else if &trtvar=1 then jointo=0;
		  	run;
		  	*Merge the datasets so that can remove people;
		  	proc sql;
		  		create table _ps2 as
				select a.*, b.minPS, b.maxPS
				from _ps as a
				left join _sumstat2 as b
				on a.&trtvar=b.jointo;
				quit;
		   	*Indicate which values to delete & remove those individuals;
		  	data _ps3 (where=(delete=0));
		  	set _ps2;
				if &trtvar=0 & ps<=minPS then delete=1;
					else if &trtvar=1 & ps>=maxPS then delete=1;
					else delete = 0;
		  	run;

			%*Subset the analysis cohort to those in the trimmed sample;
			proc sql;
				create table _anacohort_trim&g as
				select *
				from _anacohort&g
				where _id in (select distinct _id from _ps3)
				;
				quit;
				
			%*At this point, we want to retain the number of individuals retained within the gestational
			age strata, by treatment;
			
			%*Count the number of rows within each treatment;
			proc freq data= _anacohort_trim&g noprint;
				table &trtvar / out=_ga_dist_&g;
			run;
			
			%*Modify the dataset so that we have the ga strata and the round;
			data __ga_dist_&g;
			set _ga_dist_&g;
				ga_strat = &g;
				round = &i;
			run;
		
			%*Output the frequency distributions;
			%IF &i = 1 and &g=1 %THEN %DO;
				data ana.&outds_dist;
				set __ga_dist_&g;
				run;
			%END;
	
	      	%ELSE %DO; 
				proc append base=ana.&outds_dist data=__ga_dist_&g; 
				run; 
			%END;

			%*Output counts;
			proc sql noprint;
				select count(distinct _id) into :num_preg_pretrim from _anacohort&g;
				select count(distinct _id) into :num_preg_trim from _anacohort_trim&g;
				quit;
			%put Number of pregnancies prior to trimming in round &i GA strata &g : &num_preg_pretrim ;
			%put Number of pregnancies after trimming in round &i GA strata &g : &num_preg_trim ;

	      	%PUT STEP 3c: RE-FIT PS MODEL & CREATE STANDARDIZED IPT WEIGHTS;
			%*Marginal probabilities of treatment;
	      	proc logistic data=_anacohort_trim&g noprint; 
				model &trtvar (reference = '0')= ; 
				output out=_num(keep=_id n) p=n; 
			run;
			%*Conditional probabilities of treatment;
	      	proc logistic data=_anacohort_trim&g noprint; 
				class &psclassVars; 
				model &trtvar (reference = '0')=&psvars; 
				output out=_den(keep=_id d) p=d; 
			run;
			
			proc sort data=_num; by _id; run;
		  	proc sort data=_den; by _id; run;
		 	
		 	%*Output the PS and SMR weight information for subsequent evaluation;
	      	data __anacohort&g; 
	         	merge _anacohort_trim&g _num _den;  
				by _id;  
	         	if &trtvar = 1 then expwgt=1; 
					else expwgt=d/(1-d); 
	         	keep _id &trtvar expwgt &outcomevar;
	      	run;

	      	%PUT STEP 4: WEIGHTED 2x2 TABLE;

			proc sort data=__anacohort&g;
				by &trtvar;
			run;

			%*Weighted proc means by treatment;
			proc means data=__anacohort&g mean noprint;
			    class &trtvar;
			    var &outcomevar;
			    weight expwgt;
			    output out=means_output(where = (&trtvar ne .)) mean=weighted_mean;
			run;

			%IF &boot = 1 %THEN %DO;
				%*Delete unnecessary datasets;
				proc datasets lib=work nolist nodetails; 
					delete _anacohort&g _anacohort_trim&g __anacohort&g; 
				run; quit;
			%END;
	
		  	%PUT Summarize the datasets;

			%*First transpose the weighted mean datasets so that the risks are columns;
			proc transpose data=means_output out=means_transpose prefix=risk;
				id trt;
				var weighted_mean;
			run;

			data _rd&g;
			set means_transpose;
				e0_risk = risk0;
				e1_risk = risk1;
				riskDiff = risk1 - risk0;
				riskRatio = risk1 / risk0;
				lnriskRatio = log(risk1 / risk0);
				round = &i;
			run;

	      	%IF &i=1 %THEN %DO; 
		  		data &outds._&g; 
				set _rd&g; 
				run; 
		  	%END;
	      	%ELSE %DO; 
				proc append base=&outds._&g data=_rd&g; 
				run; 
			%END;

		%END;

		%IF &boot = 1 %THEN %DO;
      		proc datasets lib=work nolist nodetails; delete _surv _rd:; run; quit;
      	%END;

	%END;

%mend;


