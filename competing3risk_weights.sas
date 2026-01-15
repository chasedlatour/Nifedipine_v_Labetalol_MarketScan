/*
MACRO: competing3risk_weights
PURPOSE: To implement a competing risk analysis using an Aalen Johansen estimator with standardized inverse probability
of treatment and censoring weights, in the setting where there are 3 competing events that should be modeled separately.
In addition, this macro is intended to accomodate the following Bootstrap procedure:

1. Draw a random sample, with replacement, of size N.
2. Quantify the distribution of individuals across the strata of gestational age at enrollment.
3. Complete primary analyses within each strata.
	a. Estimate standardized inverse probability of treatment weights and censoring weights within the strata.
	b. Estimate the standardized risks at the end of follow-up. 
	c. Contrast those risks as the risk difference and risk ratio. Retain the risk, risk difference, and risk ratio estimates. 
4. Reweight the strata-specific risk estimates according to the gestational age distribution at enrollment, as calculated in Step 2. Re-calculate the risk difference and risk ratio using those reweighted risks. Retain the risk, risk difference, and risk ratio estimates.
5. Repeat all the prior steps 2,000 times.

IMPORTANT NOTE ON THE INVERSE PROBABILITY OF CENSORING WEIGHTS:
Generally, when fitting a logistic regression model to address censoring, best practice is to remove intervals during which an individual
experiences an outcome. However, this creates non-convergence issues in the censoring model when applied in teh context of pregnancy studies.
This is because everyone must experience _some_ outcome by the end of follow-up. To make concrete the implications of this, consider that follow-up
is discretized into 5 intervals (as done here), there will be no uncensored individuals who contribute person-time to the 5th interval because they
will have experienced an outcome by or during that interval. As such, there is no one in the 5th interval from which to predict censoring.

To deal with this, we have decided to not exclude those intervals where individuals experience the positive competing event--example: live birth.
This derives from the idea that this competing event is included in the model out of statistical necessity to not over-inflate risks, not because
we are interested in it preventing another outcome. We note that this is an area that needs additional exploration in subsequent work.

This macro is built under the assumption that the last competing event put into the macro is this positive competing event, and the IPCW
are modeled as such.

INPUTS:
- BOOT = indicator (1=yes, 0=no) as to whether we want to conduct a bootstrap or not.
- INDS = input dataset name
- STARTDT = variable name for date of start of follow-up
- EVENTDT = variable name for date of outcome of interest -- All of our outcomes are defined at the end of pregnancy, so
	default is dt_gapreg
- OUTCOMEVAR = variable where the outcome is stored
- EVENT = value(s) of OUTCOMEVAR that represent the study outcome
- CR1 = value(s) of OUTCOMEVAR that represent competing risk event 1 (previously: CRDT)
- CR2 = value(s) of OUTCOMEVAR that represent competing risk event 2 (added)
- CR3 = value(s) of OUTCOMEVAR that represent competing risk evetn 3 (added) - POSITIVE COMPETING EVENT (e.g., live birth)
- CENSORDT = variable name for date of censoring event
- PSVARS = list of variables to be included in the propensity score model (separated by spaces)
- DOVARS = list of variables to be included in the dropout model (separated by spaces), no interaction terms
- DOVARSMODEL = list of variables to be included in the dropout model, including their functional form (e.g., interactions)
- PSCLASSVARS = list of variables to be included in the class statement for the PS model
- DOCLASSVARS = list of variables to be included in the class statement fo rthe IPCW model
- TRTVAR = treatment variable name  (values 0,1)
- NUMITERATIONS = number of iterations (bootstraps)
- INITIALSEED = optional parameter for initial seed (so results will be reproducible
- OUTDS = name of output dataset, containing one record per iteration. This will be stratified by GA at index
- GACATVAR = The variable that stratifies the pregnancies according to their gestational age at index
- GACAT = Values of GACATVAR that define the gestational age categories at index
- OUTDS_DIST = Name of the output dataset for the distribution of GA at index
- OUTDS_SURV = Name of the output dataset for the cumulative incidence/survival estimates.
- SMT = INdicator (Yes=1/No=1) to include SMR weights
- IPCW = Indicator (Yes=1/No=1) to include IPC weights
*/



/*TESTING:
%let boot=0;
%let inds=ana.primary_cohort;
%let gacatvar = ga_index_cat;
%let startDT=dt_index;
%let outcomevar = preg_outcome_ptb;
%let eventDT=dt_GApreg;
%let event = 'PTB';
%let cr1='IAB';
%let cr2='SAB' 'UAB' 'SB' 'MLS';
%let cr3='TB';
%let psvars=ga_quartile age_at_index age_at_index_2
			year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			year_index2017*t2dmrx_post year_index2017*t1t2dmrx_post year_index2017*metforrx_post
			diabetes_simp*t2dmrx_post diabetes_simp*t1t2dmrx_post diabetes_simp*metforrx_post
			nausea_pre recurlos_pre obesity_post obesity_post*year_index2017
			chronichypertension_pre substance_use_pre recurlos_pre ckd_post
			thyroid_disorder_post thyroidrx_post thyroid_disorder_post*thyroidrx_post
			depressi_post anxiety_post antideprx_post benzorx_post 
			adhdrx_post teratrx_pre num_outptpnc num_outptpnc_2;
%let dovars= ;
%let psclassvars=ga_quartile year_index2017 t2dmrx_post t1t2dmrx_post metforrx_post diabetes_simp
			nausea_pre recurlos_pre obesity_post chronichypertension_pre substance_use_pre
			recurlos_pre ckd_post thyroid_disorder_post thyroidrx_post
			depressi_post anxiety_post antideprx_post benzorx_post adhdrx_post
			teratrx_pre ;
%let doclassvars= ;
%let dovarsmodel= ;
%let trtvar=trt ;
%let numiterations=1;
%let initialseed=23244;
%let outds=ana.primary_point_noipcw;
%let outds_dist=ga_dist_primary_point_noipcw;
%let outds_ps = ps_primary_point_noipcw;
*/

%macro competing3risk_weights(
		boot=0,
		inds=pregnancies, 
		gacatvar = ga_index_cat,
		startDT=dt_index, 
		outcomevar = preg_outcome_clean,
		eventDT=dt_GApreg, 
		event = 'SAB' 'UAB' ,
		cr1='IAB', 
		cr2='LBM' 'LBS' 'UDL',
		cr3= 'SB' 'MLS',
		psvars=age_at_index chronichypertension_pre, 
		dovars=age_at_index chronichypertension_pre, 
		dovarsmodel = age_at_index chronichypertension_pre, 
		psclassvars=chronichypertension_pre, 
		doclassvars=chronichypertension_pre,
		trtvar=trt, 
      	numiterations=4, 
		initialseed=23244, 
		outds=ana.comprisk_primary,
		outds_dist=ga_dist_primary,
		outds_ps = ps_primary,
		outds_surv = surv_primary,
		smr=1,
		ipcw=1
	   );



	%PUT STEP 1: PREPARE ANALYTIC COHORT DATASET;
	* 1a: create analytic variables for time to event and outcome;
	data _anacohort; 
	set &inds;
			
		%*Define the time to event;
		days = &eventDT - &startDT + 1;

		%*Create an indicator for the outocme experienced;
		if &outcomevar in ( &event ) then outcome = 1;
			else outcome = 0;
		if &outcomevar in ( &cr1 ) then comprisk1 = 1;
			else comprisk1 = 0;
		if &outcomevar in ( &cr2 ) then comprisk2 = 1;
			else comprisk2 = 0;
		if &outcomevar in ( &cr3 ) then comprisk3 = 1;
			else comprisk3 = 0;

		%*Define the outcome;
		if &outcomevar = "UNK" then event = 0;
			else if &outcomevar in ( &event ) then event = 1;
			else if &outcomevar in ( &cr1 ) then event = 2;
			else if &outcomevar in ( &cr2 ) then event = 3;
			else if &outcomevar in ( &cr3 ) then event = 4;

		if event in (1, 2, 3, 4) then combined = 1;
			else combined = 0;
	run;


	%PUT STEP2A: If bootstrapping, then resample with replacement.;
	%if &boot = 1 %then %do;
		%*Implement the bootstrap of the samples prior to the loop. This ensures that one stream of random samples are used.;
		%PUT STEP 2: CREATE BOOTSTRAP SAMPLES;
	    proc surveyselect noprint data=_anacohort out=_anacohorta(rename=(replicate=b) drop=numberhits) 
	    	seed=&initialseed method=urs samprate=1 outhits rep=&numIterations; 
	%end;
	%else %do;
		data _anacohorta;
		set _anacohort;
		run;
	%end;

	%PUT STEP2B: JITTER EVENT TIMES TO AVOID TIES;
    %LET jitterFlag=1;
    %LET i=0;
    
    %IF &boot = 1 %THEN %DO;
	    %DO %WHILE (&jitterFlag=1);
	      	Title "Jitter Iteration &i";
	      	proc sql noprint; 
	         	select case when max(numEvents)>1 then 1 else 0 end into :jitterFlag 
				from (select distinct days, count(*) as numEvents 
						from _anacohorta 
						where combined=1 
						group by b, days); /*Group by bootstrap replicate because we only care about ties within each repeated sample*/
	      		quit;
	
	      	%IF &jitterFlag=1 %THEN %DO;
	         	%LET i = %EVAL(&i+1);
	         	data _anacohorta; 
				set _anacohorta; 
					call streaminit(123+&i); 
					days=days+rand("uniform")*.0055-.00275; 
				run;
	      	%END;
	    %END;
    %END;
    %ELSE %DO;
    	data _anacohorta; 
		set _anacohorta; 
			call streaminit(123+&i); 
			days=days+rand("uniform")*.0055-.00275; 
		run;
    %END;

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

			%IF &smr=1 %THEN %DO;

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
					%*where combined=0; %*CDL: ADDED 6.09.2025 -- Pooled logistic regression not supposed to be fit on events.;
					model &trtvar (reference = '0')= ; 
					output out=_num(keep=_id n) p=n; 
				run;
				%*Conditional probabilities of treatment;
		      	proc logistic data=_anacohort_trim&g noprint; 
					%*where combined=0; %*CDL: ADDED 6.09.2025 -- Pooled logistic regression not supposed to be fit on events.;
					class &psclassVars; 
					model &trtvar (reference = '0')=&psvars; 
					output out=_den(keep=_id d) p=d; 
				run;
				
				proc sort data=_num; by _id; run;
			  	proc sort data=_den; by _id; run;
				
				%*Output the propensity scores if not bootstrapping;
				%IF &boot = 0 %THEN %DO;
					data ana.&outds_ps._num_&g; 
			         	merge _anacohort_trim&g _num _den; 
						by _id; 
						%*Fit SMR weights. Originally, fit stabilized IPTW.;
						if &trtvar = 1 then expwgt=1;
							else expwgt = d / (1-d);
			         	keep idxpren &trtvar n d expwgt;
			      	run;
				%END;

			 	
			 	%*Output the PS and SMR weight information for subsequent evaluation;
		      	data __anacohort&g; 
		         	merge _anacohort_trim&g _num _den;  
					by _id;  
		         	if &trtvar = 1 then expwgt=1; 
						else expwgt=d/(1-d); 
		         	keep _id &trtvar expwgt event outcome comprisk1 comprisk2 comprisk3 combined days event &dovars;
		      	run;

			%END;

			%ELSE %DO;

				%*Count the number of rows within each treatment;
				proc freq data= _anacohort&g noprint;
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

				%*Assign the exposure weights as 1;
				data __anacohort&g;
				set _anacohort&g;
					expwgt=1;
					keep _id &trtvar expwgt event outcome comprisk1 comprisk2 comprisk3 combined days event &dovars;
				run;
			%END;

	      	%PUT STEP 4: STABILIZED INVERSE PROBABILITY OF CENSORING WEIGHTS;

			%*Get maximum follow-up over the included pregnancies;
			proc sql noprint;
				select max(days)
				into :max_follow
				from __anacohort&g;
				quit;

			%*CDL: REVISED 12.4.2025 to discretize weekly;
			data ___anacohort&g; 
			set __anacohort&g(rename=(days=days1 outcome=outcome1 comprisk1=comprisk11 comprisk2=comprisk21 comprisk3=comprisk31 combined=combined1));

				%*Make sure there are no 0 days values;
				if days1 = 0 then days1 = 0.001;

				cat = ceil(ceil(days1)/7); *Round up by week;

				%*Now create a long dataset where each row represents a week of follow-up for each person;
				do j=1 to cat;

					%*Create the start time for the 1-week interval;
					in = (j-1) * 7; %*start of the interval;

					if in < days1 and days1 <= j*7 then do;
						days = days1;
						if combined1 = 0 then drop=1;
							else drop=0;
						outcome = outcome1;
						comprisk1 = comprisk11;
						comprisk2 = comprisk21;
						comprisk3 = comprisk31;
						combined = combined1;
					end;
					else do;
						days = j*7;
						drop = 0;
						outcome=0; comprisk1=0; comprisk2=0; comprisk3=0; combined=0;
					end;
					output;

				end;

	         	keep _id in days outcome comprisk1 comprisk2 comprisk3 combined drop expwgt &trtvar &dovars;
	      	run;

			%IF &boot = 1 %THEN %DO;
				%*Delete unnecessary datasets;
				proc datasets lib=work nolist nodetails; 
					delete _anacohort&g __anacohort&g; 
				run; quit;
			%END;

			%*Now fit the pooled logistic regressions to fit the stabilized inverse probability of censoring weights;
	     	proc logistic data=___anacohort&g noprint; 
				where combined = 0 ; %*CDL: ADDED 6.09.2025 -- Pooled logistic regression not supposed to be fit on events.;
				effect spl = spline(in / details naturalcubic basis=tpf(noint) knotmethod=percentiles(5)); %*CDL: ADDED 12.4.2025 restricted cubic spline;
				class &doclassvars / param=ref;  %*CDL: REMOVED in from class statement because will not converge;
	        	model drop = &trtvar spl; 
				output out=num(keep=_id in dn) p=dn; 
			run;
	      	proc logistic data=___anacohort&g noprint; 
				where combined = 0 ; %*CDL: ADDED 6.09.2025 -- Pooled logistic regression not supposed to be fit on events.;
				effect spl = spline(in / details naturalcubic basis=tpf(noint) knotmethod=percentiles(5)); %*CDL: ADDED 12.4.2025 restricted cubic spline;
				class &doclassvars / param=ref; %*CDL: REMOVED in from class statement because will not converge;
	        	model drop = &trtvar spl &dovarsmodel; 
				output out=denom(keep=_id in dd) p=dd; 
			run;

	      	proc sort data=___anacohort&g; by _id in; run;
	      	proc sort data=denom; by _id in; run;
	      	proc sort data=num; by _id in; run;

	      	data ____anacohort&g; 
			merge ___anacohort&g denom num; 
				by _id in;
	         	retain num denom lastNum lastDenom;

	         	if first._id then do; 
					num=1; 
					denom=1; 
				end;
	         	else do; 
					num=num*lastNum; 
					denom=denom*lastDenom; 
				end;

	         	lastNum=dn; 
				lastDenom=dd;
				%IF &ipcw = 1 %THEN %DO;
	         		dowgt = num / denom;
				%END;
				%ELSE %DO;
					dowgt = 1;
				%END;
	         	wgt = expwgt * dowgt;
	         	keep _id &trtvar in days outcome comprisk1 comprisk2 comprisk3 combined wgt dowgt;
	     	run;

			%IF &boot = 1 %THEN %DO;
		      	proc datasets lib=work nolist nodetails; 
					delete ___anacohort&g _ps _sumstat _sumstat2 _ps2 _ps3 _num _den _quintiles num denom; 
				run; quit;
			%END;

	      	%PUT STEP 5: RUN ANALYSIS; 

			%*Calculate the survival for ANY event;
	      	proc phreg data=____anacohort&g noprint; 
				strata &trtvar;
	        	model days*combined(0)= /entry=in; 
				weight wgt;
	        	baseline out=_combined(rename=(survival=s_combined) keep=&trtvar days survival) survival=_ALL_/method=ch; 
			run;

			%*Calculate the survival for the primary outcome;
	      	proc phreg data=____anacohort&g noprint; 
				strata &trtvar;
	         	model days*outcome(0)= /entry=in; 
				weight wgt;
	         	baseline out=_outcome(rename=(survival=s_outcome) keep=&trtvar days survival) survival=_ALL_/method=ch; 
			run;

			%*Calculate the survival for competing risk 1;
			proc phreg data=____anacohort&g noprint; 
				strata &trtvar;
	         	model days*comprisk1(0)= /entry=in; 
				weight wgt;
	         	baseline out=_comprisk1(rename=(survival=s_comprisk1) keep=&trtvar days survival) survival=_ALL_/method=ch; 
			run;

			%*Calculate the survival for competing risk 2;
			proc phreg data=____anacohort&g noprint; 
				strata &trtvar;
	         	model days*comprisk2(0)= /entry=in; 
				weight wgt;
	         	baseline out=_comprisk2(rename=(survival=s_comprisk2) keep=&trtvar days survival) survival=_ALL_/method=ch; 
			run;
			
			%*Calculate the survival for competing risk 3;
			proc phreg data=____anacohort&g noprint; 
				strata &trtvar;
	         	model days*comprisk3(0)= /entry=in; 
				weight wgt;
	         	baseline out=_comprisk3(rename=(survival=s_comprisk3) keep=&trtvar days survival) survival=_ALL_/method=ch; 
			run;
	
			%IF &boot = 1 %THEN %DO;
		      	proc datasets lib=work nolist nodetails; 
					delete ____anacohort&g; 
				run; quit;
			%END;
	      	proc sort data=_combined; by &trtvar days; run;
	      	proc sort data=_outcome;  by &trtvar days; run;
	      	proc sort data=_comprisk1; by &trtvar days; run;
			proc sort data=_comprisk2; by &trtvar days; run;
			proc sort data=_comprisk3; by &trtvar days; run;

	      	data _surva; %* (keep=&trtvar rate); 
	        merge _combined _outcome _comprisk1 _comprisk2 _comprisk3; 
				by &trtvar days;
	        	retain cum_outcome 0 cum_comprisk1 0 cum_comprisk2 0 cum_comprisk3 0 olds_combined olds_outcome olds_comprisk1 olds_comprisk2 olds_comprisk3 rate;
				if first.&trtvar then rate = .; 
	         	if days=0 then do; 
	            	olds_combined=1; olds_outcome=1; olds_comprisk1=1; olds_comprisk2 = 1; olds_comprisk3 = 1;
	            	outcome=0; comprisk1=0; comprisk2=0; comprisk3=0;
					cum_outcome=0; cum_comprisk1=0; cum_comprisk2=0; cum_comprisk3=0;
					e_outcome=0; e_comprisk1=0; e_comprisk2=0; e_comprisk3=0;
	         	end;

				if days >= &max_follow and rate = . then do; 
					rate=1-cum_outcome;
/* 					output; */
				end;

				%*Implement algorithm described by Steve Cole in EPID 722;
	         	if s_outcome ne . then do; 
	            	h_outcome=-log(s_outcome)--log(olds_outcome);  %*Calculate the (unscaled by time) hazard;
	            	e_outcome=olds_combined*h_outcome;  %*Calculate lagged survival by the hazard;
	            	cum_outcome=cum_outcome+e_outcome; %*Then add it to the risk up to that point;
	         	end;

	         	if s_comprisk1 ne . then do; 
	            	h_comprisk1=-log(s_comprisk1)--log(olds_comprisk1); 
	            	e_comprisk1=olds_combined*h_comprisk1; 
	            	cum_comprisk1=cum_comprisk1+e_comprisk1; 
	         	end;  

				if s_comprisk2 ne . then do; 
	            	h_comprisk2=-log(s_comprisk2)--log(olds_comprisk2); 
	            	e_comprisk2=olds_combined*h_comprisk2; 
	            	cum_comprisk2=cum_comprisk2+e_comprisk2; 
	         	end; 
	         	
	         	if s_comprisk3 ne . then do; 
	            	h_comprisk3=-log(s_comprisk3)--log(olds_comprisk3); 
	            	e_comprisk3=olds_combined*h_comprisk3; 
	            	cum_comprisk3=cum_comprisk3+e_comprisk3; 
	         	end; 

	         	combined=cum_outcome+cum_comprisk1+cum_comprisk2+cum_comprisk3;
	         	if s_combined ne . then olds_combined=s_combined;
	         	if s_outcome ne . then olds_outcome=s_outcome;
	         	if s_comprisk1 ne . then olds_comprisk1=s_comprisk1;
				if s_comprisk2 ne . then olds_comprisk2=s_comprisk2;
				if s_comprisk3 ne . then olds_comprisk3=s_comprisk3;

	         	if last.&trtvar then do; 
					if rate = . then do;
						rate = 1 - cum_outcome;
/* 						output; */
					end;
				end;
	      	run;
	      	
	      	%*Output the survival dataset for non-bootstrapped analyses;
	      	%IF &boot = 0 and &outds_surv ne NA %THEN %DO;
		      	data &outds_surv._&g;
		      	set _surva;
		      	run;
	      	%END;
	      	
	      	%*Output the last trt row for subsequent analyses;
	      	data _surv (keep=&trtvar rate);
	      	set _surva;
	      		by &trtvar days;
	      		
	      		if last.&trtvar then output;
	      	run;
	      		

			%IF &boot = 1 %THEN %DO;
			  	proc datasets lib=work nolist nodetails; 
					delete _combined _outcome _comprisk:; 
				run; quit;
			%END;
		
		  	%PUT Summarize the datasets;

	      	data _rd&g(drop=&trtvar); 
	        merge _surv(where=(&trtvar=0) rename=(rate=e0_surv))  
	               _surv(where=(&trtvar=1) rename=(rate=e1_surv)); 
				e0_risk = 1-e0_surv;
				e1_risk = 1-e1_surv;
				riskDiff = e1_risk - e0_risk;
				riskRatio = e1_risk / e0_risk;
				lnriskRatio = log( e1_risk / e0_risk );
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


