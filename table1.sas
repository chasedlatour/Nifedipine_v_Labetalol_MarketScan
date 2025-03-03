/*****************************************************************************************************/
/* Program: /nearline/files/datasources/references/macros/table1.sas                                 */
/* Purpose: Produce output for a typical Table 1: Baseline Covariates                                */
/*                                                                                                   */
/* Created on: November 5, 2011                                                                      */
/* Created by: Virginia Pate                                                                         */
/*                                                                                                   */
/* Inputs: INDS = one or two level name of the input dataset to be used to create the table; data    */
/*                set should have one record per ID                                                  */
/*                                                                                                   */
/*         COLVAR  = name of the variable that identifies the analysis group for the given record;   */
/*                   the final table will have one column for each distinct value of COLVAR          */
/*                                                                                                   */
/*         ROWVARS = a list of the variables to be summarized in the table; each variable will have  */
/*                   one (or multiple) rows in the final table.  The variables should be listed in   */
/*                   the order in which they should appear in the table.  Variable names should be   */
/*                   separated with a space                                                          */
/*                                                                                                   */
/*         WGTVAR = name of variable to be used for weighting (leave missing for unweighted)         */
/*                                                                                                   */
/*         MAXLEVELS = the highest number of unique values that a numeric variable can have a be     */
/*                     classified as a categorical variable - if the variable has >&MAXLEVELS        */
/*                     unique values, it will be treated as a continuous variable                    */
/*                                                                                                   */
/*         PCTONLY = Y/N to indicate whether to include Ns(N) or only percentages(Y) for categorical */
/*                   variables (default = N, meaning N(%) will be presented)                         */
/*                                                                                                   */
/*         CONTSTAT = continuous variable statistic (default = mean)                                 */
/*                    mean: mean(SD), median: median(IQR)                                            */
/*                                                                                                   */
/*         SMD_CAT = method for calculating SMD for categorical variables with >2 levels             */
/*                   default = level -- calculates SMD for each individual level                     */
/*                   to calculate overall SMD, specify overall                                       */
/*                   NOTE: to use overall, you must ALSO include the %STDDIFF macro                  */
/*                                                                                                   */
/*         OUTFILE = the file name for the final table - this will be output as an RTF file to the   */
/*                   output subdirectory (&OutPath) previously defined                               */
/*                                                                                                   */
/* Details: This macro produces a typical Table 1.  The final table will have one column for each    */
/*          distinct value of COLVAR, with the largest group in the first column, plus a total       */
/*          column as the last column.  Continuous variables are presented as mean(SD).  Two level   */
/*          categorical variables are presented as N(%) for the larger group (which is specified in  */
/*          the Characteristics column).  Categorical variables with more than two levels are        */
/*          presented as N(%) for each distinct level, listed in ascending order by value of the     */
/*          variable.                                                                                */
/*                                                                                                   */
/* Example macro call: %table1(inds=der.cohort_bl,                                                   */
/*                             colvar=trt,                                                           */
/*                             rowvars=age sex race chf hpt cancer,                                  */
/*                             maxlevels=3,                                                          */
/*                             outFile=Table 1)                                                      */
/*                                                                                                   */
/* Updates:                                                                                          */
/*    9/12/2018 - corrected SMD calculation for proportions!                                         */
/*    10/2/2019 - added options PCTONLY and CONTSTAT                                                 */
/*    5/5/2021  - added call to %STDDIFF to calculate overall SMD for categorical variables          */
/*    3/9/2022  - moved SDIFF2 variable to top line of categories with >2 levels                     */
/*****************************************************************************************************/


%macro table1(inds=, colVar=, rowVars=, wgtVar=, maxLevels=4, pctOnly=N, contStat=mean, 
				  smd_cat=level, outfile=, title=);

   ods noptitle escapechar='~'; options minoperator;

/***************************************************************************/
/**  STEP 1 - PREPARE DATASET & TABLE, USING BYGROUP VARIABLE AS COLUMNS  **/
/***************************************************************************/

/*** STEP 1A - GET DATASET AND LIBNAME ***/
/*Parse libname and dataset name from &inds.  If only one level, assign libname=WORK*/
   %IF %INDEX(&inds, .)>0 %THEN %DO;
      %LET libname = %UPCASE(%SCAN("&inds",1,"."));
      %LET ds = %UPCASE(%SCAN("&inds",2,"."));
   %END;/*end libname specified*/ %ELSE %DO;
      %LET libname = WORK;
      %LET ds = %UPCASE(&inds);
   %END;/*end libname not specified*/
   %LET colVar = %UPCASE(&colVar);

/*** STEP 1B - GET BY GROUP (COLUMN) DATA ***/
   /*Get BY GROUP variable information: variable type, label, format, distinct levels, etc.*/
   proc contents data=&libname..&ds noprint out=metadata(keep=name type label format length); run;

   /*The distinct levels of BYGROUP will be used to create variable names (and later table columns)*/
   /*Therefore, we need to apply a format if it has one and they must be character variables*/
   /*Also, get N for each level of BYGROUP to use in column headers*/
   data _null_; set metadata(where=(upcase(name)="&colVar"));
      if type=2 and format='$' then format=cats(format, put(length,8.),'.'); 
         else if type=2 and format='' then format=cats('$', put(length,8.),'.');
         else if type=1 and format='' then format='8.'; 
         else if format ne '' then format=cats(format,'.');
      if type=2 then call symput('colType','char'); else call symput('colType','num');
      call symput('colFMT', format); run;

   proc sql noprint;  select distinct 
           put(&colVar, &colFMT.), 
           case when count(*)<0 then 'NTSR' else put(count(*),comma12.) end as colN
      into :colLabel1-:colLabel10, :colN1-:colN10
      from &inds group by put(&colVar,&colFmt.);
      %LET NumCol = &SqlObs;
   quit;
   
   %DO l=1 %TO &numCol; 
        %IF &&colN&l = NTSR %THEN %DO;
           %PUT ********* All columns must include at least 11 people in order to run table *********;
          %ABORT;
        %END;     
        %IF "%SYSFUNC(COMPRESS(&&colLabel&l,,kns))" ^= "&&colLabel&l" OR 
          %SUBSTR(&&colLabel&l,1,1) IN 0 1 2 3 4 5 6 7 8 9 _ %THEN %DO;
           %PUT ********* All formatted column variable values or, if no format is applied *********;
           %PUT ********* to the column variable then unformatted column variable values, *********;
           %PUT ********* must be a valid SAS variable name *********; 
         %ABORT;
     %END;  
   %END;

   proc sql noprint; 
     select distinct  put(&colVar,&colFMT.), 
            case when count(*)<0 then 'NTSR' else  
              %IF &wgtVar ^=  %THEN put(sum(&wgtvar), comma12.0); %ELSE put(count(*), comma12.); end as colN
         into          :colLabel1-:colLabel&numCol., :colN1-:colN&numCol.
        from &inds group by put(&colVar,&colFMT.) order by colN desc;

     select case when count(*)<0 then 'NTSR' else %IF &wgtVar ^=  %THEN strip(put(sum(&wgtvar), comma12.0)); %ELSE strip(put(count(*), comma12.)); end
         into :totN from &inds;     quit;

   %DO l=1 %TO &numCol; %LET col&l = %SYSFUNC(TRANSLATE(&&colLabel&l,_,%STR( ))); %END;

/*** STEP 1C - CREATE SHELL ***/ 
   /*Create shell dataset using variable names and labels created above*/
   data final; length row $250 %DO i=1 %TO &NumCol; &&col&i %END; total $100 order rowOrder 8;
      label row = 'Characteristic'  total = "Total ~n N=&totN"
            %DO i=1 %TO &NumCol; &&col&i = "&&colLabel&i ~n N=&&colN&i" %END;;              
      set _null_;     run;

/****************************************************/
/**  STEP 2 - PROCESS VARIABLES FOR ROWS OF TABLE  **/
/****************************************************/

   /*Now process each variable in the variable list - VARLIST will serve as the ROWS*/     
   /*Loop through each variable listed in &varlist*/
   %LET NumRow = %SYSFUNC(countw(&rowVars));
   %DO r=1 %TO &NumRow;  %LET rowVar = %SCAN(&rowVars,&r);

/*** STEP 2A - GET METADATA ON VARIABLE ***/
   /*Get variable LABEL to serve as the value for the Covariate column of Table 1*/
   /*Get variable FORMAT to ensure output dataset is correct*/
   data _null_; set metadata(where=(upcase(name)="%upcase(&rowVar)"));
      if type=2 and format='$' then format=cats(format, put(length,8.),'.'); 
         else if type=2 and format='' then format=cats('$', put(length,8.),'.');
         else if type=1 and format='' then format='8.'; 
         else if format ne '' then format=cats(format,'.');
      if type=2 then call symput('rowType','char'); else call symput('rowType','num');
      if label='' then label=propcase(name);
      call symput('rowFMT', format);  call symput('rowLabel', label); run;

   /*Determine if variable is CONTINUOUS OR CATEGORICAL based on the variable type and the number of distinct values*/
   proc sql noprint; 
      select %IF &rowFMT =  %THEN  
         /*no FMT*/ count(distinct &rowVar)               ;/*end no FMT*/ %ELSE  
         /*FMT   */ count(distinct put(&rowVar,&rowFMT.)) ;/*end FMT*/ 
               into :rowLevel from &inds; quit;

   %IF &rowLevel = 1 %THEN %DO;           %LET type = ;    %END; %ELSE
   %IF &rowType = char %THEN %DO;         %LET type = cat; %END;/*end char*/ %ELSE 
   %IF (&rowLevel<=&maxLevels) %THEN %DO; %LET type = cat; %END;/*end <MaxLevels*/ %ELSE %DO; 
                                          %LET type = con; %END;/*end ow*/
                
   /*For BINARY variables, use yes if y/n variable, otherwise use the smaller group as the group to report*/
   /*For variables with >2 levels, report each level*/
   proc sql noprint;
      create table temp as select distinct   
            %IF &rowFmt =  %THEN &rowVar; %ELSE put(&rowVar., &rowFMT.); as &rowVar, 
               count(*) as numRecs from &inds group by 
            %IF &rowFmt =  %THEN &rowVar; %ELSE put(&rowVar., &rowFMT.); ;

      select distinct &rowVar into :rowLevel1-:rowLevel100 from temp; quit;

      %IF &rowLevel = 2 %THEN %DO; /*IF BINARY*/ /* if yes/no variable, set reference to yes*/ 
        %IF (&rowLevel1=0 & &rowLevel2=1) OR (&rowLevel1=1 & &rowLevel2=0) OR
             (&rowLevel1=Y & &rowLevel2=N) OR (&rowLevel1=N & &rowLevel2=Y) OR
             (&rowLevel1=No & &rowLevel2=Yes) OR (&rowLevel1=Yes & &rowLevel2=No) %THEN %DO;
           %LET rowRefLabel= ;
           data _null_; set temp; if strip(&rowVar) in ('1' 'Y' 'Yes') then do;
                call symput('rowRefValF', strip(&rowVar)); call symput('rowRefVal', strip(&rowVar)); end;
           run;
        %END;
       
        %ELSE %DO; /* otherwise use the smaller group */ proc sql noprint; 
            select   &rowVar,     &rowVar,              &rowVar
              into  :rowRefVal,  :rowRefValF,          :rowRefLabel 
              from temp having numRecs=min(numRecs); 
        %END;
      %END;/*end rowLevel=2*/ quit;

/***********************************************************************************/
/**  STEP 3 - GET MEANS, SD, & STANDARDIZED DIFFERENCES FOR CONTINUOUS VARIABLES  **/
/***********************************************************************************/
   %IF &type =  %THEN %DO;
      data &rowVar; set temp; 
         length row $250 %DO i=1 %TO &NumCol; &&col&i %END; total $100 order rowOrder 8;
         row="%SYSFUNC(strip(&rowLabel)), " || strip(&rowVar); order=&r; rowOrder=1; 
         %DO j=1 %TO &numCol; &&col&j = "%SYSFUNC(strip(&&colN&j)) (100%)"; %END;
         total = "%SYSFUNC(strip(&totN)) (100%)";
         drop numRecs &rowVar; run;
      data final; set final &rowVar; run;
   %END; /*end missing type (rowLevel=1)*/
   %ELSE %IF &type = con %THEN %DO; Title "&rowVar by &colVar";
      proc means data=&inds mean std vardef=wdf %IF &contStat=median %THEN median q1 q3;;
         class &colVar;  var &rowVar; %IF &wgtvar ^=  %THEN %DO; weight &wgtvar; %END;
         %IF &colFMT ne  %THEN %DO; format &colVar &colFMT.; %END; /*end cont, no FMT*/
         ods output summary=&rowVar;          run;

      proc means data=&inds mean std vardef=wdf %IF &contStat=median %THEN median q1 q3;; Title "&rowVar Overall";
         %IF &wgtvar ^=  %THEN %DO; weight &wgtvar; %END;
         var &rowVar;  ods output summary=&rowVar.overall;   run;

      proc datasets lib=work nolist nodetails; modify &rowVar; format &colVar; run;quit;

      data &rowVar.2 (keep=group &rowVar &colVar &rowVar._mean &rowVar._stddev); set &rowVar. (in=bygroup ) &rowVar.overall (in=all ); length group &rowVar $100 ;
         if bygroup then group =                                  %IF &colTYPE = num %THEN %DO; %IF &colFMT =  %THEN 
            /*num, no FMT*/ "grp_" || strip(put(&colVar,8.));    %ELSE 
            /*num,    FMT*/ strip(put(&colVar, &colFMT.));       %END; /*end num*/  %ELSE %DO; %IF &colFMT =  %THEN 
            /*char,no FMT*/ strip(&colVar);                      %ELSE 
            /*char,   FMT*/ strip(put(&colVar, &colFMT.));       %END;/*end char*/ ; 
         else if all then group = 'Total';
         &rowVar = %IF &contStat=mean %THEN strip(put(&rowVar._mean,8.1)) || '(' || strip(put(&rowVar._stddev,8.2)) || ')';
             %ELSE %IF &contStat=median %THEN strip(put(&rowVar._median,8.1)) || '(' || strip(put(&rowVar._q1,8.1)) || '-' || strip(put(&rowVar._q3,8.1)) || ')';;
      run;

      proc transpose data=&rowVar.2 out=&rowVar.T(drop=_name_ _label_);                                  var &rowVar;         id group;   run;
      proc transpose data=&rowVar.2 out=meanT(drop=_name_ _label_) prefix=mean; where group ne 'Total';  var &rowVar._mean;   id group; run;
      proc transpose data=&rowVar.2 out=stdT(drop=_name_ _label_)  prefix=std;  where group ne 'Total';  var &rowVar._stddev; id group; run;

     %IF &NumCol>1 %THEN %DO; 
         data &rowVar.3; merge meanT stdT; length %DO c=2 %TO &numCol; sdiff&c %END; $12; 
            %DO c=2 %TO &numCol;
               d&c=round(abs((mean&col1-mean&&col&c)/sqrt((std&col1*std&col1+std&&col&c*std&&col&c)/2)),0.001);
               sdiff&c = put(d&c, 8.3); %END;
            keep %DO c=2 %TO &numCol; sdiff&c %END;; 
            label %DO c=2 %TO &numCol; sdiff&c ="Stdz Diff ~n Col 1 vs Col &c" %END;; run; %END;

      data &rowVar.F; merge %IF &numCol>1 %THEN &rowVar.3; &rowVar.T ; length row $250; 
            row= %IF &contStat=mean %THEN "%SYSFUNC(strip(&rowLabel)), mean(SD)"; %ELSE %IF &contStat=median %THEN "%SYSFUNC(strip(&rowLabel)), median(IQR)";;
            order=&r; rowOrder=1; run;

      data final;  set final &rowvar.F; run;
   %END;/*end con*/
/*************************************************************************************/
/**  STEP 4 - GET PROPORTIONS & STANDARDIZED DIFFERENCES FOR CATEGORICAL VARIABLES  **/
/*************************************************************************************/
   %ELSE %IF &type = cat %THEN %DO;
      proc freq data=&inds;  Title "&rowVar by &colVar";
         %IF &wgtvar ^=  %THEN %DO; weight &wgtvar; %END;
         tables &colVar * &rowVar / missing;  %IF &rowFMT ^=  %THEN %DO; format &rowVar &rowFMT.; %END;/*end no FMT*/   
         ods output crosstabfreqs=&rowVar.(
               keep   = &colVar &rowVar rowPercent Percent Frequency _type_ 
               rename = (&rowVar=&rowVar.2)                                
               %IF &rowLevel = 2 %THEN where = ( strip(put(&rowVar.2,&rowFMT.)) = "&rowRefValF");
               ) ; run;
      data &rowVar.2 (keep=group &rowVar &colVar %IF &numCol>1 %THEN rowpct; %IF &rowLevel>2 %THEN row rowOrder &rowVar.2;); set &rowVar; 
         length %IF &rowLevel>2 %THEN row $250; group &rowVar $100; %IF &rowLevel>2 %THEN %DO; rowOrder=_N_; %END;/*end rowLevel>2*/;
         if _type_ = '11' then do; %IF &numCol>1 %THEN %DO; rowpct=rowpercent/100; %END;
            group =                                            %IF &colTYPE = num %THEN %DO; %IF &colFMT =  %THEN 
               /*num, no FMT*/ "grp_" || strip(put(&colVar,8.));  %ELSE 
               /*num,    FMT*/ put(&colVar, &colFMT.);            %END;/*end num*/   %ELSE %DO; %IF &colFMT =  %THEN 
               /*char,no FMT*/ &colVar;                           %ELSE 
               /*char,   FMT*/ put(&colVar, &colFMT.);            %END;/*end char*/ ; 
            if Frequency>=0 then &rowVar = %IF &pctOnly=N %THEN strip(put(Frequency,comma12.)) || ' (' ||; strip(put(rowpercent,8.1)) %IF &pctOnly=N %THEN || '%)' ;;
            else &rowVar = 'NTSR';
         end; else if _type_ = '01' then do;
         group = 'Total';
         if Frequency>=0 then &rowVar = %IF &pctOnly=N %THEN strip(put(Frequency,comma12.)) || ' (' ||; strip(put(percent,8.1)) %IF &pctOnly=N %THEN || '%)';; else &rowVar = 'NTSR';
         end;  else delete;
                    
         /*For multi-level variables, create one row per value of ROWVAR */
         %IF &rowLevel > 2 %THEN %DO;                                           %IF &rowType=char %THEN %DO; %IF &rowFMT =  %THEN %DO;
            /*char, no FMT*/row = '~R/RTF"\tab" ' || strip(&rowVar.2);                                      %END; %ELSE %DO; 
            /*char, FMT   */row = '~R/RTF"\tab" ' || strip(put(&rowVar.2,&rowFMT.));                %END; %END;/*end char*/ %ELSE %DO; %IF &rowFMT = %THEN %DO; 
            /*num, no FMT */row = '~R/RTF"\tab"Group ' || strip(put(&rowVar.2,8.));       %END; /*end num, no FMT*/ %ELSE %DO; 
            /*num, FMT    */row = '~R/RTF"\tab" ' || strip(put(&rowVar.2, &rowFMT.)); %END; /*end num, FMT*/ %END; /*end num*/ %END;/*end rowLevel>2*/
      run;

      proc datasets lib=work nolist nodetails; modify &rowVar.2; format &colVar; run;quit;

      %IF &rowLevel>2 %THEN %DO; proc sort data=&rowVar.2; by row; run; %END;/*end rowLevel>2*/
      proc transpose data=&rowVar.2 out=&rowVar.T;  %IF &rowLevel>2 %THEN %DO; by row; %END; var &rowVar; id group; run;
      %IF &NumCol>1 %THEN %DO;
         proc transpose data=&rowVar.2 out=meanT prefix=mean; where group ne 'Total'; %IF &rowLevel>2 %THEN %DO; by row; %END;/*end rowLevel>2*/
               var rowpct; id group; run;*9/12/18 - corrected to use proportion rather than percentage in sdiff calculation!;
         data &rowVar.3; set meanT; length %DO c=2 %TO &numCol; sdiff&c %END; $12; 
            %DO c=2 %TO &numCol; 
            d&c=round(abs((mean&col1-mean&&col&c)/sqrt(abs((mean&col1*(1-mean&col1)+mean&&col&c*(1-mean&&col&c)))/2)),0.001);
               sdiff&c = put(d&c, 8.3); %END; 
            keep %DO c=2 %TO &numCol; sdiff&c %END; %IF &rowLevel>2 %THEN row; ; 
            label %DO c=2 %TO &numCol; sdiff&c="Stzd Diff ~n Col 1 vs Col &c" %END;; run; %END;

      %IF &rowLevel=2 %THEN %DO; 
         data &rowVar.F; merge %IF &NumCol>1 %THEN &rowVar.3; &rowVar.T(drop=_name_); length row $250; rowOrder=1; order=&r;
              row="%SYSFUNC(strip(&rowLabel))" %IF "&rowRefLabel" ^= "" %THEN || ", " || "%SYSFUNC(strip(&RowRefLabel))";; run;
      %END;/*end rowLevel=2*/ %ELSE %DO;
         proc sort data=&rowVar.2 out=rowOrder; by row rowOrder; run; data rowOrder2; set rowOrder; by row; if first.row; keep row rowOrder;run;
         data &rowVar.F %IF &smd_cat=overall %THEN (drop=sdiff2);; 
				merge &rowVar.T(drop=_name_) %IF &NumCol>1 %THEN &rowVar.3; rowOrder2; by row; order=&r; run; proc sort data=&rowVar.F; by rowOrder; run;

         data filler; length row $250; order=&r-0.5; rowOrder=1; row = "~S={font_weight=bold}%SYSFUNC(strip(&rowLabel))"; run;
         %IF &smd_cat = overall %THEN %DO; *if overall SMD requested, run correct SMD for >2 levels;
   			%stddiff(inds=&inds, groupvar=&colVar, wtvar=&wgtvar, charvars=&rowVar, stdfmt=8.3, outds=stdiff_&rowVar)
				data filler; merge filler stdiff_&rowVar; sdiff2=put(abs(input(stddiff,8.4)),8.3); drop stddiff; run; 
         %END;
      %END;/*end rowLevel>2*/
                         
      data final; set final %IF &rowLevel>2 %THEN filler; &rowvar.F; run;
   %END;/*end cat*/
   %END;/*end numRow loop*/

   %IF &outfile ^=  %THEN %DO;
      ods noptitle; ods escapechar='~'; title "&title";
      ods rtf file="&OutPath./&OutFile..rtf" style=vpmin startpage=never bodytitle;
         proc print data=final(drop=order rowOrder) noobs label; run;
      ods rtf close;
   %END;

%mend table1;

