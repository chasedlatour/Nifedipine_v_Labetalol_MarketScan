
/********************************************************************************************************************************************
PROGRAM: 00_create_reference_files.sas
PROGRAMMER: Chase Latour
PURPOSE: The purpose of this program is to derive SAS datasets that contain reference code lists.
	
Goal: 
Output data: 

Date: 2/7/2024
********************************************************************************************************************************************/

/*run local:*/
/*SIGNOFF;*/
/*%LET server=n2.schsr.unc.edu 1234; */
/*options comamid=tcp remote=server; */
/*signon username=_prompt_;*/



/*--------------------------------------------------------------------------------------------------------------*/
/*		00 			SETUP 																						*/
/*--------------------------------------------------------------------------------------------------------------*/
*Run setup macro and define libnames;

options sasautos=(SASAUTOS "/local/projects/marketscan_preg/Latour_23_2322/programs/macros");
/*change "saveLog=" to "Y" when program is closer to complete*/
%setup(sample=random1pct, programname=00_create reference files, savelog=N);

*Set this up for the GEMs Mapping;
libname projlib "/local/projects/marketscan_preg/Latour_23_2322/data/coderef";

*Set up this library for the codelist;
libname codelist xlsx '/local/projects/marketscan_preg/Latour_23_2322/data/coderef/Codelist.xlsx';

options mprint;


******************************************************************************************************************************************;
/*Create local mirros of the libraries from the set up macro - Run locally*/
/*libname raw slibref=raw server = server;*/
/*libname der slibref=der server = server;*/
/*libname lout slibref=out server=server; */
/*libname lana slibref=ana server=server;*/
/*libname lpreg slibref=preg server=server;*/
/*libname lwork slibref=work server=server;*/
/*libname lexpref slibref=expref server=server;*/
/*libname lcovref slibref=covref server=server; */
/*libname loutref slibref=outref server=server;*/
/*libname lrxcov slibref=rxcov server=server;*/
/*libname latc slibref=atc server=server;*/
/*libname lprojlib slibref=projlib server=server;*/

/*%inc "/local/projects/marketscan_preg/Latour_23-2322/programs/FormatStatements_CDWH.sas";*/


*Macro for simple situation where only ICD-10 codes available in the Excel file - no manipulation required;

%macro simple_dwnld(covar);

	data &covar;
	format code $8.;
	length code $8; /*Manually specify the code variable to be a character variable of length 8*/
	set codelist."&covar"n;
	run;

	data covref.&covar._dx;
	set &covar;
		where include = 1;
		code = compress(code, '.');
	run;

%mend;


*Macro for simple NDC grab when Excel file only has ATC codes;

%macro simple_med(covar);

	data &covar;
	format code $8.;
	length code $8;
	set codelist."&covar"n;
	run;

	*Convert the ATCs to NDCs;
	proc sql;
		create table rxcov.&covar._rx as
		select distinct a.*
		from atc.atc_ndc as a
		inner join (select * from &covar where include = 1) as b
	on a.atc = b.code
	;
	quit;

%mend;




********CHRONIC HYPERTENSION code list;

%simple_dwnld(chronic_htn);


%*******SECONDARY HYPERTENSION code list;

%simple_dwnld(second_htn);


********GESTATIONAL HYPERTENSION code list;

%simple_dwnld(gest_htn);

********PREECLAMPSIA code list;

%simple_dwnld(preeclampsia);

********UNSPECIFIED HYPERTENSION code list;

%simple_dwnld(unsp_htn);

********HELLP SYNDROME code list;

%simple_dwnld(hellp);

********ECLAMPSIA code list;

%simple_dwnld(eclampsia);



********OPIOID USE DISORDER code list;

data oud;
format code $8.;
length code $8;
set codelist."opioid_use_disorder"n;
run;

*Get a list of all diagnosis and procedure codes;
data covref.oud_dx;
set oud;
	where codetype ne "ATC" and include = 1;
	code = compress(code, '.');
run;

*Now get all of the NDC codes for the ATC codes.;
proc sql;
	create table rxcov.oud_rx as
	select distinct a.*
	from atc.atc_ndc as a
	inner join (select * from oud where codetype = "ATC" and include = 1) as b
	on a.atc = b.code
	;
	quit;

/**Do a quick check to make sure it looks as expected;*/
/*proc freq data=rxcov.oud_rx;*/
/*	table atc_label drug_name;*/
/*run;*/




************NAUSEA-VOMITING code list;

data nausea_vomiting;
format code $8.;
length code $8;
set codelist."nausea_vomiting"n;
run;

*Get a list of all diagnosis and procedure codes;
data covref.nausea_vomiting_dx;
set nausea_vomiting;
	where codetype ne "ATC" and include = 1;
	code = compress(code, '.');
run;

*Now get all of the NDC codes for the ATC codes.;
proc sql;
	create table rxcov.nausea_vomiting_rx as
	select distinct a.*
	from atc.atc_ndc as a
	inner join (select * from nausea_vomiting where codetype = "ATC" and include = 1) as b
	on a.atc = b.code
	;
	quit;


***********DIABETES code list;

%simple_dwnld(diabetes);


************ANTIDIABETICS for type-1 and type-2 diabetes;

%simple_med(t1t2dm_antidiabetics);



************ANTIDIABETICS for type-2 diabetes only;

data t2dm_rx;
	format group $25.;
	length code $25; 
	set codelist."t2dm_antidiabetics"n;
run;

*Select all the distint group values from t2dm_rx;
proc sql;
	create table t2dm_group as
	select distinct upcase(group) as group
	from t2dm_rx
	;
	quit;


proc sql;
    create table t2dm_rx2 as 
    select * 
    from red.redbook r
    where exists (
        select 1 
        from t2dm_group t
        where upcase(r.gennme) like cats('%', t.group, '%')
    );
	quit;
	
	
*Create the variables that we will need to match those from the First Data Bank file;
data t2dm_rx3;
set t2dm_rx2;
	format ndc9 $9. atc_label $56.;
	ndc9 = substr(ndcnum, 1, 9);
	atc_label = upcase(gennme);
run;	

data rxcov.t2dm_antidiabetics_rx;
set t2dm_rx3;
	where prodnme not in ("WEGOVY" "SAXENDA");
run;


********OBESITY medications;

*First get the GLP-1 agonists;
data rxcov.wgtloss_glp1_rx;
set t2dm_rx3;
	where prodnme in ("WEGOVY" "SAXENDA");
run;

*Now get the rest of the weightloss drugs. Identify by generic name and product name;
data others;
	format group $25.;
	set codelist."obesity_meds"n;
	where type="Other";
run;

*Select all the distinct description and product values from others;
proc sql;
	create table generic as
	select distinct upcase(description) as generic 
	from others
	;
	quit;
	
proc sql;
	create table product as
	select distinct upcase(product) as product 
	from others
	having product ne "NA"
	;
	quit;

*Grab all redbook rows based on generic and product then stack;
proc sql;
    create table generic2 as 
    select * 
    from red.redbook r
    where exists (
        select 1 
        from generic t
        where upcase(r.gennme) like cats('%', t.generic, '%')
    );
    
    create table product2 as
    select * 
    from red.redbook r
    where exists (
        select 1 
        from product t
        where upcase(r.prodnme) like cats('%', t.product, '%')
    );
	quit;
	
*Now stack them using union so no duplicate rows;
proc sql;
	create table both as
	select * from generic2
	union
	select * from product2;
	quit;
	
data rxcov.wgtloss_other_rx;
set both;
	format ndc9 $9. atc_label $56.;

	ndc9 = substr(ndcnum, 1, 9);
	atc_label = upcase(gennme);

	if prodnme in ("ADIPEX" "ADIPEX-P") then output;
	if prodnme in ("BELVIQ" "BELVIQ XR") then output;
	if prodnme in ("BONTRIL" "BONTRIL PDM" "BONTRIL SLOW-RELEASE") then output;
	if upcase(prodnme) = "CONTRAVE" or
		upcase(prodnme) = "DIDREX" or
		upcase(prodnme) = "LOMAIRA" or
		upcase(prodnme) = "MERIDIA" or
		upcase(prodnme) = "QSYMIA" or
		upcase(prodnme) = "REGIMEX" or
		upcase(prodnme) = "SUPRENZA" or 
		upcase(prodnme) = "XENICAL" then output
	;
	if gennme in ("Diethylpropion Hydrochloride" "Diethylpropion Hydrochloride/Tartaric Acid") then output;
run;
	
	
	
	
********BARIATRIC surgery;

*Primary surgery;
%simple_dwnld(bariatric_surgery);

*Revisions;
%simple_dwnld(rev_bariatric);




********METFORMIN code list;

%simple_med(metformin);


*********SMOKING and tobacco;

data smoking;
format code $8.;
length code $8;
set codelist."smoking"n;
run;

*Get the medications to identify smoking;
proc sql;
	create table rxcov.smoking_rx as
	select distinct a.*
	from atc.atc_ndc as a
	inner join (select * from smoking where include = 1 and codetype = "ATC") as b
	on a.atc = b.code
	;
	quit;

*Identify the procedure codes;
data smoking_proc smoking_dx9;
set smoking;
	where include = 1;
	code = compress(code, '.');
	if codetype in ('HCPCS' 'CPT') then output smoking_proc;
		else if codetype = 'DX9' then output smoking_dx9;
run;

/**GEMs for the ICD-9 diagnnosis codes;*/
/*proc sql;*/
/*	create table icd9_dx_codelist as*/
/*	select distinct a.source as icd9_dx_self*/
/*	from projlib.icd9to10dx as a*/
/*	inner join smoking_dx as b*/
/*	on a.source = compress(b.code, '.')*/
/*	;*/
/*	quit;*/
/**/
/*%gemsmap(smk, dx);*/
/**This gave incomplete list, so manually reviewed and added codes manually.;*/

data smoking_dx10;
format code $8.;
length code $8;
set codelist."smoking_dx10"n;
run;

data covref.smoking_dx;
format code $8.;
set smoking_proc (keep = code include codetype description) smoking_dx10;
	where include = 1;
	code = compress(code, '.');
run;





*********ALCOHOL code list;

data alcohol;
format code $8.;
length code $8;
set codelist."alcohol"n;
run;

*Output the DX10 and DX9 codes. DX9 codes need to be forward-backward mapped;
data alcohol_dx10 alcohol_dx9 problem;
set alcohol;
	code = compress(code, '.');
	if codetype = "DX10" then output alcohol_dx10;
		else if codetype = "DX9" then output alcohol_dx9;
		else output problem;
run;


*FBM the DX-9 codes;
proc sql;
	create table icd9_dx_codelist as
	select distinct a.source as icd9_dx_self
	from projlib.icd9to10dx as a
	inner join alcohol_dx9 as b
	on a.source = b.code
	;
	quit;

%gemsmap(alc, dx);
*Manually reviewed this list and decided that a few additional codes should have been added;

data alcohol_dx10_manual;
input code $8.;
cards;
	F10130
	F10131
	F10132
	F10139
	F1090
	F10930
	F10931
	F10932
	F10939
	;
run;

*Stack the DX10 datasets together and output;
proc sql;
	create table covref.alcohol_dx as
	select min(include) as include, min(codetype) as codetype, code,
		max(Description) as Description
	from (
		select 1 as include, "DX10" as codetype, code, "" as Description
		from alcohol_dx10_manual
		union corr
		select 1 as include, "DX10" as codetype, icd10dx as code, "" as Description
		from alc_dx_fbm_final
		union corr
		select *
		from alcohol_dx10
	)
	group by code /*Make sure that we only have unique codes*/
	;
	quit;

proc datasets gennum = all;
	delete alc: smok: smk_:;
run; quit; run;




*********MIGRAINE code list;

data migraine_dx;
format code $8.;
length code $8;
set codelist."migraine"n;
run;

*Output Dx10 and ATC codes separately;
data covref.migraine_dx migraine_atc problem;
set migraine_dx;
	where include = 1;
	code = compress(code, '.');
	if codetype = "DX10" then output covref.migraine_dx;
		else if codetype = "ATC" then output migraine_atc;
		else output problem;
run;

*Now, grab all the NDC codes where the ATC code contains N02CC;
data migraine_rx;
    set atc.atc_ndc;
    if index(atc, 'N02CC') > 0;
run;

*Make sure all rows are distinct;
proc sql;
	create table rxcov.migraine_rx as
	select distinct *
	from migraine_rx
	;
	quit;

/*proc freq data=migraine_rx;*/
/*	table atc_label  / missing;*/
/*run;*/





**********RECURRENT PREGNANCY LOSS code list;

%simple_dwnld(recur_preg_loss);




**********CHRONIC KIDNEY DISEASE code list;

data ckd;
format code $8.;
length code $8;
set codelist."ckd"n;
run;

*Output Dx10 and Dx 9 codes separately;
data ckd_dx9 ckd_dx10 problem;
set ckd;
	where include = 1;
	code = compress(code, '.');
	if codetype = "DX9" then output ckd_dx9;
		else if codetype = "DX10" then output ckd_dx10;
		else output problem;
run;

*Forward-backwards GEMSMAP the DX9 codes;
proc sql;
	create table icd9_dx_codelist as
	select distinct a.source as icd9_dx_self
	from projlib.icd9to10dx as a
	inner join ckd_dx9 as b
	on a.source = b.code
	;
	quit;

%gemsmap(ckd, dx);

*Now stack the datasets;
proc sql;
	create table covref.ckd_dx as
	select min(include) as include, min(codetype) as codetype, code,
		max(Description) as Description
	from (
		select distinct include, codetype, code, description
		from ckd_dx10
		union corr
		select distinct 1 as include, "DX10" as codetype, icd10dx as code, "" as description
		from ckd_dx_fbm_final
		)
	group by code
	;
	quit;


proc datasets gennum = all;
	delete ckd:;
run; quit; run;







*******OTHER SUBSTANCE USE DISORDER;

data sud;
format code $8.;
length code $8;
set codelist."non_opioid_sud"n;
run;

*Forward-backwards GEMSMAP the DX9 codes;
proc sql;
	create table icd9_dx_codelist as
	select distinct a.source as icd9_dx_self
	from projlib.icd9to10dx as a
	inner join (select * from sud where include = 1) as b
	on a.source = compress(b.code, '.')
	;
	quit;

%gemsmap(sud, dx);

proc sql;
	create table covref.other_sud_dx as
	select distinct 1 as include, "DX10" as codetype, icd10dx as code, "" as description
	from sud_dx_fbm_final
	;
	quit;

proc datasets gennum = all;
	delete sud:;
run; quit; run;





*******OBESITY;

%simple_dwnld(obese);

*******STROKE - TIA;

%simple_dwnld(stroke_tia);
	


*******MYOCARDIAL INFARCTION;

%simple_dwnld(mi);





*******DIABETIC RETINOPATHY;

data retinopathy;
format code $8.;
length code $8;
set codelist."retinopathy"n;
run;

*Forward-backwards GEMSMAP the DX9 codes;
proc sql;
	create table icd9_dx_codelist as
	select distinct a.source as icd9_dx_self
	from projlib.icd9to10dx as a
	inner join (select * from retinopathy where include = 1) as b
	on a.source = compress(b.code, '.')
	;
	quit;

%gemsmap(ret, dx);
*Spot checking looked good;

*Create final dataset;
proc sql;
	create table covref.retinopathy_dx as
	select distinct 1 as include, "DX10" as codetype, icd10dx as code, "" as description
	from ret_dx_fbm_final
	;
	quit;

proc datasets gennum = all;
	delete ret:;
run;





*******ANTIPHOSPHOLIPID;

%simple_dwnld(antiphospholipid);


*******LUPUS;

data lupus;
format code $8.;
length code $8;
set codelist."lupus"n;
run;

*Forward-backwards GEMSMAP the DX9 codes;
proc sql;
	create table icd9_dx_codelist as
	select distinct a.source as icd9_dx_self
	from projlib.icd9to10dx as a
	inner join (select * from lupus where include = 1) as b
	on a.source = compress(b.code, '.')
	;
	quit;

%gemsmap(lup, dx);
*Spot checking looked good;

*Create final dataset;
proc sql;
	create table covref.lupus_dx as
	select distinct 1 as include, "DX10" as codetype, icd10dx as code, "" as description
	from lup_dx_fbm_final
	;
	quit;

proc datasets gennum = all;
	delete lup:;
run;






*******ASTHMA;

data asthma;
format code $8.;
length code $8;
set codelist."asthma"n;
run;

*Output diagnoses and atc separately;
data asthma_dx asthma_atc problem;
set asthma;
	where include = 1;
	code = compress(code, '.');

	if codetype = "DX9" then output asthma_dx;
		else if codetype = "ATC" then output asthma_atc;
		else output problem;
run;


*Forward-backwards GEMSMAP the DX9 codes;
proc sql;
	create table icd9_dx_codelist as
	select distinct a.source as icd9_dx_self
	from projlib.icd9to10dx as a
	inner join asthma_dx as b
	on a.source = b.code
	;
	quit;

%gemsmap(ast, dx);
*Spot checking looked good;

*Create final dataset;
proc sql;
	create table covref.asthma_dx as
	select distinct 1 as include, "DX10" as codetype, icd10dx as code, "" as description
	from ast_dx_fbm_final
	;
	quit;

*Finally, output the relevant NDCs;

proc sql;
	create table rxcov.asthma_rx as
	select distinct a.*
	from atc.atc_ndc as a
	inner join asthma_atc as b
	on a.atc = b.code
	;
	quit;


proc datasets gennum = all;
	delete ast:;
run;




***********CORONARY ARTERY DISEASE;

%simple_dwnld(cad);




***********HEART VALVE DISEASE;

data heart_valve_disease;
format code $8.;
length code $8;
set codelist."heart_valve_disease"n;
run;

*Forward-backwards GEMSMAP the DX9 codes;
proc sql;
	create table icd9_dx_codelist as
	select distinct a.source as icd9_dx_self
	from projlib.icd9to10dx as a
	inner join (select * from heart_valve_disease where include = 1) as b
	on a.source = compress(b.code, '.')
	;
	quit;

%gemsmap(hea, dx);

*Create final dataset;
proc sql;
	create table covref.heart_valve_dx as
	select distinct 1 as include, "DX10" as codetype, icd10dx as code, "" as description
	from hea_dx_fbm_final
	;
	quit;

proc datasets gennum = all;
	delete hea:;
run;






**************CARDIOMYOPATHY;

%simple_dwnld(cardiomyopathy);



**************HEART FAILURE;

data hf;
format code $8.;
length code $8;
set codelist."heartfailure"n;
run;

*Forward-backwards GEMSMAP the DX9 codes;
proc sql;
	create table icd9_dx_codelist as
	select distinct a.source as icd9_dx_self
	from projlib.icd9to10dx as a
	inner join (select * from hf where include = 1) as b
	on a.source = compress(b.code, '.')
	;
	quit;

%gemsmap(hf, dx);

*Create final dataset;
proc sql;
	create table covref.hf_dx as
	select distinct 1 as include, "DX10" as codetype, icd10dx as code, "" as description
	from hf_dx_fbm_final
	;
	quit;

proc datasets gennum = all;
	delete hf:;
run;





**************OTHER HEART DISEASE;

%simple_dwnld(other_heart_disease);




**************ARRHYTHMIA;

data arrhythmia;
format code $8.;
length code $8;
set codelist."arrhythmia"n;
run;

*Forward-backwards GEMSMAP the DX9 codes;
proc sql;
	create table icd9_dx_codelist as
	select distinct a.source as icd9_dx_self
	from projlib.icd9to10dx as a
	inner join (select * from arrhythmia where include = 1) as b
	on a.source = compress(b.code, '.')
	;
	quit;

%gemsmap(arr, dx);

*Create final dataset;
proc sql;
	create table covref.arrhythmia_dx as
	select distinct 1 as include, "DX10" as codetype, icd10dx as code, "" as description
	from arr_dx_fbm_final
	;
	quit;

proc datasets gennum = all;
	delete arr:;
run;






**************ENDOCARDITIS;

%simple_dwnld(endocarditis);


**************MYOCARDITIS PERICARDITIS;

%simple_dwnld(myo_pericarditis);

**************CONGENITAL HEART DEFECT;

%simple_dwnld(congenital_heart_defect);




**************HYPERTHYROID;

data hyperthyroid;
format code $8.;
length code $8;
set codelist."hyperthyroid"n;
run;

*Now, output the diagnosis codes and atc codes separately;
data covref.hyperthyroid_dx hyper_atc problem;
set hyperthyroid;
	where include = 1;
	code = compress(code, '.');
	if codetype = "DX10" then output covref.hyperthyroid_dx;
		else if codetype = "ATC" then output hyper_atc;
		else output problem;
run;

*Now get the relevant ATC codes;
proc sql;
	create table rxcov.hyperthyroid_rx as
	select distinct a.*
	from atc.atc_ndc as a
	inner join hyper_atc as b
	on a.atc = b.code
	;
	quit;



**************HYPOTHYROID;

data hypothyroid;
format code $8.;
length code $8;
set codelist."hypothyroid"n;
run;

*Now, output the diagnosis codes and atc codes separately;
data covref.hypothyroid_dx hypo_atc problem;
set hypothyroid;
	where include = 1;
	code = compress(code, '.');
	if codetype = "DX10" then output covref.hypothyroid_dx;
		else if codetype = "ATC" then output hypo_atc;
		else output problem;
run;

*Now get the relevant ATC codes;
proc sql;
	create table rxcov.hypothyroid_rx as
	select distinct a.*
	from atc.atc_ndc as a
	inner join hypo_atc as b
	on a.atc = b.code
	;
	quit;





**************DEPRESSION;

%simple_dwnld(depression);


**************ANXIETY;

data anxiety;
format code $8.;
length code $8;
set codelist."anxiety"n;
run;

*Now, output the diagnosis codes and atc codes separately;
data covref.anxiety_dx anxiety_atc problem;
set anxiety;
	where include = 1;
	code = compress(code, '.');
	if codetype = "DX10" then output covref.anxiety_dx;
		else if codetype = "ATC" then output anxiety_atc;
		else output problem;
run;

*Now get the relevant ATC codes;
proc sql;
	create table rxcov.anxiety_rx as
	select distinct a.*
	from atc.atc_ndc as a
	inner join anxiety_atc as b
	on a.atc = b.code
	;
	quit;




**************BIPOLAR DISORDER;

%simple_dwnld(bipolar_disorder);



**************MOOD STABILIZERS;

%simple_med(mood_stabilizers);


**************SCHIZOPHRENIA;

%simple_dwnld(schizophrenia);


**************ANTIPSYCHOTICS;

%simple_med(antipsychotics);




**************ADHD;

data adhd;
format code $8.;
length code $8;
set codelist."adhd"n;
run;

*Forward-backwards GEMSMAP the DX9 codes;
proc sql;
	create table icd9_dx_codelist as
	select distinct a.source as icd9_dx_self
	from projlib.icd9to10dx as a
	inner join (select * from adhd where include = 1) as b
	on a.source = compress(b.code, '.')
	;
	quit;

%gemsmap(adh, dx);

*Create final dataset;
proc sql;
	create table covref.adhd_dx as
	select distinct 1 as include, "DX10" as codetype, icd10dx as code, "" as description
	from adh_dx_fbm_final
	;
	quit;

proc datasets gennum = all;
	delete adh:;
run;





**************ADHD MEDICATIONS;

%simple_med(adhd_meds);

**************PTSD;

%simple_dwnld(ptsd);


**************PTSD MEDICATIONS;

%simple_med(ptsd_meds);

**************ANGINA;

%simple_dwnld(angina);

**************ATHEROSCLEROSIS;

%simple_dwnld(atherosclerosis);


**************PERIPHERAL VASCULAR DISEASE;

%simple_dwnld(peripheral_vasc_disease);

**************BREAST CANCER;

%simple_dwnld(breast_cancer);

**************COLORECTAL CANCER;

%simple_dwnld(colorectal_cancer);

**************ENDOMETRIAL CANCER;

%simple_dwnld(endometrial_cancer);

**************LUNG CANCER;

%simple_dwnld(lung_cancer);

**************UROLOGIC CANCER;

%simple_dwnld(urologic_cancer);

**************HYPERLIPIDEMIA;

%simple_dwnld(hyperlipidemia);

**************ANEMIA;

%simple_dwnld(anemia);

**************TERATOGENIC MEDICATIONS;

%simple_med(teratogenic_meds);

**************STATINS;

%simple_med(statins);

**************BENZODIAZEPINES;

%simple_med(benzodiazepines);

**************ANTIDEPRESSANTS;

%simple_med(antidepressants);

**************ANTICONVULSANTS;

%simple_med(anticonvulsants);

**************ANTIHYPERTENSIVES;

%simple_med(antihypertensives);





**************SICKLE CELL DISEASE;

data sickle_disease;
format code $8.;
length code $8;
set codelist."sickle_disease"n;
run;

*Forward-backwards GEMSMAP the DX9 codes;
proc sql;
	create table icd9_dx_codelist as
	select distinct a.source as icd9_dx_self
	from projlib.icd9to10dx as a
	inner join (select * from sickle_disease where include = 1) as b
	on a.source = compress(b.code, '.')
	;
	quit;

%gemsmap(sic, dx);

*Create final dataset;
proc sql;
	create table covref.sickle_disease_dx as
	select distinct 1 as include, "DX10" as codetype, icd10dx as code, "" as description
	from sic_dx_fbm_final
	;
	quit;

proc datasets gennum = all;
	delete sic:;
run;




**************SICKLE CELL TRAIT;

data sickle_trait;
format code $8.;
length code $8;
set codelist."sickle_trait"n;
run;


*Forward-backwards GEMSMAP the DX9 codes;
proc sql;
	create table icd9_dx_codelist as
	select distinct a.source as icd9_dx_self
	from projlib.icd9to10dx as a
	inner join (select * from sickle_trait where include = 1) as b
	on a.source = compress(b.code, '.')
	;
	quit;

%gemsmap(sic, dx);

*Create final dataset;
proc sql;
	create table covref.sickle_trait_dx as
	select distinct 1 as include, "DX10" as codetype, icd10dx as code, "" as description
	from sic_dx_fbm_final
	;
	quit;

proc datasets gennum = all;
	delete sic:;
run;





****Finally, we want to identify relevant codes for our exposures;

*******LABETALOL;

*ATC for only labetalol is C07AG01;

data labetalol;
set atc.atc_ndc;
	where atc = 'C07AG01';
run;

proc freq data=labetalol;
	table route;
run;

*We want to limit to only the oral route of labetalol;
data expref.labetalol;
set labetalol;
	where route = "oral";
run;

/**Check it out;*/
/*proc freq data=labetalol;*/
/*	table drug_name;*/
/*run;*/


*******NIFEDIPINE;

data nifedipine;
set atc.atc_ndc;
	where atc = 'C08CA05';
run;

/*proc freq data=nifedipine;*/
/*	table route form;*/
/*run;*/

*Now, limit to oral forms that are modified - extended - release nifedipine;
data expref.nifedipine;
set nifedipine;
	where route = "oral" and index(form, "extended release");
run;

/*proc freq data=expref.nifedipine;*/
/*	table form;*/
/*run;*/

********NON-EXPOSURE ANTIHYPERTENSIVES;

*We want to identify all antihypertensives that are not included in the exposures;

/*proc contents data=rxcov.antihypertensives_rx; run;*/
proc sql;
	create table antihypertensives as
	select a.*, not missing(b.ndc11) as labetalol, not missing(c.ndc11) as nifedipine
	from rxcov.antihypertensives_rx as a
	left join expref.labetalol as b
	on a.ndc11 = b.ndc11
	left join expref.nifedipine as c
	on a.ndc11 = c.ndc11
	;
	quit;

/**Confirm that looks as expected;*/
/*proc freq data=antihypertensives (where = (labetalol = 1));*/
/*	table atc_label;*/
/*run;*/
/*proc freq data=antihypertensives (where = (nifedipine = 1));*/
/*	table atc_label;*/
/*run;*/

data expref.other_antihypertensives;
set antihypertensives;
	where labetalol = 0 and nifedipine = 0;
run;

/*data others labetalol nifedipine problem;*/
/*set antihypertensives;*/
/*	if labetalol = 0 and nifedipine = 0 then output others;*/
/*		else if labetalol = 1 then output labetalol;*/
/*		else if nifedipine = 1 then output nifedipine;*/
/*		else output problem;*/
/*run;*/












	














	




