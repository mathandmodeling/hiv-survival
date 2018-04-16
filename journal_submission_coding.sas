/********************************************************************/
/*																	*/
/*																	*/
/*		     Data subsetting and formatting					     	*/
/*																	*/
/*																	*/
/********************************************************************/
*formmatting a graph;
proc import datafile="E:\Research\Manuscript\HIVprogression\DATA_FOR_JOURNAL\02_Merged\Journal_submission_Af_demog_vst_hivcare_cd4_s.xls" dbms=xls out=Af_demog_vst_hivcare_cd4_s replace;
run;
data af_demog_vst_hivcare_cd4_m;
	retain RHSP_ID afb_geo ln_afb_geo Time aflatoxinB1 age date;
	array afb[3] aflatoxin_r1 aflatoxin_r2 aflatoxin_r3;
	array ageyrs[3] age1 age2 age3;
	array date_samp[3] date1 date2 date3;
	afb_geo=.;
	ln_afb_geo=.;
	set af_demog_vst_hivcare_cd4_s;
	
	do i=1 to 3;
	Time=i;
	aflatoxinB1=afb[i];
	age=ageyrs[i];
	date=date_samp[i];
	output; 
	end;
	drop afb_geo ln_afb_geo aflatoxin_r1 aflatoxin_r2 aflatoxin_r3 age1 age2 age3 date1 date2 date3;
run;
data af_demog_vst_hivcare_cd4_m;
	retain RHSP_ID ln_aflatoxinB1;
	set af_demog_vst_hivcare_cd4_m;
	ln_aflatoxinB1=log(aflatoxinB1);*log transformation of aflatoxin B1-Lys;
run;

/*subsetting data*/
*HIV positive data;
data positive_m;
	set af_demog_vst_hivcare_cd4_m;
	if hiv=2 then delete;*HIV=2 (HIV negative);
run;
*HIV negative data;
data negative_m;
	set af_demog_vst_hivcare_cd4_m;
	if hiv=1 then delete; *HIV=1 (HIV positive);
run;

*HIV positive data;
data positive_s;
	retain RHSP_ID afb_geo ln_afb_geo hiv ;
	set af_demog_vst_hivcare_cd4_s;
	if hiv=2 then delete;
run;

/*classify HIV positive participants based on ln(geomean of aflatoxin-lysine), lower than median and higher than median*/
	
proc means data=positive_s median;
var ln_afb_geo;
output out=median_dataset median(ln_afb_geo)=med;
run;

proc sql noprint;
select med into :median_value
from median_dataset;
quit;

%put &median_value;

data positive_s;
set positive_s;
if ln_afb_geo > &median_value then af_class=2;
else af_class=1;
run;


/*creat variables for survival analysis*/
*creat variable for survival analysis (incubation time (incubation), earliest date of HIV confirm (hiv_confirm_d1), date of AIDS confirm(event_date), AIDS progression (event);

data positive_s;
	retain RHSP_ID afb_geo ln_afb_geo aflatoxin_r1 date1 age1 aflatoxin_r2 date2 age2 aflatoxin_r3 date3 age3 incubation event hiv_confirm ;
	format hiv_confirm date1-date3 mmddyy10.;
	set positive_s;

	*to define event (AIDS);
	if cd4_conver_date^=. then event=1;*AIDS (event=1);
	else event=0;*NO-AIDS(event=0);

	
	hiv_confirm=min(first_hiv_confirm, enroldate);*earliest date of hiv confirmation;
			
	*last visit before 01AUG2011 and the date of event (AIDS);
	if cd4_conver_date^=. then conver_date=cd4_conver_date;
	if cd4_conver_date=. then conver_date=last_visit;

	*calculating incubation time;
	incubation=(conver_date-hiv_confirm)/365;
	
	*calculating the age at seroconversion;
	if date3^=. then 
	do;
		diff3=(date3-hiv_confirm)/365;
		age_at_seroconversion=age3-diff3;
	end;
	else if date2^=. then 
	do;
		diff2=(date2-hiv_confirm)/365;
		age_at_seroconversion=age2-diff2;
	end;
	else if date1^=. then 
	do;
		diff1=(date1-hiv_confirm)/365;
		age_at_seroconversion=age1-diff1;
	end;
run;



/********************************************************************/
/*																	*/
/*																	*/
/*	 		Logistic regression (HIV(+) vs. HIV(-))					*/
/*																	*/
/*																	*/
/********************************************************************/


/************Logistic regression (Odds ratio)***********************/
proc logistic data=af_demog_vst_hivcare_cd4_s;*multivariate logistic regression;
	class hiv education(ref='0') occupation(ref='1') religion(ref='1') tribe(ref='1') area(ref='1') sex(ref='1') marital_dis(ref='0')/param=ref;
	model hiv= ln_afb_geo age1 education occupation religion tribe area sex marital_dis/selection=backward slstay=0.1;
run;
proc logistic data=af_demog_vst_hivcare_cd4_s ;*univariate logistic regression;
	class hiv;
	model hiv= ln_afb_geo;
run;

/*********************Case-Control (Table)**************************/
proc freq data=af_demog_vst_hivcare_cd4_s;
	table hiv*sex/chisq;
	table hiv*education/chisq;
	table hiv*occupation/chisq;
	table hiv*marital_dis/chisq;
	table hiv*tribe/chisq;
	table hiv*area/chisq;
run;

proc npar1way data=af_demog_vst_hivcare_cd4_s;
	class hiv;
	var ln_afb_geo age1;
run;

proc print data=af_demog_vst_hivcare_cd4_m;
where hiv=2;
run;

proc freq data=af_demog_vst_hivcare_cd4_m;
table time;
run;

/***********************Box Plot across rounds***********************/

*the effect of time on the level of aflatoxin B1-lysine using GEE in the HIV positives; 
proc genmod data=af_demog_vst_hivcare_cd4_m;
    class RHSP_ID time (ref='1')/param=ref;
    model  ln_aflatoxinB1 = time;
    repeated  subject=RHSP_ID / type=ind corrw;
	contrast "time effect" time 1 -1 0, time 0 1 -1;
	where hiv=1;
run;
proc genmod data=af_demog_vst_hivcare_cd4_m;
    class RHSP_ID time/param=ref;
    model  aflatoxinB1 = time;
    repeated  subject=RHSP_ID / type=ind corrw;
	contrast "time effect" time 1 -1 0, time 0 1 -1;
	where hiv=1;
run;

*the effect of time on the level of aflatoxin B1-lysine using GEE in the HIV negatives; 
proc genmod data=af_demog_vst_hivcare_cd4_m;
    class RHSP_ID time/param=ref;
    model  ln_aflatoxinB1 = time;
    repeated  subject=RHSP_ID / type=ind corrw;
	contrast "time effect" time 1 -1 0 time 0 1 -1;
	where hiv=2;
run;
proc genmod data=af_demog_vst_hivcare_cd4_m;
    class RHSP_ID time/param=ref;
    model  aflatoxinB1 = time;
    repeated  subject=RHSP_ID / type=ind corrw;
	contrast "time effect" time 1 -1 0, time 0 1 -1;
	where hiv=2;
run;
*Drawing the box plots of aflatoxinB1-lysine across times;
*in HIV positive;
proc sort data=positive_m;
	by time;
run;
proc boxplot data =positive_m;
  plot ln_aflatoxinB1*time;
run;

*in HIV negative;
proc sort data=negative_m;
	by time;
run;
proc boxplot data =negative_m;
  plot ln_aflatoxinB1*time;
run;

/********************************************************************/
/*																	*/
/*																	*/
/*	      	HIV positive dataset (Survival Analysis)	     		*/
/*																	*/
/*																	*/
/********************************************************************/
*Multiple Cox Proportional Hazard model;

proc phreg data=positive_s;
      class education(ref='0') occupation(ref='1') religion(ref='1') tribe(ref='1') sex (ref='1')   area(ref='1')     marital_dis(ref='0')/param=ref;
model incubation*event(0)=ln_afb_geo ln_cd41 age_at_seroconversion education occupation religion tribe area sex marital_dis/selection=backward slstay=0.1;
      assess ph/resample seed=385;
run;


proc freq data=positive_s;
	tables af_class*sex/chisq;
	tables af_class*education/chisq;
	tables af_class*marital_dis/chisq;
	tables af_class*tribe/chisq;
	tables af_class*area/chisq;
run;

proc npar1way data=positive_s;
	class af_class;
	var afb_geo ln_afb_geo ln_cd41 age_at_seroconversion;
run;


proc lifetest data=positive_s atrisk plots=survival(atrisk cb) outs=outrural;
	strata af_class;
	time incubation*event(0);
run; 

/****************************************************************************/
/*																			*/
/*																			*/
/*					Testting proportional hazard assumption 				*/
/*  	        	(cumulative sums of martingale residuals)			    */
/*																			*/
/****************************************************************************/

proc phreg data=positive_s;
      class education(ref='0') occupation(ref='1') religion(ref='1') tribe(ref='1') sex (ref='1')   area(ref='1')     marital_dis(ref='0')/param=ref;
model incubation*event(0)=ln_afb_geo ln_cd41 age_at_seroconversion education occupation religion tribe area sex marital_dis/selection=backward slstay=0.1;
      assess ph/resample seed=385;
run;

