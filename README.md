# Longitudinal-Data-Analysis

## Introduction 

The project  aims to investigate the relationship between nitric oxide (NO) concentration and flow of exhalation. The project analyzes a dataset containing clustered and longitudinal data, specifically focusing on the differences between an asthma group and a healthy group.

## Data Analysis

The initial descriptive analysis reveals that the missing observations in the dataset are due to attrition or dropout. The variability in the NO concentration is influenced by both individual differences and growth rates.

To model the concentration of NO over time, a linear mixed model with fixed effects for flow, diagnosis, and their interaction is fitted. The model selection process using different covariance structures reveals that the model with unstructured covariance provides the best fit based on the AIC criterion. 

Further, separate random effects models are fitted for the asthma group and the healthy group. The models estimate the variances of the random intercept and slope for each group. Confidence intervals are calculated for these parameters, providing insights into the variability in baseline NO concentration and flow change for both groups.

Finally, a combined random effects model is fitted for both groups, allowing for a different covariance structure between them. Likelihood ratio tests confirm that the covariance structure differs significantly between the two groups. 

## Discussion 

The plot suggests the presence of random intercepts and slopes, indicating variations in NO concentration with increasing flow of exhalation. The asthma group exhibits higher variability compared to the healthy group. While both groups show a decrease in NO concentration as the flow of exhalation increases, it is not clear from the plot whether one group has higher NO concentration than the other. The fixed effects analysis shows that the impact of diagnosis (asthma or healthy) is not statistically significant, but both flow and the interaction term are significant. Further the result reveals a significant interaction effect between flow and diagnosis, suggesting that the average NO concentration differs between the two groups at different levels of flow.


## Conclusion

In conclusion, the project provides a comprehensive analysis of the relationship between NO concentration and flow of exhalation using clustered and longitudinal data. The findings highlight the importance of considering random intercepts and slopes, as well as the differences between the asthma and healthy groups, in understanding the dynamics of NO concentration.

## Project Code_SAS

```

project2.sas

LIBNAME ACLD 'H:/2020-21/ACLD/Project2';
data ACLD.new;
set ACLD.no;
cflow=flow;
run;
proc print data=ACLD.new;
run;
*Separate dataset for the healthy and asthma group;
DATA OnlyHealthy;
SET ACLD.new;
IF (diagnose=1) then
delete;
run;
DATA OnlyAsthma;
SET ACLD.new;
IF (diagnose=0) then
delete;
run;
proc print DATA=OnlyHealthy;
run;
proc print DATA=OnlyAsthma;
run;
title 'Profile Plot of concentration of NO';
PROC SGPLOT DATA=ACLD.new NOAUTOLEGEND;
format diagnose diagnosefmt.;
SERIES X=flow Y=no / GROUP=num GROUPLC=diagnose NAME='diagnose' BREAK
LINEATTRS=(PATTERN=1);
LABEL no='concentration of NO' flow='flow of exhalation';
KEYLEGEND 'diagnose' / NOBORDER TYPE=linecolor LOCATION=inside;
RUN;
title;
*profile Plot by group;
title 'Paneled Profile Plot by group';
proc sgpanel data=ACLD.new noautolegend;
panelby diagnose / columns=2;
series X=flow Y=no / group=num break lineattrs=(pattern=1);
LABEL pi='concentration of NO' time='flow of exhalation';
run;
title;
*smoother plot;
title 'Smoother Plot of concentration of NO ';
PROC SGPLOT DATA=ACLD.new;
LOESS X=flow Y=no / GROUP=diagnose NAME='diagnose';
LABEL no='concentration of NO' flow='flow of exhalation';
KEYLEGEND 'diagnose' / NOBORDER TYPE=linecolor LOCATION=inside;
RUN;

title;
*checking missing data;
proc sgpanel data=ACLD.new noautolegend;
panelby num / columns=7 rows=6;
series x=flow y=no / group=num break lineattrs=(pattern=1);
label day="flow of exhalation" y="concentration of NO";
run;
***********
Model building covariance structure
*********;
**********************************************************************
*unstructured;
PROC MIXED DATA=ACLD.new;
CLASS num;
MODEL no=diagnose flow flow*diagnose /s chisq;
REPEATED / TYPE=UN SUBJECT=num;
ODS OUTPUT Fitstatistics=fit_un(rename=(value=UN));
ODS OUTPUT Dimensions=Parm_UN(rename=(value=Num_UN));
RUN;
* Compund symmetry;
PROC MIXED DATA=ACLD.new;
CLASS num;
MODEL no=diagnose flow flow*diagnose /s chisq;
REPEATED / TYPE=CS SUBJECT=num;
ODS OUTPUT Fitstatistics=fit_CS(rename=(value=CS));
ODS OUTPUT Dimensions=Parm_CS(rename=(value=Num_CS));
RUN;
* Heterogeneous Compound Symmetry;
PROC MIXED DATA=ACLD.new;
CLASS num;
MODEL no=diagnose flow flow*diagnose /s chisq;
REPEATED / TYPE=CSH SUBJECT=num;
ODS OUTPUT Fitstatistics=fit_CSH(rename=(value=CSH));
ODS OUTPUT Dimensions=Parm_CSH(rename=(value=Num_CSH));
RUN;
* combine all loglikelihood, AIC and BIC between models and compare;
Data compare;
merge Fit_CS Fit_CSH Fit_UN;
run;
title 'Output Summary comparing the models';
proc print data=compare label noobs;
run;
title;
*From the AIC's UN does best compared to the others;
* Just a review of the number of parameter estimate per model;
Data all_parmnr;
merge Parm_UN Parm_CS Parm_CSH;

title "Number of parameter estimate per model";
proc print data=all_parmnr (obs=1) label noobs;
run;
title;
/* Now we construct LRT tests for nested models to see what
model does better; */
*a. Unstructured (full) vs compound symmetry (reduced);
data UN_CS;
merge Fit_UN Fit_CS Parm_UN Parm_CS;
if _n_=1 then
do;
Chi_UN_CS=CS - UN;
df_UN_CS=Num_UN - Num_CS;
p_UN_CS=1-probchi(Chi_UN_CS, df_UN_CS);
output;
stop;
end;
run;
title 'Likelihood Ratio Test: Unstructured vs Compound symmetry';
proc print data=UN_CS label noobs;
var UN CS Chi_UN_CS df_UN_CS p_UN_CS;
label UN="-2 loglik UN" CS="-2 loglik CS" Chi_UN_CS="Chi-Square" df_UN_CS="DF"
p_UN_CS="Pr > ChiSq";
run;
title;
*b. Unstructured (full) vs heterogeneous compound symmetry (reduced);
data UN_CSH;
merge Fit_UN Fit_CSH Parm_UN Parm_CSH;
if _n_=1 then
do;
Chi_UN_CSH=CSH - UN;
df_UN_CSH=Num_UN - Num_CSH;
p_UN_CSH=1-probchi(Chi_UN_CSH, df_UN_CSH);
output;
stop;
end;
run;
title 'Likelihood Ratio Test: Unstructured vs Hetero Compound symmetry';
proc print data=UN_CSH label noobs;
var UN CSH Chi_UN_CSH df_UN_CSH p_UN_CSH;
label UN="-2 loglik UN" CSH="-2 loglik CSH" Chi_UN_CSH="Chi-Square"
df_UN_CSH="DF" p_UN_CSH="Pr > ChiSq";
run;
title;
**********************************************************************
combine all loglikelihood, AIC and BIC between models and compare;

Data compare;
merge Fit_CS Fit_CSH Fit_UN;
run;
title 'Table 1. Output Summary comparing the models';
proc print data=compare label noobs;
run;
title;
* testing if covariance structure for treatment groups differ or not ;
proc mixed data=ACLD.new;
class num diagnose ;
MODEL no=diagnose flow flow*diagnose /s chisq;
random Intercept flow / SUBJECT=num G GCORR TYPE=un V VCORR GROUP=diagnose;
repeated / type=simple subject=num;
ods output Fitstatistics=fit_f;
run;
*reduced;
proc mixed data=ACLD.new;
class num diagnose ;
MODEL no=diagnose flow flow*diagnose /s chisq;
random Intercept flow / SUBJECT=num G GCORR TYPE=un V VCORR;
repeated / type=simple subject=num;
ods output Fitstatistics=fit_r;
run;
*LR test use REML, Ho: no different covariance structure for treatment groups ;
data lrt;
set fit_f(obs=1);
if Descr='-2 Res Log Likelihood' then
lrf=Value;
set fit_r(obs=1);
if Descr='-2 Res Log Likelihood' then
lrr=Value;
lr=lrr-lrf;
df=3;
pvalue=1-probchi(lr, df);
keep lr df pvalue;
run;
proc print data=lrt;
title 'LR Test';
run;
*LR test Ho: H0: random intercepts only vs
H1: random intercepts and random slopes;
proc mixed data=ACLD.new;
class num diagnose;
model no=diagnose flow flow*diagnose/s;
random Intercept/ SUBJECT = num G GCORR TYPE= un V VCORR GROUP = diagnose;
repeated / type=simple subject=num;
ods output Fitstatistics=fit1(rename=(value=Random_int));
run;
proc mixed data=ACLD.new;
class num diagnose;
model no=diagnose flow flow*diagnose/s;

random Intercept flow/ SUBJECT = num G GCORR TYPE= un V VCORR GROUP = diagnose;
repeated / type=simple subject=num;
ods output Fitstatistics=fit2(rename=(value=Random_intslope));
run;
data lrt;
set fit1(obs=1);
if Descr='-2 Res Log Likelihood' then lrr=Random_int;
set fit2(obs=1);
if Descr='-2 Res Log Likelihood' then lrf=Random_intslope;
lr=lrr-lrf;
df1=3;
df2=7;
df=df2-df1;
p1=1-probchi(lr,df1);
p2=1-probchi(lr,df2);
p=0.5*p1+0.5*p2;;
keep lr df p;
run;
proc print data=lrt;
title 'LR Test';
run;
*Final Model;
proc mixed data=ACLD.new;
class num diagnose;
model no=diagnose flow flow*diagnose/s;
random Intercept flow/ SUBJECT = num G GCORR TYPE= un V VCORR GROUP = diagnose;
repeated / type=simple subject=num;
run;
*Random effect model for healthy group;
proc mixed data=OnlyHealthy;
class num;
model no=flow / s;
random intercept flow / type=un subject=num g gcorr v=2 vcorr=2;
repeated / type=simple subject=num r=2;
run;
*Random effect model for Asthma group;
proc mixed data=OnlyAsthma;
class num;
model no=flow / s;
random intercept flow / type=un subject=num g gcorr v=2 vcorr=2;
repeated / type=simple subject=num r=2;
run;
* model for both healthy and asthma groups;
proc mixed data=ACLD.new;
class num diagnose;
model no=diagnose flow flow*diagnose/s cl;
random Intercept flow/ SUBJECT = num g gcorr type= un v vcorr group = diagnose;
repeated / type=simple subject=num;
```

