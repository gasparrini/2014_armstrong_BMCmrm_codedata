*******************************************************************************
* R code to ilustrate conditional Poisson regression analysis as described in 
*   "Conditional Poisson models: a flexible alternative to conditional logistic
*    case cross-over analysis"
*   Ben Armstrong, Antonio Gasparrini, Aurelio Tobias
*   BMC Medical Research Methodology - 2014
*   http://www.ag-myresearch.com/bmcmrm2014b.html
*
* Uses data from: "Time series regression studies in environmental epidemiology"
*   Bhaskaran et al International Journal of Epidemiology - 2013
* THE ANALYSIS IS AN EXERCISE ONLY. IN PARTICULAR THERE IS POOR CONTROL FOR TEMPERATURE 
*
* 09 05 2014
* For any problem with this code, please contact the authors
*******************************************************************************

use londondataset2002_2006, clear

* DIVIDE THE OZONE VARIABLE BY 10 SO THAT MODEL ESTIMATES REFER TO A 
*   "PER 10ug/m3 INCREASE" (AS PER CONVENTION)
replace ozone = ozone/10
rename ozone ozone10
label var ozone10 "Ozone level in ug/m3 divided by 10"

* CREATE YEAR X MONTH X DOW STRATUM VARIABLE
gen month=month(date)
gen year=year(date)
gen dow=dow(date)
egen stratum_YMD=group(year month dow) 

*** FIT CONDITIONAL POISSON MODEL
xtset stratum_YMD
xtpoisson numdeaths  ozone10 temperature, fe

******************************************************************************
*** PROGRAM TO CORRECT ESTIMATES FOR OVERDISPERSION AFTER USING XTPOISSON, FE
capture program drop xtpoisson_addOD
program def xtpoisson_addOD, eclass
dis _n(1) "Estimate and standard errors corrected for over-dipersion"
tempvar ppred  nonmissxY stratumsumY stratumsumpred pred x2 
qui predict `ppred', nu0       // GIVES PRED COUNT WITHOUT STRATUM EFFECT
local Y `e(depvar)'
local i `e(ivar)'               // STRATUM INDEX VARIABLE
local dfres=e(N)-e(df_m)-e(N_g) // DF OF THE RESIDUALS
qui gen `nonmissxY'=`Y'*(`ppred'!=.)
qui egen `stratumsumY'=sum(`nonmissxY'), by(`i')
qui egen `stratumsumpred'=sum(`ppred'), by(`i')
qui gen `pred'=`ppred'*`stratumsumY'/`stratumsumpred' // RESCALES PRED COUNTS TO MATCH STRATUM SUMS
qui gen `x2'=(`Y'-`pred')^2/(`pred')
qui summ `x2'
local dispers=r(sum)/`dfres'
dis "df: `dfres' ; pearson x2:" %8.1f r(sum) " ; dispersion: " %8.2f `dispers'

matrix B=get(_b)
matrix V=get(VCE)
matrix corrV=V*`dispers'
ereturn scalar dispers=`dispers'
ereturn post B corrV
ereturn display

*STORE PREDICTED COUNTS AND PEARSON RESIDUALS
capture  drop _xtp_pred_count
qui gen _xtp_pred_count = `pred'
capture drop _xtp_pearsonres 
qui gen _xtp_pearsonres= (`x2'^.5)*sign(`Y'-`pred')
capture drop _xtp_devianceres 
qui gen _xtp_devianceres= sqrt( 2*(`Y'*log(`Y'/`pred')-( `Y'-`pred') )) *sign(`Y'-`pred')
end
******************************************************************************

*** NOW ALLOW FOR OVERDISPERSION (function xtpoisson_addOD below, used after xtpoisson)
xtpoisson_addOD 

*** ADD BRUMBACK AUTOCORRELATION ADJUSTMENT (NEED TO HAVE USED xtpoisson_addOD BEFORE)
gen devreslag1=_xtp_devianceres[_n-1] 
xtpoisson numdeaths devreslag1 ozone10 temperature , fe

 * FINALLY ADD ALLOWANCE FOR OVERDISPERSION TO THAT FOR AUTOCORRELATION
xtpoisson_addOD

*** ILLUSTRATION OF ALLOWING FOR VARYING RATE DENOMINATORS
**  FOR THIS WE HAVE IMAGINED AVAILABILITY OF A RELEVANT POPULATION MEASURE CHANGING
**  AT SHORT TIME SCALES (THOUGH ARFICIALLY SPECIFIED HERE AS A CONSTANT, TO DEMONSTRATE CODE)
gen population = 3000000
xtpoisson numdeaths  ozone10 temperature, exp(population) fe


** FURTHER CODE FOR UNCONDITIONAL POISSON AND CONDITIONAL LOGISTIC (CASE CROSSOVER)
**  ANALYSES REPORTED IN THE TEXT

*** FIT UNCONDITIONAL POISSON MODEL
set matsize 500
glm numdeaths i.stratum_YMD  ozone10 temperature , f(poisson)

*** FIT CONDITIONAL LOGISTIC MODEL
* FIRST EXPAND DATA 
save temp, replace
use temp, clear
sort stratum date
gen one=1                        // convenience variable
by stratum: gen origdos=sum(one) // numbers days in strata 1-4 or 1-5
by stratum: egen n_in_stratum = max(origdos)
expand n_in_stratum
sort stratum origdos
by stratum origdos: gen dos=sum(one)    // distribute duplicated days across case-ref sets 
gen caseday=(dos==origdos)              // set indicator for case day
egen ccset=group(year month dow dos) , label
	
* WEIGHT OBSERVATIONS BY N OF DEATHS ON INDEX DAY
gen tempweight=numdeaths*caseday
egen weight=max(tempweight), by(ccset)
drop if weight==0

* CLOGIT ANALYSIS 
clogit  caseday ozone10 temperature [fweight=weight], group(ccset)
