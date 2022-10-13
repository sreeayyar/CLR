/*
PROGRAM: qobsims_reverse_engineer.do
PURPOSE: calibrate dataset to AK91 

DESCRIPTION: 
1) Adapted the replication code from Mikusheva and Sun (2022)
2) Changed the dgp generating lwklywge to be linear and homoskedastic, 
3) QOB, YOB, POB are obtained from the original sample as we treat the instruments are fixed
*/

set seed 1
set matsize 11000


* load AK91 data (obtain "NEW7080.dta" and apply cleaning code for cohort 30-39 from the Angrist Data Archive):
use "ak91.dta",clear

capture program drop geninsts
program define geninsts
	* generate regressors
	xi i.yob*i.pob,prefix(_YP)

	* generate instruments
	* qob main effects (3)
	xi i.qob,prefix(_Q)
	* qob x yob interactions (27)
	xi i.qob*i.yob,prefix(_QY)
	drop _QYqob_* _QYyob*
	* qob x pob interactions (150)
	xi i.qob*i.pob,prefix(_QP)
	drop _QPqob_* _QPpob_*
	
	* qob x yob x pob interactions (1350)
	tab pob,matrow(pobs)
	local numpobs=rowsof(pobs)
	forvalues q = 2/4 {
		forvalues y = 31/39 {
			forvalues pval = 2/`numpobs' {
				local p = pobs[`pval',1]
				gen _QXYPqobXyobXpob_`q'_`y'_`p' = qob==`q' & yob==`y' & pob==`p'
			}
		}
	}
end

qui geninsts

* define variable sets
local X _YP*
global X _YP*
* define instrument sets
local Z1 _Qqob* _QY* _QP*
local Z2 _Q*
global Z1 _Qqob* _QY* _QP*
global Z2 _Q*

	qui regress educ $Z1 $X // confirm we replicate the first-stage F 
	testparm $Z1
egen cell=group(yob pob)
xtset cell

* run LIML to get coeffs on exogenous regressors
ivregress liml lwklywge (educ=`Z1') `X'
local constant=_b[_cons]
local xb0 `constant'
global xb0 `constant'
local xb1 0
local xb2 0
local ii=0

foreach var of varlist `e(exogr)' {
	local ++ ii
	local j = floor(`ii'/200)
	local beta`var' = _b[`var']
	local xb`j' `xb`j'' + `beta`var''*`var'

}
display `j'


predict uhat,resid
sum uhat
local sigma=r(sd)
global sigma=r(sd)

 * now collapse
collapse (mean) educbar=educ (sd) sigz=uhat,by(qob yob pob)

tempfile fs
save `fs'

use "ak91.dta"
*use `subsample'
keep qob pob yob
* generate instruments
qui geninsts
egen cell=group(yob pob)
xtset cell
	
	qui merge m:1 qob pob yob using `fs'
	qui keep if _merge==3
	*xtset pob
	* generate omitted variable
	gen v = invnormal(uniform())

	* generate schooling
	* gen educmean = educbar+max(-educbar+1,1.7*v)
	* gen educ = rpoisson(educmean)
	scalar pi2 = 0.1
	gen educ = pi2*educbar + v
	
	* generate outcome
	gen epsilon = invnormal(uniform())
	*gen lwklywge =  .1*educ + `xb0' + sigz/`sigma'*(v + .1*epsilon)
	gen lwklywge =  educ + `xb0' + (v + .1*epsilon)
	forvalues ii=1/`j' {
		replace lwklywge = lwklywge + `xb`ii''
	}
	

bys qob yob pob: gen cell_count = _N

// reverse engineer the controls

gen sigz_normalized = sigz/${sigma}

*gen control = lwklywge - 0.1*educ - sigz_normalized * (v + 0.1*epsilon)
gen control = lwklywge - 0.1*educ - (v + 0.1*epsilon)

// empty cells
replace control = 0 if control == .
replace sigz_normalized = 0 if sigz_normalized == .
 
outsheet educbar sigz_normalized control cell_count qob yob pob using "simulation_AK91_no_ZW_homolinear_tenth.dat" if sigz_normalized != 0 , nonames replace
