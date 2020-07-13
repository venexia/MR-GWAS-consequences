log using MR-GWAS-consequences

* 21 simulations for different effects of phenotype A on phenotype B

forvalues ab = -1(0.1)1 {
		
	* Specify intitial values for parameters
	
	local at = 1 // Effect of phenotype A on treatment
	local ama = 1 // Effect of phenotype A on measured phenotype A
	local tma = -0.2 // Effect of treatment on measured phenotype A
	local ca = 1 // Effect of covariate on phenotype A
	local cb = 1 // Effect of covariate on phenotype B
	local n = 1000000 // Sample size
	local nsnps = 200 // Number of SNPs
	
	* Set number of observations equal to number of inidviduals in analysis
	
	clear
	set obs `n'

	* Generate SNPs
	
	forvalues i = 1/`nsnps' {
		
		** Each SNP has a random binomial probability
		local a = runiform(0,1)
		qui gen z_`i' = rbinomial(2,`a')
		
		** Each SNP is then multiplied by a random normal variable
		** This essentially weights each SNP by its effect on phenotype A, so when summed it gives a PRS
		local b = rnormal(0,0.1)
		qui replace z_`i' = z_`i' * `b'
		
	}

	* Generate the covariate
	
	gen c = rnormal(0,1)

	* Generate phenotype A 
	
	egen a = rowtotal(z_*)
	replace a = a + rnormal(0,2.5) + c*`ca'

	* Generate a continuous treatment variable
	
	gen logodds_t = rnormal(0,1) + a*`at'
	
	* Convert the continuous treatment variable to a probability
	
	gen p_t = exp(logodds_t) / (1 + exp(logodds_t))

	* Generate a binary treatment indicator using the probability
	gen t = 0
	replace t = 1 if runiform(0,1) < p_t

	* Generate measured phenotype A
	
	gen ma = a*`ama' + t*`tma' + rnormal(0,1)
	
	* Calculate the treatment correction
	
	qui reg ma t a
	di _b[t]
	
	* Correct measured phenotype A for treatment
	
	gen ma_t = ma - t*_b[t]

	* Generate phenotype B
	
	gen b = a*`ab' + c*`cb' + rnormal(0,1)

	* Split sample to prevent sample overlap for phenotype A and B genetic associations
	
	gen split = 0
	replace split = 1 if runiform(0,1) < 0.5

	* Calculate correlation coefficients for all variables with initial values
	
	** Phenotype A - treatment
	corr a t
	gen at = r(rho)
	
	** Phenotype A - measured phenotype A
	corr a ma
	gen ama = r(rho)
	
	** Treatment - measured phenotype A
	corr t ma
	gen tma = r(rho)
	
	** Measured phenotype A - phenotype B
	corr ma b
	gen mab = r(rho)
	
	** Phenotype A - phenotype B	
	corr a b
	gen ab = r(rho)
	
	** Covariate - phenotype A
	corr c a
	gen ca = r(rho)
	
	** Covariate - phenotype B
	corr c b
	gen cb = r(rho)
	
	* Generate variable to store SNP number alongside genetic associations
	
	gen snp = _n in 1/`nsnps'

	* Generate variables to store genetic associations for phenotype A with model number
	
	forvalues i = 0/5 {
		gen za`i' = .
		gen za`i'_se = .
	}

	* Generate variable to store genetic associations for phenotype B

	gen zb = .
	gen zb_se = .
	
	* Calculate and record genetic associations for phenotype B
	
	qui reg b z_* if split == 1
	forvalues i = 1/200 {
		qui replace zb = _b[z_`i'] in `i'
		qui replace zb_se = _se[z_`i'] in `i'
	}

	* Calculate and record genetic associations for phenotype A in model 0

	qui reg ma z_* if split == 0
	forvalues i = 1/`nsnps' {
		qui replace za0 = _b[z_`i'] in `i'
		qui replace za0_se = _se[z_`i'] in `i'
	}

	* Calculate and record genetic associations for phenotype A in model 1

	qui reg ma z_* c if split == 0
	forvalues i = 1/200 {
		qui replace za1 = _b[z_`i'] in `i'
		qui replace za1_se = _se[z_`i'] in `i'
	}

	* Calculate and record genetic associations for phenotype A in model 2

	qui reg ma_t z_* if split == 0
	forvalues i = 1/200 {
		qui replace za2 = _b[z_`i'] in `i'
		qui replace za2_se = _se[z_`i'] in `i'
	}
	
	* Calculate and record genetic associations for phenotype A in model 3
	
	qui reg ma z_* if split == 0 & t == 0
	forvalues i = 1/200 {
		qui replace za3 = _b[z_`i'] in `i'
		qui replace za3_se = _se[z_`i'] in `i'
	}

	* Calculate and record genetic associations for phenotype A in model 4

	qui reg ma_t z_* c if split == 0
	forvalues i = 1/200 {
		qui replace za4 = _b[z_`i'] in `i'
		qui replace za4_se = _se[z_`i'] in `i'
	}
	
	* Calculate and record genetic associations for phenotype A in model 5

	qui reg ma z_* c if split == 0 & t == 0
	forvalues i = 1/200 {
		qui replace za5 = _b[z_`i'] in `i'
		qui replace za5_se = _se[z_`i'] in `i'
	}
	
	* Perform MR-IVW and record results at the bottom of the dataset

	local res = `nsnps' + 1
	
	forvalues i = 0/5 {
		qui mregger zb za`i' [aw=1/(zb_se^2)], ivw heterogi
		qui replace za`i' = _b[za`i'] in `res'
		qui replace za`i'_se = _se[za`i'] in `res'
	}
	
	* Restrict to relevant variables and the MR results
	
	keep za0-za5_se at-cb
	keep in `res'
	
	* Record initial values used
	
	foreach var in at ama tma ab ca cb {
		gen init_`var' = ``var''
	}
		
	* Append this simulation to the previous simulations and save
	
	capture confirm file "data/simulation.dta"
	if !_rc {
		append using data/simulation.dta"
		save "data/simulation.dta", replace
	}
	else {
		save "data/simulation.dta", replace
	}
}

* Load results

use "data/simulation.dta", clear

* Drop when `ab' == 0 as this estimate cannot be scaled

drop if init_ab > -0.01 & init_ab < 0.01

* Generate scaled differences between each model and the initial value

forvalues i = 0/5 {

	gen z`i'_ab = (za`i' - init_ab)/init_ab
	gen z`i'_ab_se = (za`i'_se)/init_ab
	
}

* Meta-analyse scaled differences

local obs = c(N) + 1
local obs2 = c(N) + 6
set obs `obs2'
qui gen model = ""
local i = `obs'
foreach var of varlist z0_ab z1_ab z2_ab z3_ab z4_ab z5_ab {
	metan `var' `var'_se, nograph
	replace z1_ab = r(ES) in `i'
	replace z1_ab_se = r(seES) in `i'
	replace model = "`var'" in `i'
	local i = `i'+1
}

* Plot results

metan z1_ab z1_ab_se in `obs'/`obs2', nooverall nowt lcols(model)

* Save results

outsheet using "output/simulation.csv", comma replace

log off
