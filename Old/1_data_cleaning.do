* Merge Datasets
cd "S:\IOECLS_Research\Bann\46c"

use "46c_ht_cog.dta", clear

merge 1:1 nshdid_db1120  using "magels_fagels.dta", nogen

rename nshdid_db1120 NSHDID_DB1120 
merge 1:1 NSHDID_DB1120 using "bmi_khera2019_prs_withlabels_SCRAMBLED.dta", nogen
merge 1:1 NSHDID_DB1120 using "2021_06_NSHD_two_childBMI_prs_scrambled_seqn340.dta", nogen
rename NSHDID_DB1120 id

cd "F:\IOE\Paper 1\Data"

* Polygenic Risk Scores
rename bmi_prs_khera prs_k
rename childbmi_prs_vogelezang2020 prs_v
rename childbmi_prs_richardson2020 prs_r
// egen zprs_v =std(childbmi_prs_vogelezang2020)
// ADD QUANTILES

* Probability Weight
rename inf survey_weight

* Female
gen female = sex - 1

* Height & Weight
capture program drop make_ages
program define make_ages
	args stub
	
	foreach var of varlist `stub'*{
		local year = subinstr("`var'", "`stub'", "", .)
		local age = cond(`year' < 46, 54 + `year', `year' - 46)
		local age = cond(`age' < 10, "0`age'", "`age'")
		rename `var' `stub'_`age'
	}
	recode `stub'_* (7777/max = .) (min/-1 = .)
end

foreach l in w h{
	rename `l'tn* `l'eight*
	rename `l'eight15x `l'eight15
	rename `l't*u `l'eight*
	make_ages `l'eight
}
replace height_63 = height_63*100

* BMI
rename bmi*x bmi*
rename bmi*u bmi*
make_ages bmi
// Standardize; Z-Score; Log

* Socio-Economic Position
numlabel fsc50, add
tab fsc50
recode fsc50 (10 = 1) (20 = 2) (30 = 3) (35 = 4) (40 = 5) (50 = 6) (60 = .), gen(sep)
label define sep 	1 "I Professional" 2 "II Intermediate" ///
					3 "III Skilled Non-Manual" 4 "III Skilled Manual" ///
					5 "IV" 6 "V Unskilled"
label values sep sep
// ssc install egenmore
egen sep_ridit = ridit(sep)

* Format
keep id prs_* bmi_* height_* weight_* sep sep_ridit female
order id prs_* bmi_* height_* weight_* sep sep_ridit female
foreach var of varlist *{
	label variable `var'
}
compress
