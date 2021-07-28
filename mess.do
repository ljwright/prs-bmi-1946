cd "F:\COVID-19 Data\CLS COVID-19\stata\stata13"

* Load Files
tempfile Temp

forval i = 1/3{
	use covid-19_wave`i'_survey_cls.dta, clear
	
	qui ds, not(type string)
	recode `r(varlist)' (min/-1 = .)
	
	rename *, lower
	rename cw`i'_* *
	gen wave = `i'
	
	capture append using `Temp'
	save `Temp', replace
} 


capture program drop prog_total
program define prog_total
	args newvar vlist
	
	capture drop miss
	egen `newvar' = rowtotal(`vlist')
	egen miss = rowmiss(`vlist')
	replace `newvar' = . if miss > 0
	drop miss
end

prog_total gadphq gad2phq2*


use "F:\UKHLS 1-9 and BHPS 1-18, Standard Licence\stata\stata11_se\ukhls_w1\a_indresp.dta", clear
rename a_* *

gen cohort = .
replace cohort = 1 if inrange(birthy, 1928, 1945)
replace cohort = 2 if inrange(birthy, 1946, 1964)
replace cohort = 3 if inrange(birthy, 1965, 1980)
replace cohort = 4 if inrange(birthy, 1981, 1996)
drop if cohort == .
