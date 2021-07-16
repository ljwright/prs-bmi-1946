**Polygenic risk for body mass index: life course data from infancy to older age
**analysis file

global output "C:\Users\DAB\OneDrive - University College London\cross-cohort\prs_bmi\analysis\output"
use  "$output\46c_bmi_ht_bmiprs.dta", clear

*exposures
sum zprs_k zprs_r zprs_v
 
*outcomes
sum bmi*

*desc stats by prs group

*table 1 bmi
cap erase "$output/t1.csv"

*15y + 53y
foreach exposure	in  zprs_k zprs_r zprs_v {
foreach outcome in  bmi15 bmi53 {
estpost tabstat `outcome' if !missing(bmi15), by(q`exposure') stats(mean sd median cv skewness )  nototal
esttab . using "$output/t1.csv", cells("mean(fmt(%9.3g)) sd skewness") append noobs title("`outcome'"_"`exposure'") 
}
}


*output estimates to enable forestplots
foreach outcome in  logbmi {
foreach exposure	in  zprs_k zprs_r zprs_v {
cap erase "$output\b`outcome'_`exposure'.csv"
cap erase "$output\r`outcome'_`exposure'.csv"
}
}
foreach outcome in  logbmi {
foreach exposure	in  zprs_k zprs_r zprs_v {
foreach age in  2 4 6 7 11 15 20 26 36 43 53 63 69 {

reg `outcome'`age' `exposure' i.sex 

est store `outcome'`age'
esttab  `outcome'`age'  using "$output/b`outcome'_`exposure'.csv", ///
cells("b ci_l ci_u")  stats() modelwidth(20) ///
plain nolabel nogaps varwidth (15) nolines  compress append   nolabel  ///
mlabels("`age'") collabels("") nocons keep(`exposure') noobs 

est store `outcome'`age'
esttab  `outcome'`age'  using "$output/r`outcome'_`exposure'.csv", ///
cells("")  stats(r2) modelwidth(20) ///
plain nolabel nogaps varwidth (15) nolines  compress append   nolabel  ///
mlabels("`age'") collabels("") nocons keep("") noobs 

}
}
}



**imports linear regression models in analyses do file and plots them in foresplots
*loop file name unique parts logbmi_zprs_v
foreach x in logbmi_zprs_k logbmi_zprs_r   logbmi_zprs_v   {
*absolute diff
import delimited "$output\b`x'.csv", colrange(1) clear

gen id = _n
rename (v2 v3 v4) (es lci uci)

*save stata format
save  "$output\b`x'_.dta", replace

*import r squared file
import delimited "$output\r`x'.csv", colrange(1) clear
gen id = _n
rename v2 rsquared
merge 1:1 _n using   "$output\b`x'_.dta"

destring , replace
drop if missing(lci)

destring , replace
keep rsquared es lci uci
replace es = es *100
replace lci = lci *100
replace uci = uci *100


replace rsquared = rsquared  *100
format rsquared %9.1fc

cap drop age
gen Age= _n
label define age 1 "2 years" 2 "4 years"	3 "6 years" 4 "7 years" 5 "11 years" 6 "15 years" 7 "20 years" 8 "26 years" 9 "36 years" 10  "43 years" 11 "53 years" 12 "63 years" 13 "69 years"			, replace

label values Age age

cap drop mean_sd 
gen mean_sd = _n
tostring mean_sd , replace
replace mean_sd =  "17.8 (2.5)" if mean_sd=="1" //need check/update
replace mean_sd =  "16.2 (1.6)" if mean_sd=="2"
replace mean_sd =  "15.8 (1.4)" if mean_sd=="3"
replace mean_sd =  "15.8 (1.4)" if mean_sd=="4"
replace mean_sd =  "17.4 (2.3)" if mean_sd=="5"
replace mean_sd =  "20.1 (2.6)" if mean_sd=="6"
replace mean_sd =  "22.2 (2.7)" if mean_sd=="7"
replace mean_sd =  "22.8 (3.0)" if mean_sd=="8"
replace mean_sd =  "24.1 (3.6)" if mean_sd=="9"
replace mean_sd =  "25.4 (4.1)" if mean_sd=="10"
replace mean_sd =  "27.4 (4.7)" if mean_sd=="11"
replace mean_sd =  "27.9 (4.9)" if mean_sd=="12"
replace mean_sd =  "28.2 (5.2)" if mean_sd=="13"

metan es lci uci, nooverall  group1(mean_sd) rcols( mean_sd rsquared )  lcols(Age) graphregion(color(white)) xline(0,lcolor(black)) ///
					xlabel(  0  , 1,  2,  3 ,4 , 5, 6)  force texts(120) dp(1)  nobox pointopt(msize(medsmall)) nohet title("", size(small) color(black) )

*BMI polygenic risk score (z-score) in relation to measured BMI (kg/m{superscript:2}) across life					
*output
graph export "$output\fplot_`x'.tif", replace width(1000)
}			
