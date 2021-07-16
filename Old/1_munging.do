*Polygenic risk for body mass index: life course data from infancy to older age
*cleaning and derivation file

**1946c merge genetic and phenotypic data

*home
*global data "S:\IOEQSS_Main\Bann\crosscinequality"
// global data "c:\users\dab\onedrive - university college london\stats_cls\cohorts"
global output "C:\Users\liamj\OneDrive - University College London\prs_bmi\analysis\output"

*global output  "S:\IOEQSS_Main\Bann\crosscinequality\cognition"
use  "$output\46c_bmi_ht_bmiprs.dta", clear
*clean up 1946c
global data "c:\users\dab\onedrive - university college london\stats_cls\cohorts"
use "$data\46c\46c_ht_cog.dta", clear
save "$data\46c\46c_ht_wt_cog_pared.dta", replace

merge 1:1 nshdid_db1120  using "$data\46c\magels_fagels.dta"

clonevar NSHDID_DB1120 = nshdid_db1120  
cap drop _merge
merge 1:1 NSHDID_DB1120  using "$data\46c\bmi_khera2019_prs_withlabels_SCRAMBLED.dta"

*merge scores supplied july 2021 - childhood prs for BMI
cap drop _merge
merge 1:1 NSHDID_DB1120  using "$data\46c\2021_06_NSHD_two_childBMI_prs_scrambled_seqn340.dta"

*bmi prs
*khera et al 2019 - ~adult BMI
desc bmi_prs_khera nmiss_allele_ct named_allele_dosage_sum 
sum bmi_prs_khera nmiss_allele_ct named_allele_dosage_sum 
corr bmi_prs_khera nmiss_allele_ct named_allele_dosage_sum 
 
*br 4167880 //all values appear 4167880
clonevar prs_k = bmi_prs_khera
egen zprs_k =std(prs_k) 
*egen bmiprsg =cut(bmi_prs_khera) , group(5)

*childhood
*vogelezang et al plos gen
desc childbmi_prs_vogelezang2020 effect_allele_sum_vogelezang2020 zchildbmi_prs_vogelezang2020
corr childbmi_prs_vogelezang2020 effect_allele_sum_vogelezang2020 zchildbmi_prs_vogelezang2020 //as expected high correlation

clonevar prs_v = childbmi_prs_vogelezang2020
egen zprs_v =std(childbmi_prs_vogelezang2020) 

*richardson et al 2020 bmj
desc childbmi_prs_richardson2020 effect_allele_sum_richardson2020 zchildbmi_prs_richardson2020 
sum childbmi_prs_richardson2020 effect_allele_sum_richardson2020 zchildbmi_prs_richardson2020
corr childbmi_prs_richardson2020 effect_allele_sum_richardson2020 zchildbmi_prs_richardson2020 //[unexpected weak -ive correlation w/ effect allele count?]

clonevar prs_r = childbmi_prs_richardson2020
egen zprs_r =std(prs_r) 

*label vars
label var prs_k "PRS for BMI, Khera et al"
label var prs_r "PRS for BMI, Richardson et al"
label var prs_v "PRS for BMI, Vogelezang et al"

*label vars
label var zprs_k "zPRS for BMI, Khera et al"
label var zprs_r "zPRS for BMI, Richardson et al"
label var zprs_v "zPRS for BMI, Vogelezang et al"

*quintiles for prs

foreach exposure	in  zprs_k zprs_r zprs_v {
cap drop  q`exposure' 
egen q`exposure' = cut(`exposure'), group(5)
tabstat `exposure' , by(q`exposure' ) stats(n mean min max)
}

*label vars
label var qzprs_k "PRS for BMI quintile, Khera et al"
label var qzprs_r "PRS for BMI quintile, Richardson et al"
label var qzprs_v "PRS for BMI quintile, Vogelezang et al"

*parental height
lookfor height
sum mht52 fht52
mvdecode mht52 fht52, mv(-9999/-9000=.)
sum mht52 fht52

gen mht = mht52 * 0.0254
gen pht = fht52 * 0.0254

sum mht pht
corr mht pht

*med
desc magels
mvdecode magels , mv(-9999=. \ 88/99=. \ 0/9=.)
tab magels

*convert years to 1-10 to match other cohorts 
tab magels , nolab
gen medy = magels -9
tab medy
replace  medy = 10 if medy>10 & !missing(medy)
tab medy

mvdecode fagels , mv(-9999=. \ 88/99=. \ 0/9=.)
tab fagels

tab fagels , nolab
gen fedy = fagels -9

replace  fedy = 10 if fedy>10 & !missing(fedy)
tab fedy

*gen binary 
**based on # years, consistent with the 58 + 70 cohorts 
tab magels //largely concordant 14y age in sam's documentation
tab magels , nolab
recode magels (10/14= 0 "low") (15/23 = 1 "high, stayed on") (0/9=.) (99=.) (-9999=.) (88=.), gen(medb) //DON: 0-9 deemed missing as their labels are not year values
																							//DON: recoded 0/1 values flipped for consistency with NCDS & BCS70
tab magels medb, mi
tab medb
label variable medb "Mother's education level (stayed in school or left)"

recode fagels (10/14= 0 "low") (15/30 = 1 "high, stayed on") (0/9=.) (99=.) (-9999=.) (88=.), gen(fedb) //DON: 0-9 deemed missing as their labels are not year values
tab fagels fedb  

*qualifications
mvdecode med fed , mv(-9899=.)
tab med , nolab
tab med
tab fed

*qualifications
*clonevar  medq = med
*drop med

*clonevar fedq = fed
*drop fed


*tab med medb 

*for now (given can't seem to merge w closer data *sigh*)
*make binary var for this

*fsc
desc fsc57 fsc50

sum fsc50 fsc57

tab  fsc50 
tab  fsc50 , nolab
tab  fsc57 
tab  fsc57 , nolab

cap drop fsc11
recode fsc57  (50=0 "v unskilled") (40=1 ) (35=2 ) (30=3 )  (20=4 )   (10=5 "i professional") (60=.), gen(fsc11)
tab fsc11 fsc57 , mi

cap drop fsc4
recode fsc50  (50=0 "v unskilled") (40=1 ) (35=2 ) (30=3 )  (20=4 )   (10=5 "i professional") (60=.), gen(fsc4)
tab fsc4

*reverse coding so higher score -> higher BMI
vreverse fsc4, generate(fsc4r)
tab fsc4r

vreverse fsc11, generate(fsc11r)
tab fsc11r

drop fsc4  fsc11
rename fsc4r fsc4
rename fsc11r fsc11

cap drop fsc4_11
gen fsc4_11 = fsc4
replace fsc4_11 = fsc11 if missing(fsc4)

sum fsc4_11 if missing(fsc11)
sum fsc11 fsc4

*gen ridit
tab fsc4_11
egen fsc4_11_ridit = ridit(fsc4_11)
egen fsc4_ridit = ridit(fsc4)
egen fsc11_ridit = ridit(fsc11)

*sex
tab sex
tab sex, nolab
rename sex x1
recode x1 (1=1 "men") (2=0 "women"), gen(sex)
tab sex x1, mi
tab sex

*cohort
gen cohort =0

*ht - all in cm
desc ht**u
sum ht**u
mvdecode ht**u, mv(7777/9999=.)
mvdecode htn09, mv(7777/9999=.)
mvdecode htn15x, mv(-99/-1=.)

sum ht*
clonevar ht2 = ht48u
clonevar ht4 = ht50u
clonevar ht6 = ht52u
clonevar ht7 = ht53u
clonevar ht11 = ht57u
clonevar ht15 = ht61u
clonevar ht20 = ht66u
clonevar ht26 = ht72u
clonevar ht36 = ht82u
clonevar ht43 = ht89u
clonevar ht53 = ht99u

gen  ht63 = htn09 * 100
clonevar ht69 = htn15x

sum ht4-ht69

*adult ht
*bmi
desc bmi**u
sum bmi**u
mvdecode bmi**u bmi09, mv(7777/9999=.)
mvdecode bmi15x, mv(-9/-1=.)

clonevar bmi2 = bmi48u
clonevar bmi4 = bmi50u
clonevar bmi6 = bmi52u
clonevar bmi7 = bmi53u
clonevar bmi11 = bmi57u
clonevar bmi15 = bmi61u
clonevar bmi20 = bmi66u
clonevar bmi26 = bmi72u
clonevar bmi36 = bmi82u
clonevar bmi43 = bmi89u
clonevar bmi53 = bmi99u
clonevar bmi63 = bmi09
clonevar bmi69 = bmi15x
sum bmi**

*ages
lookfor age
desc date**c 
sum date**c 
mvdecode date**c , mv(888=.)

cap drop xage*
gen xage2 = 2 //note seemingly no exact age at interview at 2 years, thus assume all interviewed at same age
gen xage4= date50c/12
gen xage6= date52c/12
gen xage7= date53c/12
gen xage11= date57c/12
gen xage15= date61c/12

sum xage*

*ht in cm -> z scores (adj for age and sex)
foreach x in 2 4 6 7 11 15 {
egen zht`x' = zanthro(ht`x', ha, UK) , xvar(xage`x') gender(sex)  gencode(male=1, female=0) 
}
sum zht*

*bmi z scores (adj for age and sex)
foreach x in 2 4 6 7 11 15 {
egen zbmi`x' = zanthro(bmi`x', ba, UK) , xvar(xage`x') gender(sex)  gencode(male=1, female=0) 
}
sum zbmi*


*dxa
***DXA measures
*****Whole body measures - excluding the head
desc dxawbft09 dxawbln09
sum dxawbft09 dxawbln09

mvdecode dxawbft09 dxawbln09, mv(777777=. \ 999999=.)

*whole body - fat
gen fm = dxawbft09*0.001
sum fm

*whole body - lean 
gen lm = dxawbln09*0.001
sum lm

*other components
*android/gynoid
*reg fmi zbmiprs
*reg lmi zbmiprs

*height-adjusted indices
mvdecode htn09 , mv(7777/9999=.)
sum htn09 //ht in m

gen fmi = fm / htn09 ^2
*gen fmi2 = fm / htn09^2

gen lmi = lm / htn09^2

sum fmi lmi

lookfor android

sum dxaangyr09
tab dxaangyr09, nolab
mvdecode dxaangyr09, mv(7/9=.)
clonevar agratio = dxaangyr09
sum agratio

*relative difference - percent?
foreach outcome in 2 4 6 7 11 15 20 26 36 43 53 63 69 {
gen logbmi`outcome' = ln(bmi`outcome')
gen loght`outcome' = ln(ht`outcome')
}

*z score
foreach outcome in 2 4 6 7 11 15 20 26 36 43 53 63 69 {
egen zibmi`outcome' = std(bmi`outcome')
egen ziht`outcome' = std(ht`outcome')
}

**make minumal dataset for analyses
cap drop bmi**u
cap drop ht**u

keep *ht* *bmi* fm* lm* sex    fsc* cohort  mht pht inf med* fed*  *prs* agratio
save  "$output\46c_bmi_ht_bmiprs.dta", replace
