**# 1. Introdiction to Casual inference

**# 1.1. Problem of regression : Omitted confounders

clear
set seed 123
set obs 1000

gen x1=rnormal()
gen x2=2*x1+rnormal()

gen y=1+1*x1+2*x2+rnormal()

reg y x1 x2

reg y x1
reg y x2

**# 1.2. Strategies of causal inference - Casual diagram

*bad control 1(collider)
* T->Y true coef. : 0
clear
set obs 10000

gen T=rnormal()
gen Y=rnormal()

reg Y T

gen X=2*T+3*Y+rnormal()

reg Y T X

*bad control 2
* T->Y true coef. : 2

clear
set obs 10000

gen T=rnormal()
gen Y=2*T+rnormal()

reg Y T

gen X=3*Y+rnormal()

reg Y T X
 

*backdoor path
* x->y true coef. : 30

clear
set seed 123
set obs 100000
gen c=rnormal()
gen a=rnormal()
gen x=1*c+2*a+rnormal()
gen k=3*a+rnormal()
gen d=5*x+rnormal()
gen y=6*d+4*k+rnormal()
gen f=7*x+rnormal()
gen g=8*d+rnormal()
gen h=9*y+rnormal()

reg y x
reg y x a
reg y x k
reg y x c
reg y x a k
reg y x a k c
reg y x k f
reg y x k g
reg y x k h

ivreg y (x=c)

reg x c
predict res, res
reg y x res

**# 1.3. Randomized Controlled Trials(RCT) - simulation

clear all
set obs 100000
gen x=rnormal()
gen t=.

gen t1=runiform()
gen n=_n

sort t1
replace t=1 in 1/1000

gen t0=runiform() if t!=1
sort t0
replace t=0 in 1/1000

twoway (kdensity x) (kdensity x if t==1) (kdensity x if t==0) ///
		, legend(label(1 "population") label(2 "treatment") ///
		label(3 "control"))

		
**# 2. Regression adjustment
use TE, clear

reg lnsalary i.job_s

tabstat job_yexp, by(job_s)

reg lnsalary job_yexp i.job_s
margins job_s, at(job_yexp=(0(1)48))
marginsplot, noci plotopt(m(i)) addplot( ///
		(scatter lnsalary job_yexp if job_s==0, mc(orange%20) msize(small)) ///
	    (scatter lnsalary job_yexp if job_s==1, mc(green%20) msize(small)))
		
reg lnsalary job_yexp if job_s==0
est store j0
reg lnsalary job_yexp if job_s==1
est store j1

est table j0 j1

twoway ///
lfit lnsalary job_yexp if job_s==0 || lfit lnsalary job_yexp if job_s==1 || ///
scatter lnsalary job_yexp  if job_s==0, mc(orange%20) msize(small) || ///
scatter lnsalary job_yexp if job_s==1,  mc(green%20) msize(small) 

reg lnsalary c.job_yexp##i.job_s
margins job_s, at(job_yexp=(0(1)48))
marginsplot, noci plotopt(m(i)) addplot( ///
		(scatter lnsalary job_yexp if job_s==0, mc(orange%20) msize(small)) ///
	    (scatter lnsalary job_yexp if job_s==1, mc(green%20) msize(small)))

reg lnsalary c.job_yexp##i.job_s
margins r.job_s, at(job_yexp=(0(1)48))
mat t=r(table)
mata t=st_matrix("t")
mata mean(t[1,.]')
marginsplot	, horizontal plotopt(connect(i))


**# 2.1. potential outcome(marginal model)

reg lnsalary job_yexp if job_s==0 , vce(ro)
predict xb0
reg lnsalary job_yexp if job_s==1 , vce(ro)
predict xb1

mean xb0 xb1
nlcom _b[xb1]-_b[xb0]
gen te=xb1-xb0
mean te

reg lnsalary c.job_yexp##i.job_s, vce(ro)
margins r.job_s

*separated regression using gsem
reg lnsalary job_yexp if job_s==0 , vce(ro)
reg lnsalary job_yexp if job_s==1 , vce(ro)

gen lnsalary1=lnsalary if job_s==1
gen lnsalary0=lnsalary if job_s==0

gsem (lnsalary0 <- job_yexp) (lnsalary1 <- job_yexp)
margins , predict(outcome(lnsalary0)) predict(outcome(lnsalary1)) post
margins, coefl
nlcom  _b[2._predict]-_b[1bn._predict]


**# 2.2. teffects ra
teffects ra (lnsalary job_yexp) (gender)

teffects ra (lnsalary job_yexp i.co_type i.incen_policy) ///
	        (gender), aequations
						
reg lnsalary job_yexp i.co_type i.incen_policy  if gender==0
reg lnsalary job_yexp i.co_type i.incen_policy  if gender==1

**# 3. Inverse probability weighting

use TE, clear

**# 3.1. logit model
logit job_s i.edu age i.gender i.co_type
estat classification
lroc
predict pr

**# 3.1. non-normalilized IPW

*non-normalilized IPW
logit job_s i.edu age i.gender i.co_type
predict p

gen ipw=1/p if job_s==1
replace ipw=1/(1-p) if job_s==0

reg lnsalary i.job_s, vce(ro)
reg lnsalary i.job_s [pw=ipw], vce(ro)

*normalilized IPW
logit job_s
predict p_un

logit job_s i.edu age i.gender i.co_type
predict p_norm

gen ipw_norm=p_un/p_norm if job_s==1
replace ipw_norm=(1-p_un)/(1-p_norm) if job_s==0

sum ipw ipw_norm

reg lnsalary i.job_s, vce(ro)
reg lnsalary i.job_s [pw=ipw], vce(ro)
reg lnsalary i.job_s [pw=ipw_norm], vce(ro)

**# 3.2. teffects ipw

teffects ipw (lnsalary) (job_s i.edu age i.gender i.co_type)
teffects ipw (lnsalary) (job_s i.edu age i.gender i.co_type, probit)

teffects ipw (lnsalary) (job_s i.edu age i.gender i.co_type), atet

teffects ipw (lnsalary) (job_s i.edu age i.gender i.co_type), aequations
logit job_s i.edu age i.gender i.co_type

**# 3.3. diagnostics

*Balance check
teffects ipw (lnsalary) (job_s job_s i.edu age i.gender i.co_type)
gen s=e(sample)

tebalance summarize, baseline
tabstat age if s==1, by(job_s) s(mean variance) columns(s) long

tebalance summarize

tebalance density age
tebalance density gender

tebalance overid

* overlap
teoverlap

**# 3.4. trimmed IPW
use "TE.dta", clear

teffects ipw (lnsalary) (job_s i.edu age i.co_type i.gender i.job_status job_yexp)
predict ps*

teoverlap
tebalance summarize

teffects ipw (lnsalary) (job_s i.edu age i.co_type i.gender i.job_status job_yexp) ///
 if ps1<.9 & ps1>.1 & ps2<.9 & ps2>.1 
teoverlap
tebalance summarize
gen s=e(sample)

dtable lnsalary job_s i.edu age i.co_type i.gender i.job_status job_yexp
dtable lnsalary job_s i.edu age i.co_type i.gender i.job_status job_yexp if s==1


**# 3.5. Trade-off between CI and Overlap assumption

logit job_s i.edu age i.gender i.co_type
estat classification
logit job_s i.edu age i.co_type i.gender i.job_status job_yexp
estat classification

teffects ipw (lnsalary) (job_s i.edu age i.gender i.co_type)
predict ps1, ps

teoverlap
graph copy g1, replace

twoway (histogram ps1 if job_s==0, col(stc1%30)) (histogram ps1 if job_s==1, col(stc2%30))
graph copy g2, replace

teffects ipw (lnsalary) (job_s i.edu age i.co_type i.gender i.job_status job_yexp)
predict ps2, ps

teoverlap
graph copy g3, replace

twoway (histogram ps2 if job_s==0, col(stc1%30)) (histogram ps2 if job_s==1, col(stc2%30))
graph copy g4, replace

sum ps1 ps2

graph combine g1 g2 g3 g4


**# 3.  Propensity-score matching(psmatch)
use TE, clear

**# 3.1. teffects psmatch
gsort -job_s
gen i=_n

teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type)
teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type), atet

**# 3.2. finding nearest neighbors
teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type), gen(mat)
predict ps*, ps

logit job_s i.edu age i.gender i.co_type
predict p

list ps1 ps2 p in 1/3

list mat* in 1
list job_s i.edu age i.gender i.co_type mat1 ps* in 1
list job_s i.edu age i.gender i.co_type ps* in 5291
list job_s i.edu age i.gender i.co_type ps* in 5351

drop mat*
drop ps*

teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type), caliper(0.01) osample(omit1)
tab omit1

teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type), caliper(0.001) osample(omit2)
tab omit2

teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type) if omit2==0 

* using sample
teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type)

gen i=_n
teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type), gen(mat)

gen u=.

des mat*, varlist
local v="`r(varlist)'"
dis "`v'"

mata : mat=st_data(. ,tokens("`v'"))
mata : mat[1..10,.]
mata : mat=colshape(mat, 1)
mata : mat
mata : mat=uniqrows(mat)
mata : mat
mata : mat=mat[1..(rows(mat)-1),.]
mata : mat
mata : st_view(u=.,.,"u")
mata : u[mat,1]=J(rows(u[mat,1]),1,1)
mata : mat[1..60,.]
replace u=0 if u==.

tab u job_s



cap program drop psm_use
program define psm_use
	syntax varlist , GENerate(string) [replace]
	if "`replace'"!="" {
		cap drop `generate'
	}
	cap gen `generate'=.
	mata : mat=st_data(. ,tokens("`varlist'"))
	mata : mat=colshape(mat, 1)
	mata : mat=uniqrows(mat)
	mata : mat=mat[1..(rows(mat)-1),.]
	mata : st_view(u=.,.,"`generate'")
	mata : u[mat,1]=J(rows(u[mat,1]),1,1)
	dis "매칭에 사용된 변수를 구분하는 변수 `generate'를 생성함"
end

**# 3.3. Nearest-neighbor(1~k) matching
use TE, clear

teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type), nn(2) gen(mat)
drop mat*

teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type), nn(1) gen(mat)
drop mat*


**# 3.4. Balance, overlap assumption
use TE, clear

teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type)

tebalance summarize, baseline
tebalance summarize

tebalance density
tebalance box

tebalance density age
tebalance box age

teoverlap

logit job_s i.edu age i.gender i.co_type
estat classification
lroc

* model sensitivity
teffects psmatch (lnsalary) (job_s i.edu)
teffects psmatch (lnsalary) (job_s i.edu age)
teffects psmatch (lnsalary) (job_s i.edu age i.co_type i.gender)

**# 3.5. Heterogeneous Treatment effect
use TE, clear

teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type), gen(mat)
predict ps*, ps
predict te, te

twoway (scatter te ps2) (lfit te ps2) (fpfit te ps2)

sum ps2, de

xtile ps_q=ps2 , n(4)
sum  ps2 if ps_q==1
sum ps2, de

teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type)

reg te if ps_q==1, vce(ro)
reg te if ps_q==2, vce(ro)
reg te if ps_q==3, vce(ro)
reg te if ps_q==4, vce(ro)

dtable i.edu age i.gender i.co_type te, by(ps_q, nototal)

logit job_s i.edu age i.gender i.co_type if ps_q==2

teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type) if ps_q==3
tebalance density
tebalance summarize
teoverlap

teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type) if ps_q==4
tebalance density
tebalance summarize
teoverlap


**# 4. augmented IPW(AIPW)
use TE, clear

teffects aipw (lnsalary job_yexp i.co_type i.incen_policy) ///
			  (job_s i.edu age i.gender i.co_type)

teoverlap
tebalance summarize
tebalance overid
tebalance density age
tebalance density edu

logit job_s i.edu age i.gender i.co_type
predict p

gen ipw1=1.job_s/p
gen ipw0=0.job_s/(1-p)

reg lnsalary job_yexp i.co_type i.incen_policy if job_s==1
predict poms1

reg lnsalary job_yexp i.co_type i.incen_policy if job_s==0
predict poms0

gen pom1_adj=poms1+ipw1*(lnsalary-poms1)
gen pom0_adj=poms0+ipw0*(lnsalary-poms0)

mean pom1_adj pom0_adj
nlcom _b[pom1_adj]-_b[pom0_adj]

teffects aipw (lnsalary job_yexp i.co_type i.incen_policy) ///
			  (job_s i.edu age i.gender i.co_type), aequations

* model sensitivity
teffects aipw (lnsalary job_yexp i.co_type i.gender i.edu i.incen_policy) ///
			  (job_s i.edu)
teffects aipw (lnsalary job_yexp i.co_type i.gender i.edu i.incen_policy) ///
			  (job_s i.edu age)
teffects aipw (lnsalary job_yexp i.co_type i.gender i.edu i.incen_policy) ///
			  (job_s i.edu age i.co_type i.gender)
	
	
teffects aipw (lnsalary job_yexp) (job_s i.edu age i.co_type i.gender)
teffects aipw (lnsalary job_yexp i.co_type) (job_s i.edu age i.co_type i.gender)
teffects aipw (lnsalary job_yexp i.co_type i.gender i.edu i.incen_policy) ///
			  (job_s i.edu age i.co_type i.gender)
	

**# 5. IPW regression adjustment(IPWRA)
use TE, clear

teffects ipwra (lnsalary job_yexp i.co_type i.incen_policy) ///
			   (job_s i.edu age i.gender i.co_type)


teoverlap
tebalance summarize
tebalance overid
tebalance density age
tebalance density edu


logit job_s i.edu age i.co_type i.gender
predict p

gen ipw1=1.job_s/p
gen ipw0=0.job_s/(1-p)

reg lnsalary job_yexp i.co_type i.incen_policy [pw=ipw1] if job_s==1
predict xb1

reg lnsalary job_yexp i.co_type i.incen_policy [pw=ipw0] if job_s==0
predict xb0

mean xb1 xb0
nlcom _b[xb1]-_b[xb0]


* model sensitivity
teffects ipwra (lnsalary job_yexp i.co_type i.gender i.edu i.incen_policy) ///
			   (job_s i.edu)
teffects ipwra (lnsalary job_yexp i.co_type i.gender i.edu i.incen_policy) ///
			   (job_s i.edu age)
teffects ipwra (lnsalary job_yexp i.co_type i.gender i.edu i.incen_policy) ///
			   (job_s i.edu age i.co_type i.gender)
	
	
teffects ipwra (lnsalary job_yexp) (job_s i.edu age i.co_type i.gender)
teffects ipwra (lnsalary job_yexp i.co_type) (job_s i.edu age i.co_type i.gender)
teffects ipwra (lnsalary job_yexp i.co_type i.gender i.edu i.incen_policy) ///
			   (job_s i.edu age i.co_type i.gender)

			   
**# [Appendix] model comparision
use TE, clear

* outcome & treatment variable
teffects ra (lnsalary job_yexp i.co_type i.incen_policy) (job_s)
est store RA	   
teffects ipw (lnsalary) (job_s i.edu age i.gender i.co_type)			
est store IPW
teffects psmatch (lnsalary) (job_s i.edu age i.gender i.co_type)			   
est store PS
teffects aipw (lnsalary job_yexp i.co_type i.incen_policy) ///
			  (job_s i.edu age i.gender i.co_type)			   
est store AIPW
teffects ipwra (lnsalary job_yexp i.co_type i.incen_policy) ///
			   (job_s i.edu age i.gender i.co_type)			   
est store IPWRA

est table RA IPW PS AIPW IPWRA, se


**# [Appendix] Mutivalued treatment effects
use TE, clear

tab job_status

teffects ipw (lnsalary) (job_status i.edu age i.gender i.co_type), osample(ipw_o)

teffects ipw (lnsalary) (job_status i.edu age i.gender i.co_type) if ipw_o==0, aequations
mlogit job_status i.edu age i.gender i.co_type if ipw_o==0

teffects ipw (lnsalary) (job_status i.edu age i.gender i.co_type) if ipw_o==0
teoverlap
tebalance summarize
tebalance density job_yexp
tebalance density edu

teffects aipw (lnsalary job_yexp i.co_type i.incen_policy) ///
			  (job_status i.edu age i.gender i.co_type) if ipw_o==0

teffects ipwra (lnsalary job_yexp i.co_type i.incen_policy) ///
			  (job_status i.edu age i.gender i.co_type) if ipw_o==0

**# 6. lasso for treatment effects

**# 6.1. telasso
webuse lung, clear
vl set
vl create cvars = vlcontinuous - (fev1p)
vl create fvars = vlcategorical - (transtype)
vl sub allvars = c.cvars i.fvars
telasso (fev1p $allvars) (transtype $allvars)
lassoinfo

lassoknots, for(fev1p) tlevel(0)
lassoknots, for(fev1p) tlevel(1)
lassoknots, for(transtype)

teoverlap

logit transtype walkdist lungpo2 ischemict 

lassocoef (., for(transtype)) (., for(fev1p) tlevel(0)) (., for(fev1p) tlevel(1)), dis(coef, postselection)
lassocoef (., for(transtype)) (., for(fev1p) tlevel(0)) (., for(fev1p) tlevel(1)), dis(coef)

**# 6.2. telasso for lnsalary
use TE, clear

vl clear
vl set

vl list vlcontinuous
vl list vlcategorical
vl list vluncertain


vl move (KSIC) vlcategorical
vl move (age job_yexp birth_year org_involvement) vlcontinuous

vl create indvars_d=(edu age gender cohort)
vl create covars_d=(KSIC co_size)
vl create jobvars_d=(incen_policy)
vl create jobvars_c=(job_yexp)

vl substitute inter1=i.indvars_d##i.jobvars_d
vl substitute inter2=i.indvars_d##c.jobvars_c
vl substitute inter3=i.jobvars_d##c.jobvars_c
vl substitute co=i.covars_d

telasso (lnsalary $inter1 $inter2 $inter3 $co) ///
		(job_s $inter1 $inter2 $inter3 $co)

lassocoef (., for(lnsalary) tlevel(0)) (., for(lnsalary) tlevel(1))	(., for(job_s))
		
tebalance summarize
teoverlap

*double machine learning
telasso (lnsalary $inter1 $inter2 $inter3 $co) ///
		(job_s $inter1 $inter2 $inter3 $co) , xfold(5) resample(3)
		
tebalance summarize
teoverlap

lassocoef (., for(job_s) xfold(1) resample(1)) ///
		  (., for(job_s) xfold(2) resample(1)) ///
		  (., for(job_s) xfold(3) resample(1)) ///
		  (., for(job_s) xfold(4) resample(1)) ///
		  (., for(job_s) xfold(5) resample(1))

lassocoef (., for(lnsalary) tlevel(1) xfold(1) resample(1)) ///
		  (., for(lnsalary) tlevel(1) xfold(2) resample(1)) ///
		  (., for(lnsalary) tlevel(1) xfold(3) resample(1)) ///
		  (., for(lnsalary) tlevel(1) xfold(4) resample(1)) ///
		  (., for(lnsalary) tlevel(1) xfold(5) resample(1))
		  
**# 7. Heckman-style models

* 7.1. heckman selection model  
webuse womenwk, clear

heckman wage educ age, select(married children educ age)
heckman wage educ age, select(married children educ age) two

heckman, coefl
dis 5.9473529*0.67284

gen selection =  wage < .
qui probit select married children educ age
predict xb, xb

gen pdf=normalden(xb)
gen cdf=normal(xb)
gen imr=pdf/cdf
reg wage educ age imr
reg wage educ age


generate selected = 0 if wage < .
generate notselected = 0 if wage >= .
gsem (wage <- educ age L) (selected <- married children educ age L@1, ///
	 family(gaussian, udepvar(notselected))), var(L@1 e.wage@a e.selected@a)

nlcom (sigma: sqrt(_b[/var(e.wage)] +_b[wage:L]^2)) ///
      (rho: _b[wage:L]/(sqrt((_b[/var(e.wage)]+1)*(_b[/var(e.wage)]+_b[wage:L]^2))))

nlcom (married: _b[selected:married]/sqrt(_b[/var(e.wage)]+1)) ///
	  (children: _b[selected:children]/sqrt(_b[/var(e.wage)]+1)) ///
	  (educ: _b[selected:educ]/sqrt(_b[/var(e.wage)]+1)) ///
	  (age: _b[selected:age]/sqrt(_b[/var(e.wage)]+1))
	  
* 7.2. endogenous treatment effects
use TE, clear

* Joint normal
teffects ipw (salary) (job_s i.edu age i.gender)

etregress salary job_yexp i.co_type i.incen_policy ///
	     , treat(job_s=i.edu age i.gender)
margins r.job_s

etregress salary i.job_s##(c.job_yexp i.co_type i.incen_policy) ///
	     , treat(job_s=i.edu age i.gender)
margins r.job_s

eregress salary c.job_yexp i.co_type i.incen_policy ///
	     , entreat(job_s=i.edu age i.gender)
estat teffects

etregress salary i.job_s##(c.job_yexp i.co_type i.incen_policy) ///
	     , treat(job_s=i.edu age i.gender)
margins r.job_s

		 
*Control function approach
eteffects (salary job_yexp i.co_type i.incen_policy) ///
		  (job_s i.edu age i.gender), aequations

etregress salary c.job_yexp i.co_type i.incen_policy ///
	     , treat(job_s=i.edu age i.gender) cfunction

etregress salary c.job_yexp i.co_type i.incen_policy ///
	     , treat(job_s=i.edu age i.gender) cfunction  hazard(hr)
gen e=e(sample)
reg salary c.job_yexp i.co_type i.incen_policy i.job_s hr if e==1

probit job_s i.edu age i.gender if e==1
predict xb, xb

gen pdf=normalden(xb)
gen cdf=normal(xb)
gen imr1=pdf/cdf
gen imr0=(-pdf)/(1-cdf)

gen imr=imr1 if job_s==1
replace imr=imr0 if job_s==0

reg salary c.job_yexp i.co_type i.incen_policy i.job_s hr if e==1
reg salary c.job_yexp i.co_type i.incen_policy i.job_s imr if e==1
		  
**# 7.3. extended regression
webuse class10re, clear

heckman gpa income, select(graduate=income)
eregress gpa income, select(graduate=income)

etregress gpa income, treat(program = income i.hscomp)
eregress gpa c.income, entreat(program = income i.hscomp,nointeract)
estat teffects

eregress gpa c.income, entreat(program = income i.hscomp,povariance)
estat teffects
eregress gpa c.income, entreat(program = income i.hscomp,povariance pocorrelation)
estat teffects

eregress gpa income, select(graduate=income) entreat(program = income i.hscomp)
estat teffects


etregress gpa income, treat(program = income i.hscomp)

use https://www.stata-press.com/data/r18/gsem_union3, clear
generate llunion = 0 if union == 1
generate ulunion = 0 if union == 0

gsem (wage <- age grade i.smsa i.black tenure 1.union L) ///
	 (llunion <- i.black tenure i.south L@1, ///
	 family(gaussian, udepvar(ulunion))), ///
	 var(L@1 e.wage@a e.llunion@a)

eregress wage age grade i.smsa i.black tenure, entreat(union = i.black tenure i.south, nointeract)

	 
*xteregress 
webuse class10re, clear
eregress gpa income, entreat(program = income i.hscomp, nointeract)
gen e=e(sample)

gen lprogram = 0 if program == 1
gen uprogram = 0 if program == 0

gsem (gpa <- income 1.program L) (lprogram <- income i.hscomp L@1, fam(gaussian, udepvar(uprogram))) ///
	 if e==1 , var(L@1 e.gpa@a e.lprogram@a)


xteregress gpa income, select(graduate=income i.program)
xtheckman gpa income, select(graduate=income i.program)


gen gpa1=gpa if program==1
gen gpa0=gpa if program==0

gsem (gpa1 <- c.income U1[collegeid]) (gpa0 <- c.income U1[collegeid]) ///
	 (program <- income i.hscomp U[collegeid]@1, probit), ///
	 cov(U1[collegeid]*U[collegeid]) vce(ro)
margins, predict(outcome(gpa1)) predict(outcome(gpa0)) vce(unconditional) post
margins, coefl
nlcom _b[2._predict]-_b[1bn._predict]

xteregress gpa c.income, entreat(program = income i.hscomp)
estat teffects

 
**# 8. Difference-in-differences estimation

**# 8.1. DID for 2X2 data
net get diff.pkg, from(http://fmwww.bc.edu/RePEc/bocode/d/)

use cardkrueger1994.dta, clear

reg fte i.treated##i.t

*pre
dis "E(Y|D=0,T=0)="_b[_cons]
dis "E(Y|D=1,T=0)="_b[_cons]+_b[1.treated]

*pre diff.
dis "E(Y|D=1,T=0)-E(Y|D=0,T=0)="_b[1.treated]

*post
dis "E(Y|D=0,T=1)="_b[_cons]+_b[1.t]
dis "E(Y|D=1,T=1)="_b[_cons]+_b[1.treated]+_b[1.t]+_b[1.treated#1.t]

*post diff.
dis "E(Y|D=1,T=1)-E(Y|D=0,T=1)="_b[1.treated]+_b[1.treated#1.t]

*did
dis "[E(Y|D=1,T=1)-E(Y|D=0,T=1)]-[E(Y|D=1,T=0)-E(Y|D=0,T=0)]="_b[1.treated#1.t]

net install st0424.pkg
diff fte, t(treated) p(t)

**# 8.2. DID for Panel data

webuse parallelt, clear

tab t1 treated1
bys id1 : egen d1=max(treated1)
line treated1 t1, sort xlabel(1(1)10) ylabel(0 1 2) by(d1)

xtdidregress (y1) (treated1), group(id1) time(t1)

xtset id1 t1
xtreg y1 i.treated1 i.t1, fe vce(cl id1)

xtdidregress (y1) (treated1), group(id1) time(t1) aequations
estat trendplots
estat ptrends
estat granger

estat grangerplot, verbose lagview post
test  _lead2  _lead3  _lead4 _lead5

xtset id3
xtdidregress (y3) (treated3), group(id3) time(t3)
estat trendplots
estat ptrends
estat granger


**# 8.3. Heterogeneous DID for Panel Data
webuse akc, clear
xtset breed year

drop tyear

bys breed (year) : gen dsum=sum(movie)
gen tyear=year if dsum==movie & dsum!=0

bys breed (tyear) : replace tyear=tyear[1] if tyear==.
replace tyear=0 if tyear==.
tab tyear

line movie year, sort by(tyear) xlabel(2031(1)2040) ylabel(0 1 2)

xtdidregress (registered) (movie) , group(breed) time(year)
estat bdecomp, graph

xthdidregress twfe (registered) (movie) , group(breed)
estat atetplot
estat aggregation
estat aggregation, cohort graph
estat aggregation, time graph
estat aggregation, dynamic graph

xthdidregress aipw (registered) (movie) , group(breed)

**# 8.4. Heterogeneous DID for Repeated-cross section data

webuse hhabits, clear

hdidregress aipw (bmi medu i.girl i.sports) (hhabit parksd), group(schools) time(year)
estat aggregation, cohort graph
estat aggregation, time graph
estat aggregation, dynamic graph
estat ptrends

hdidregress twfe (bmi medu i.girl i.sports) (hhabit), group(schools) time(year)
estat ptrends

**# [Appendix] synth2
net install st0722.pkg
net get st0722.pkg

use smoking

xtset state year
synth2 cigsale lnincome age15to24 retprice beer cigsale ///
	   , trunit(3) trperiod(1989) xperiod(1980(1)1988) nested allopt





