set more off
import excel homgovdat.xlsx, sheet("govdat") firstrow clear
gen qdate=q(1947q1)+_n-1
tsset qdate, q
gen t=_n
gen t2=t^2
gen rgdp=100*ngdp/pgdp
gen rgov=100*ngov/pgdp
gen tax=nfedreceipts/ngdp
gen lrgdp=ln(rgdp)
reg lrgdp t t2
predict lyquad
gen yquad=exp(lyquad)
gen newsy=100*rameynews/(L.yquad*L.pgdp)
gen y=rgdp/yquad
gen g=rgov/yquad
reg newsy L(1/2).newsy L(1/2).y L(1/2).g L(1/2).tax t t2
predict gov_ramey, resid
su gov_ramey
replace gov_ramey=gov_ramey/r(sd)
gen max0=max(gov_ramey,0) if gov_ramey<.
reg max0 gov_ramey, robust
local w : di %4.3f _b[gov_ramey]
su gov_ramey, meanonly
levelsof gov_ramey if gov_ramey<., local(x)
local x "`=`r(min)'-0.01' `x' `=`r(max)'+0.01'"
tempfile f
postfile h x b using `f', replace
foreach z of local x {
    gen ind=gov_ramey>=`z' if gov_ramey<.
    reg ind gov_ramey, robust
    post h (`z') (_b[gov_ramey])
    drop ind
}
postclose h
use `f', clear
sort x
line b x, connect(stairstep) lwidth(thick) title("Ramey (2011) Shock: ({&omega}>0: `w')") xtitle("Standard Deviations") ytitle("") graphregion(fcolor(255 255 245))
