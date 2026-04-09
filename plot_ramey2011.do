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
regrweight gov_ramey gov_ramey
