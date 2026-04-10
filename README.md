# lptools: A package to expand and unpack local projections

This repository is for a Stata package in development to provide useful tools to estimate local projections. Two of the commands have their own dedicated repos 

* [Smooth LP](https://github.com/paulbousquet/SmoothLP)
* [LP Decomposition](https://github.com/paulbousquet/lpdecomp)

The third command, `regreweight.ado`, is straightforward to implement. 

## Dependencies

```
ssc install moremata; ssc install bspline
```

## regrweight

The codes replicate the plot of the weights on [Ramey (2011)](https://www.aeaweb.org/articles?id=10.1257/jel.49.3.673) shown in [Kolesár and Plagborg-Møller (2025)](https://www.mikkelpm.com/files/nonlinear_causal.pdf)

```
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
regrweight gov_ramey, title("Ramey (2011) Shock")
```

<img width="1044" height="696" alt="plot_ramey2011_check" src="https://github.com/user-attachments/assets/5b4ecd77-189d-4a6a-ad4f-8ba4deeb67de" />

