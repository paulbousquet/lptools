*! regrweight v1.1 - Compute and plot regression weights
*! Syntax: regrweight shock_var [feature1 ... featureK] [if] [in] [, nostand title(string) controls(varlist)]

capture program drop regrweight
program define regrweight
    version 14.0
    syntax varlist(min=1 numeric) [if] [in] [, NOSTAND TITle(string asis) CONTROLS(varlist ts numeric)]
    marksample touse, novarlist

    gettoken shock plotfeat : varlist
    local features `plotfeat' `controls'
    local features : list uniq features

    if "`plotfeat'" == "" {
        local plotfeat `shock'
        local k = 1
    }
    else {
        local k : word count `plotfeat'
        if `k' >= 10 {
            di "Only first 9 will be plotted"
            local k = 9
            local short_plotfeat
            forvalues i = 1/9 {
                local fine : word `i' of `plotfeat'
                local short_plotfeat `short_plotfeat' `fine'
            }
            local plotfeat `short_plotfeat'
        }
    }

    preserve
    keep if `touse'

    quietly count
    if r(N) == 0 {
        di as error "No observations"
        exit 2000
    }

    tempvar shock_work
    quietly gen double `shock_work' = `shock'
    if "`nostand'" == "" {
        tempvar nonzero_only standardized_nonzero standardized_temp random_noise
        quietly gen double `nonzero_only' = `shock' if `shock' != 0
        quietly summarize `nonzero_only'
        if r(N) == 0 {
            di as error "No non-zero values found in `shock'"
            exit 198
        }
        if missing(r(sd)) | r(sd) == 0 {
            di as error "Non-zero standard deviation is zero"
            exit 198
        }
        local mean = r(mean)
        local sd = r(sd)
        quietly gen double `standardized_nonzero' = (`nonzero_only' - `mean') / `sd'
        quietly gen double `standardized_temp' = `shock'
        quietly replace `standardized_temp' = `standardized_nonzero' if `shock' != 0
        quietly replace `shock_work' = `standardized_temp'
        quietly gen double `random_noise' = rnormal(0, .00001)
        quietly replace `shock_work' = `shock_work' + `random_noise' if !missing(`shock_work') & abs(`shock_work') > 0
    }

    local plotfeat_work
    foreach var of local plotfeat {
        if "`var'" == "`shock'" & "`nostand'" == "" {
            local plotfeat_work `plotfeat_work' `shock_work'
        }
        else {
            local plotfeat_work `plotfeat_work' `var'
        }
    }

    local features_work
    foreach var of local features {
        local mapped `var'
        if "`nostand'" == "" {
            if "`var'" == "`shock'" {
                local mapped `shock_work'
            }
            else if regexm("`var'", "\\.`shock'$") {
                local mapped = regexr("`var'", "`shock'$", "`shock_work'")
            }
        }
        local features_work `features_work' `mapped'
    }
    local features_work : list uniq features_work

	local feat_num = 0
    foreach target_var of local plotfeat_work {
        local feat_num = `feat_num' + 1
        di as text _n "Processing feature `feat_num' of `k': `target_var'"
        
        * Create list of all OTHER features (excluding current target)
        local other_features ""
        foreach var of local features_work {
            if "`var'" != "`target_var'" {
                local other_features "`other_features' `var'"
            }
        }
        
        * FWL
        tempvar sample_`feat_num' residuals_`feat_num'
        quietly reg `target_var' `other_features' if !missing(`shock_work')
        quietly gen byte `sample_`feat_num'' = e(sample)
        quietly predict double `residuals_`feat_num'' if `sample_`feat_num'', residuals
		        
        tempfile results_`feat_num'
        quietly summarize `shock_work' if `sample_`feat_num''
        local min_x = r(min) - 0.01
        local max_x = r(max) + 0.01
        quietly postfile handle_`feat_num' level coef using `results_`feat_num'', replace
		
		mata: twirl("`shock_work'", "`residuals_`feat_num''", `max_x', `min_x',"handle_`feat_num'")

		
		quietly postclose handle_`feat_num'
        
        * Drop residuals
        quietly drop `sample_`feat_num'' `residuals_`feat_num''
		
		
	}
	
    
    restore
    
    di as text _n "Creating plots..."
    
    * Determine subplot layout
    if `k' == 1 {
        local rows = 1
        local cols = 1
    }
    else if `k' == 2 {
        local rows = 1
        local cols = 2
    }
    else if `k' == 3 {
        local rows = 1
        local cols = 3
    }
    else if `k' == 4 {
        local rows = 2
        local cols = 2
    }
    else if `k' == 5 | `k' == 6 {
        local rows = 2
        local cols = 3
    }
	else if `k' > 7 {
		local rows = 3
		local cols = 3
	}
    
    local title_len : length local title

    * Create individual plots and combine
    local plot_list ""
    forvalues i = 1/`k' {
        preserve
        
        quietly use `results_`i'', clear

        if `k' == 1 & `title_len' > 0 {
            local plot_title `title'
        }
        else {
            local plot_title `"Weight on x{sub:t} in {&beta}{sub:`i'}"'
        }

        local plot_cmd "line coef level, connect(stairstep) lwidth(thick) lcolor(black)"
        local plot_cmd "`plot_cmd' yline(0, lcolor(black) lpattern(dash) lwidth(thin)) "
        local plot_cmd `"`plot_cmd' title(`plot_title', size(medium))"'
        local plot_cmd "`plot_cmd' xtitle("x{sub:t}", size(small)) ytitle("", size(small)) "
        local plot_cmd "`plot_cmd' graphregion(fcolor(255 255 255)) plotregion(fcolor(255 255 255)) "
        local plot_cmd "`plot_cmd' aspectratio(0.6) "
        local plot_cmd "`plot_cmd' name(plot_`i', replace) nodraw"
        
        quietly `plot_cmd'
        local plot_list "`plot_list' plot_`i'"
        
        restore
    }
    
    if `k' == 1 {
        graph display plot_1
    }
    else {
        if `title_len' == 0 {
            local combine_title `"Regression Weight Functions"'
        }
        else {
            local combine_title `title'
        }
        graph combine `plot_list', rows(`rows') cols(`cols') ///
            graphregion(fcolor(255 255 255)) ///
            title(`combine_title', size(large))
    }
    
end

mata: 

void function twirl( string scalar xvar,
		     string scalar residvar,
		     real scalar xmax,
		     real scalar xmin,
			 string scalar handle_name) {
	 x = st_data(., xvar)
     resid = st_data(., residvar)
     keep = (x :< .) :& (resid :< .)
     x = select(x, keep)
     resid = select(resid, keep)
     ord = order(x, 1)
     x = x[ord]
     resid = resid[ord]
	 n = rows(x)
    
    unique_x = J(0, 1, .)
    unique_pos = J(0, 1, .)  
    last_x = .
    
    for (i = 1; i <= n; i++) {
        if (x[i] != last_x) {
            unique_x = unique_x \ x[i]
            unique_pos = unique_pos \ i
            last_x = x[i]
        }
    }
    
    // Add extended values
    unique_x = xmin \ unique_x \ xmax
    unique_pos = 1 \ unique_pos \ (n + 1) 
    k = rows(unique_x)
    
    // Exploting univariate OLS for efficency 
    cumsum_resid = runningsum(resid)
    total_sum = cumsum_resid[n]
    
    var_resid = variance(resid)
    
    for (j = 1; j <= k; j++) {
        pos = unique_pos[j]
        
        if (pos == 1) {
            sum_resid_above = total_sum
            n_above = n
        } else if (pos > n) {
            sum_resid_above = 0
            n_above = 0
        } else {
            sum_resid_above = total_sum - cumsum_resid[pos - 1]
            n_above = n - pos + 1
        }

        coef = (sum_resid_above / n) / var_resid
        
        stata(sprintf("post %s (%f) (%f)", handle_name, unique_x[j], coef))
    }
			 }


end 
