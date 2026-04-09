*! regrweight v1.0 - Compute and plot regression weights
*! Syntax: regrweight shock_var feature1 [feature2 ... featureK]

capture program drop regrweight
program define regrweight
    version 14.0
    
    gettoken shock features : 0
	unab features : `features'
    local k : word count `features'
    
	
    * Check number of features
    if `k' == 0 {
        di as error "Error: At least one feature variable must be specified"
        exit 198
    }
	local plotfeat `features'
    if `k' >= 10 {
        di "Only first 9 will be plotted"
		local k = 9
		local plotfeat 
		forvalues i=1/9 {
			local fine : word `i' of `features'
			local plotfeat `plotfeat' `fine'
		}

    }
	
	
    
    preserve
    
    * allign obs with actual regression 
    keep `shock' `features'
    quietly keep if !missing(`shock')
    foreach var of local features {
        quietly keep if !missing(`var')
    }
	
	sort `shock'
	

    quietly summarize `shock'
    local min_x = r(min) - 0.01
    local max_x = r(max) + 0.01
    
	
	local feat_num = 0
    foreach target_var of local plotfeat {
        local feat_num = `feat_num' + 1
        di as text _n "Processing feature `feat_num' of `k': `target_var'"
        
        * Create list of all OTHER features (excluding current target)
        local other_features ""
        foreach var of local features {
            if "`var'" != "`target_var'" {
                local other_features "`other_features' `var'"
            }
        }
        
        * FWL
        quietly reg `target_var' `other_features'
        quietly predict residuals_`feat_num', residuals
		        
        tempfile results_`feat_num'
        quietly postfile handle_`feat_num' level coef using `results_`feat_num'', replace
		
		mata: twirl("`shock'", "residuals_`feat_num'", `max_x', `min_x',"handle_`feat_num'")

		
		quietly postclose handle_`feat_num'
        
        * Drop residuals
        quietly drop residuals_`feat_num'
		
		
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
    
    * Create individual plots and combine
    local plot_list ""
    forvalues i = 1/`k' {
        preserve
        
        quietly use `results_`i'', clear
        
        local plot_cmd "line coef level, connect(stairstep) lwidth(thick) lcolor(black)"
        local plot_cmd "`plot_cmd' yline(0, lcolor(black) lpattern(dash) lwidth(thin)) "
        local plot_cmd "`plot_cmd' title("Weight on x{sub:t} in {&beta}{sub:`i'}", size(medium))"
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
        graph combine `plot_list', rows(`rows') cols(`cols') ///
            graphregion(fcolor(255 255 255)) ///
            title("Regression Weight Functions", size(large))
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