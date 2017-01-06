function [valid_fit_indices] = get_valid_fit_indices(converged, chisquare)

    mean_chisq = mean(chisquare);
    std_chisq = std(chisquare);
    chisq_upper_lim = mean_chisq + 5.0*std_chisq;
    
    valid_fit_indices = find(converged==1 & chisquare <= chisq_upper_lim);
    
end