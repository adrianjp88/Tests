function [valid_fit_indices] = get_valid_fit_results(converged, data_parameters, fit_results, chisquare)

    discard_frac = 0.1;
    n_std_dev_limit = 5.0;

    quant_chisq = quantile(chisquare, 1.0-discard_frac);
    filtered_chisq = chisquare(find(chisquare <= quant_chisq));
    
    quant_a = quantile(fit_results.a, [discard_frac 1.0-discard_frac]);
    filtered_a = fit_results.a(find(fit_results.a >= quant_a(1) & ...
                                    fit_results.a <= quant_a(2)));
                                    
    quant_x0 = quantile(fit_results.x0, [discard_frac 1.0-discard_frac]);
    filtered_x0 = fit_results.x0(find(fit_results.x0 >= quant_x0(1) & ...
                                      fit_results.x0 <= quant_x0(2)));
    
    quant_y0 = quantile(fit_results.y0, [discard_frac 1.0-discard_frac]);
    filtered_y0 = fit_results.y0(find(fit_results.y0 >= quant_y0(1) & ...
                                      fit_results.y0 <= quant_y0(2)));
                                      
    quant_s = quantile(fit_results.s, [discard_frac 1.0-discard_frac]);
    filtered_s = fit_results.s(find(fit_results.s >= quant_s(1) & ...
                                    fit_results.s <= quant_s(2)));
                                      
    quant_b = quantile(fit_results.b, [discard_frac 1.0-discard_frac]);
    filtered_b = fit_results.b(find(fit_results.b >= quant_b(1) & ...
                                    fit_results.b <= quant_b(2)));
                                      
    mean_chisq = mean(filtered_chisq);
    std_chisq = std(filtered_chisq);
        
    std_a  = std(filtered_a);
    std_x0 = std(filtered_x0);
    std_y0 = std(filtered_y0);
    std_s  = std(filtered_s);
    std_b  = std(filtered_b);
    
    mean_a  = mean(filtered_a);
    mean_x0 = mean(filtered_x0);
    mean_y0 = mean(filtered_y0);
    mean_s  = mean(filtered_s);
    mean_b  = mean(filtered_b);
    
    dev_mean_a  = abs(mean_a - fit_results.a)';
    dev_mean_x0 = abs(mean_x0 - fit_results.x0)';
    dev_mean_y0 = abs(mean_y0 - fit_results.y0)';
    dev_mean_s  = abs(mean_s - fit_results.s)';
    dev_mean_b  = abs(mean_b - fit_results.b)';
    
    chisq_upper_lim = mean_chisq + n_std_dev_limit*std_chisq;
    
    valid_fit_indices = find(converged==1 & ...
                             dev_mean_a  <= n_std_dev_limit*std_a & ...
                             dev_mean_x0 <= n_std_dev_limit*std_x0 & ...
                             dev_mean_y0 <= n_std_dev_limit*std_y0 & ...
                             dev_mean_s  <= n_std_dev_limit*std_s & ...
                             dev_mean_b  <= n_std_dev_limit*std_b & ...
                             chisquare <= chisq_upper_lim);
    
end