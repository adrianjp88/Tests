function [valid_fit_results, abs_precision, mean_n_iterations, valid_indices] = ...
    process_gaussian_fit_results(true_data_parameters, raw_fit_results, converged_flag, ...
                                 chi_square, n_iterations, alternate_data_organization)

    if nargin == 6
        chk_gpulmfit = alternate_data_organization;
    else
        chk_gpulmfit = 0;
    end
                             
    n_parameters = 5;

    if chk_gpulmfit == 0 then
    
        fit_results.a  = raw_fit_results(1:n_parameters:end).';
        fit_results.x0 = raw_fit_results(2:n_parameters:end).';
        fit_results.y0 = raw_fit_results(3:n_parameters:end).';
        fit_results.s  = raw_fit_results(4:n_parameters:end).';
        fit_results.b  = raw_fit_results(5:n_parameters:end).';
        
    else
       
        % gpulmfit data organization
        fit_results.a  = raw_fit_results(2:n_parameters:end).';
        fit_results.x0 = raw_fit_results(3:n_parameters:end).';
        fit_results.y0 = raw_fit_results(4:n_parameters:end).';
        fit_results.s  = raw_fit_results(5:n_parameters:end).';
        fit_results.b  = raw_fit_results(1:n_parameters:end).';
        
    end
    
    valid_indices = get_valid_fit_results(converged_flag, true_data_parameters, fit_results, chi_square);

    valid_n_iterations = n_iterations(valid_indices);
    mean_n_iterations = mean(valid_n_iterations);
    
    valid_fit_results.a = fit_results.a(valid_indices);
    valid_fit_results.x0 = fit_results.x0(valid_indices);
    valid_fit_results.y0 = fit_results.y0(valid_indices);
    valid_fit_results.s = fit_results.s(valid_indices);
    valid_fit_results.b = fit_results.b(valid_indices);
    
    abs_precision.a  = std(valid_fit_results.a  - true_data_parameters.a);
    abs_precision.x0 = std(valid_fit_results.x0 - true_data_parameters.x0(valid_indices));
    abs_precision.y0 = std(valid_fit_results.y0 - true_data_parameters.y0(valid_indices));
    abs_precision.s  = std(valid_fit_results.s  - true_data_parameters.s);
    abs_precision.b  = std(valid_fit_results.b  - true_data_parameters.b);
    
end