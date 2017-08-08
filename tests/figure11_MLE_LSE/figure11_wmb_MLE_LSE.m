function [] = figure11_wmb_MLE_LSE()

%% test data

n_points = 20;
amp_min = 5;
amp_max = 200;
log_min = log10(amp_min);
log_max = log10(amp_max);
log_amp = linspace(log_min, log_max, n_points);
amp = 10.^log_amp;

%% number of fits per test point
n_fits = 10000;

%% parameters determining the data to be fit
fit_size = 15;
gauss_width = 2.0;
gauss_baseline = 5;
noise = 'poisson';

%% parameters determining how the fit is carried out
weights = [];
max_iterations = 20;
model_id = 1; %GAUSS_2D
n_parameters = 5;
parameters_to_fit = ones(1,n_parameters,'int32')';
user_info = [];
tolerance = 0.0001;

%% parameters determining the randomness of the data
gauss_pos_offset_max = 1.0;
initial_guess_offset_frac = 0.3;
snr = 1.0;

for i = 1:numel(amp)

    tmp_ampl = amp(i);
    
    [data, data_parameters, initial_guess_parameters] = ...
        generate_gauss_fit_test_data(n_fits, fit_size, tmp_ampl, gauss_width, ...
                                     gauss_baseline, noise, gauss_pos_offset_max, ...
                                     initial_guess_offset_frac, snr);
    
    % set the weights array
    weights = [];
                                 
    %% run gpufit MLE
    
    tmp_estimator = 1; 
             
    [parameters_MLE, converged_MLE, chisquare_MLE, n_iterations_MLE, time_MLE]...
        = gpufit(data, weights, model_id, initial_guess_parameters, tolerance, ...
                 max_iterations, parameters_to_fit, tmp_estimator, user_info);
               
             
    converged_MLE = converged_MLE + 1;
             
    chk_gpulmfit = 0;

    [valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(data_parameters, parameters_MLE, converged_MLE, ...
                                     chisquare_MLE, n_iterations_MLE, chk_gpulmfit);
                       
    speed_MLE(i) = n_fits/time_MLE;
    precision_MLE(i) = gpufit_abs_precision;
    mean_iterations_MLE(i) = mean_n_iterations;

    print_fit_info(gpufit_abs_precision, time_MLE, 'MLE', numel(valid_indices)/n_fits, mean_n_iterations);

    %% run gpufit weighted LSE

    % set the weights array
    weights = [];
    
    tmp_estimator = 0; 
             
    [parameters_LSE, converged_LSE, chisquare_LSE, n_iterations_LSE, time_LSE]...
        = gpufit(data, weights, model_id, initial_guess_parameters, tolerance, ...
                 max_iterations, parameters_to_fit, tmp_estimator, user_info);
    
    converged_LSE = converged_LSE + 1;
             
    chk_gpulmfit = 0;

    [valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(data_parameters, parameters_LSE, converged_LSE, ...
                                     chisquare_LSE, n_iterations_LSE, chk_gpulmfit);
                       
    speed_unweighted_LSE(i) = n_fits/time_LSE;
    precision_unweighted_LSE(i) = gpufit_abs_precision;
    mean_iterations_unweighted_LSE(i) = mean_n_iterations;

    print_fit_info(gpufit_abs_precision, time_LSE, 'unweighted LSE', numel(valid_indices)/n_fits, mean_n_iterations);

    
    %% run gpufit weighted LSE

    % set the weights array
    weights = 1.0 ./ max(data,1.0);
    
    tmp_estimator = 0; 
             
    [parameters_LSE, converged_LSE, chisquare_LSE, n_iterations_LSE, time_LSE]...
        = gpufit(data, weights, model_id, initial_guess_parameters, tolerance, ...
                 max_iterations, parameters_to_fit, tmp_estimator, user_info);
    
    converged_LSE = converged_LSE + 1;
             
    chk_gpulmfit = 0;

    [valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(data_parameters, parameters_LSE, converged_LSE, ...
                                     chisquare_LSE, n_iterations_LSE, chk_gpulmfit);
                       
    speed_weighted_LSE(i) = n_fits/time_LSE;
    precision_weighted_LSE(i) = gpufit_abs_precision;
    mean_iterations_weighted_LSE(i) = mean_n_iterations;

    print_fit_info(gpufit_abs_precision, time_LSE, 'weighted LSE', numel(valid_indices)/n_fits, mean_n_iterations);

end

%% output filename
filename = 'figure11_wmb_MLE_LSE';

%% write file precision x0
xlsfilename = [filename '_x0.xls'];

xlscolumns = {'amplitude','MLE x0','unweighted LSE x0','weighted LSE x0'};
xlswrite(xlsfilename,xlscolumns,1,'A1')

xlsmat(:,1) = amp;
xlsmat(:,2) = [precision_MLE.x0];
xlsmat(:,3) = [precision_unweighted_LSE.x0];
xlsmat(:,4) = [precision_weighted_LSE.x0];

xlswrite(xlsfilename,xlsmat,1,'A2')
clear xlsmat

end
