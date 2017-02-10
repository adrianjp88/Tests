function [] = figure9_gpufit_GPULMFit_cminpack_precision()

%% test parameters
n_graph_points = 30;
snr_min = 2;
snr_max = 100000000;
log_min = log10(snr_min);
log_max = log10(snr_max);
log_snr = linspace(log_min, log_max, n_graph_points);
snr = 10.^log_snr;

%% number of fits per test point
n_fits = 10000;

%% parameters determining the data to be fit
fit_size = 15;
gauss_amplitude = 500;
gauss_width = 1.0;
gauss_baseline = 10;
noise = 'gauss';

%% parameters determining how the fit is carried out
weights = [];
max_iterations = 20;
model_id = 1; %GAUSS_2D
estimator_id = 0; %LSE
n_parameters = 5;
parameters_to_fit = ones(1,n_parameters)';
user_info = [];
tolerance = 0.0001;

%% parameters determining the randomness of the data
gauss_pos_offset_max = 1.0;
initial_guess_offset_frac = 0.3;



%% loop
for i = 1:numel(snr)
    
    fprintf('SNR = %8.3f \n', snr(i));
    
    [data, data_parameters, initial_guess_parameters] = ...
        generate_gauss_fit_test_data(n_fits, fit_size, gauss_amplitude, gauss_width, ...
                                     gauss_baseline, noise, gauss_pos_offset_max, ...
                                     initial_guess_offset_frac, snr(i));
    

    
    %% run gpufit
    
    [parameters_gpufit, converged_gpufit, chisquare_gpufit, n_iterations_gpufit, time_gpufit]...
        = gpufit(data, weights, model_id, initial_guess_parameters, tolerance, ...
                 max_iterations, parameters_to_fit, estimator_id, user_info);
    
    converged_gpufit = converged_gpufit + 1;
             
    chk_gpulmfit = 0;

    [valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(data_parameters, parameters_gpufit, converged_gpufit, ...
                                     chisquare_gpufit, n_iterations_gpufit, chk_gpulmfit);
                       
    speed_gpufit(i) = n_fits/time_gpufit;
    precision_gpufit(i) = gpufit_abs_precision;
    mean_iterations_gpufit(i) = mean_n_iterations;

    print_fit_info(gpufit_abs_precision, time_gpufit, 'Gpufit', numel(valid_indices)/n_fits, mean_n_iterations);

                                 
    %% run GPU-LMFit
    
    [parameters_GPULMFit, info_GPULMFit, time_GPULMFit] = GPULMFit(...
        data,...
        estimator_id,...
        data_parameters.s,...
        fit_size);
    
    
    converged_GPULMFit = ones(1,n_fits);
    chisquare_GPULMFit = ones(1,n_fits);
    n_iterations_GPULMFit = ones(1,n_fits);
    chk_gpulmfit = 1;

    [valid_gpulmfit_results, gpulmfit_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(data_parameters, parameters_GPULMFit, converged_GPULMFit, ...
                                     chisquare_GPULMFit, n_iterations_GPULMFit, chk_gpulmfit);

    speed_GPULMFit(i) = n_fits/time_GPULMFit;
    precision_GPULMFit(i) = gpulmfit_abs_precision;

    print_fit_info(gpulmfit_abs_precision, time_GPULMFit, 'GPU-LMFit', numel(valid_indices)/n_fits, mean_n_iterations);

    
    %% run cminpack
    [parameters_cminpack, info_cminpack, n_iterations_cminpack, time_cminpack]...
        = cminpack(data, initial_guess_parameters, model_id, tolerance);
    
    converged_cminpack = (info_cminpack > 0) & (info_cminpack <= 3);
    chisquare_cminpack = ones(1,n_fits);

    [valid_cminpack_results, cminpack_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(data_parameters, parameters_cminpack, converged_cminpack, ...
                                     chisquare_cminpack, n_iterations_cminpack);

    %% save test results
    speed_cminpack(i) = n_fits/time_cminpack;
    precision_cminpack(i) = cminpack_abs_precision;
    mean_iterations_cminpack(i) = mean_n_iterations;
    print_fit_info(cminpack_abs_precision, time_cminpack, 'C Minpack', numel(valid_indices)/n_fits, mean_n_iterations);

end


%% output filename
filename = 'figure9_gpufit_GPULMFit_cminpack_precision';

%% write .xls file
xlsfilename = [filename '.xls'];

Raw(1:100, 1:100)=deal(NaN);
xlswrite(xlsfilename,Raw,1)

xlscolumns = {'SNR' 'gpufit' 'GPULMFit' 'Minpack'};
xlswrite(xlsfilename,xlscolumns,1,'A1')

xlsmat(:,1) = snr;
xlsmat(:,2) = [precision_gpufit.x0];
xlsmat(:,3) = [precision_GPULMFit.x0];
xlsmat(:,4) = [precision_cminpack.x0];
xlswrite(xlsfilename,xlsmat,1,'A2')


%% plot
Plot_gpufit_GPULMFit_cmpipack_precision(...
    snr, ...
    [precision_gpufit.x0],...
    [precision_GPULMFit.x0],...
    [precision_cminpack.x0])
savefig(filename)

end

function [] = Plot_gpufit_GPULMFit_cmpipack_precision(...
    snr,...
    precision_gpufit,...
    precision_GPULMFit,...
    precision_cminpack)

figure('Name','gpufit vs GPULMFit vs Minpack','NumberTitle','off');
loglog(...
    snr, precision_gpufit, 'red+-', ...
    snr, precision_GPULMFit, 'blues-', ...
     snr, precision_cminpack, 'greeno-', ...
    'LineWidth', 6, ...
    'LineStyle', 'none', ...
    'MarkerSize', 30)
xlabel('snr')
ylabel('relative standard deviation')
legend('gpufit', 'GPU-LMFit', 'Minpack')
xlim([-inf inf])

grid on;
box off;
current_figure = gca;
current_figure.FontSize = 37;
current_figure.LineWidth = 4;
end
