function [] = figure5_gpufit_cminpack()

%% test parameters
n_graph_points = 3;
snr_min = 3.0;
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
parameters_to_fit = ones(1,n_parameters,'int32')';
user_info = [];
tolerance = 0.0001;

%% parameters determining the randomness of the data
gauss_pos_offset_max = 1.0;
initial_guess_offset_frac = 0.5;

for i = 1:n_graph_points
    
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
    
    
    %% save test results
    results_gpufit_mean_iterations(i) = mean_n_iterations;
    results_gpufit_precision(i) = gpufit_abs_precision;
    
    print_fit_info(results_gpufit_precision(i), time_gpufit, 'Gpufit', numel(valid_indices)/n_fits, results_gpufit_mean_iterations(i));

    %% run cminpack
    [parameters_cminpack, info_cminpack, n_iterations_cminpack, time_cminpack]...
        = cminpack(data, initial_guess_parameters, model_id, tolerance);
    
    converged_cminpack = (info_cminpack > 0) & (info_cminpack <= 3);
    chisquare_cminpack = ones(1,n_fits);

    [valid_cminpack_results, cminpack_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(data_parameters, parameters_cminpack, converged_cminpack, ...
                                     chisquare_cminpack, n_iterations_cminpack);
    
    %% save test results
    results_cminpack_mean_iterations(i) = mean_n_iterations;
    results_cminpack_precision(i) = cminpack_abs_precision;
    
    print_fit_info(results_cminpack_precision(i), time_cminpack, 'cminpack', numel(valid_indices)/n_fits, results_cminpack_mean_iterations(i));

end

%% output filename
filename = 'figure5_gpufit_cminpack';

%% write file
xlsfilename = [filename '.xls'];

Raw(1:100, 1:100)=deal(NaN);
xlswrite(xlsfilename,Raw,1)

xlscolumns = {'SNR' 'gpufit_Precision' 'Minpack_Precision' 'gpufit_Iterations' 'Minpack_Iterations'};
xlswrite(xlsfilename,xlscolumns,1,'A1')

xlsmat(:,1) = snr;
xlsmat(:,2) = [results_gpufit_precision.x0];
xlsmat(:,3) = [results_cminpack_precision.x0];
xlsmat(:,4) = results_gpufit_mean_iterations;
xlsmat(:,5) = results_cminpack_mean_iterations;
xlswrite(xlsfilename,xlsmat,1,'A2')
clear xlsmat

%% plot
Plot_gpufit_cmpipack_precision(snr, [results_gpufit_precision.x0], [results_cminpack_precision.x0])
savefig(filename)
Plot_gpufit_cmpipack_iterations(snr, results_cminpack_mean_iterations, results_gpufit_mean_iterations)
savefig([filename '_iterations'])

end

function [] = Plot_gpufit_cmpipack_iterations(...
    snr,...
    iterations_cminpack,...
    iterations_gpufit)

figure('Name','gpufit vs cminpack, snr presision','NumberTitle','off');
semilogx(...
    snr, iterations_gpufit, 'red+', ...
    snr, iterations_cminpack, 'blueo', ...
    'LineWidth', 4, ...
    'LineStyle', '-', ...
    'MarkerSize', 20)
xlabel('SNR')
ylabel('number of iterations')
legend('gpufit', 'Minpack')
xlim([-inf inf])
grid on;
box off;
current_figure = gca;
current_figure.FontSize = 25;
current_figure.LineWidth = 4;
end

function [] = Plot_gpufit_cmpipack_precision(...
    snr,...
    results_gpufit_precision,...
    precision_cminpack)

figure('Name','gpufit vs cminpack, snr presision','NumberTitle','off');
loglog(...
    snr, results_gpufit_precision, 'red+', ...
    snr, precision_cminpack, 'blues', ...
    'LineWidth', 4, ...
    'LineStyle', '-', ...
    'MarkerSize', 20)
xlabel('SNR')
ylabel('relative standard deviation')
legend(...
    'gpufit center x',...
    'Minpack center x')
xlim([-inf inf])
grid on;
box off;
current_figure = gca;
current_figure.FontSize = 25;
current_figure.LineWidth = 4;

end