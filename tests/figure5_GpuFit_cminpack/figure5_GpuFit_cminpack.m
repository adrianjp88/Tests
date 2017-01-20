function [] = figure5_GpuFit_cminpack()

%% test parameters
n_graph_points = 30;
snr_min = 2;
snr_max = 100000000;
log_min = log10(snr_min);
log_max = log10(snr_max);
log_snr = linspace(log_min, log_max, n_graph_points);
snr = 10.^log_snr;

%% number of fits per test point
n_fits = 1000;

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
parameters_to_fit = ones(1,n_parameters);
user_info = 0;
tolerance = 0.0001;

%% parameters determining the randomness of the data
gauss_pos_offset_max = 2.0;
initial_guess_offset_frac = 0.5;

for i = 1:n_graph_points
    
    fprintf('SNR = %8.3f \n', snr(i));
    
    [data, data_parameters, initial_guess_parameters] = ...
        generate_gauss_fit_test_data(n_fits, fit_size, gauss_amplitude, gauss_width, ...
                                     gauss_baseline, noise, gauss_pos_offset_max, ...
                                     initial_guess_offset_frac, snr(i));
    
    %% run GpuFit
    [parameters_GpuFit, converged_GpuFit, chisquare_GpuFit, n_iterations_GpuFit, time_GpuFit]...
        = GpuFit(n_fits, data, model_id, initial_guess_parameters, weights, tolerance, ...
                 max_iterations, parameters_to_fit, estimator_id, user_info);
    
    converged_GpuFit = converged_GpuFit + 1;
             
    chk_gpulmfit = 0;

    [valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(tmp_data_params, parameters_GpuFit, converged_GpuFit, ...
                                     chisquare_GpuFit, n_iterations_GpuFit, chk_gpulmfit);
    
    
    %% save test results
    results_gpufit_mean_iterations(i) = mean_n_iterations;
    results_gpufit_precision(i) = gpufit_abs_precision;
    
    print_fit_info(results_gpufit_precision(i), time_GpuFit, 'Gpufit', numel(valid_indices)/n_fits, results_gpufit_mean_iterations(i));

    %% run cminpack
    [parameters_cminpack, info_cminpack, n_iterations_cminpack, time_cminpack]...
        = cminpack(data, initial_guess_parameters, model_id, tolerance);
    
    converged_cminpack = (info_cminpack > 0) & (info_cminpack <= 3);
    chisquare_cminpack = ones(1,tmp_n_fits);

    [valid_cminpack_results, cminpack_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(tmp_data_params, parameters_cminpack, converged_cminpack, ...
                                     chisquare_cminpack, n_iterations_cminpack);
    
    %% save test results
    results_cminpack_mean_iterations(i) = mean_n_iterations;
    results_cminpack_precision(i) = cminpack_abs_precision;
    
    print_fit_info(results_cminpack_precision(i), time_cminpack, 'cminpack', numel(valid_indices)/n_fits, results_cminpack_mean_iterations(i));

end

%% output filename
filename = 'figure5_GpuFit_cminpack';

%% write file
xlsfilename = [filename '.xls'];

Raw(1:100, 1:100)=deal(NaN);
xlswrite(xlsfilename,Raw,1)

xlscolumns = {'SNR' 'GpuFit_Precision' 'Minpack_Precision' 'GpuFit_Iterations' 'Minpack_Iterations'};
xlswrite(xlsfilename,xlscolumns,1,'A1')

xlsmat(:,1) = snr;
xlsmat(:,2) = [results_gpufit_precision.x0];
xlsmat(:,3) = [results_cminpack_precision.x0];
xlsmat(:,4) = results_gpufit_mean_iterations;
xlsmat(:,5) = results_cminpack_mean_iterations;
xlswrite(xlsfilename,xlsmat,1,'A2')
clear xlsmat

%% plot
Plot_GpuFit_cmpipack_precision(snr, [results_gpufit_precision.x0], [results_cminpack_precision.x0])
savefig(filename)
Plot_GpuFit_cmpipack_iterations(snr, results_cminpack_mean_iterations, results_gpufit_mean_iterations)
savefig([filename '_iterations'])

end

function [] = Plot_GpuFit_cmpipack_iterations(...
    snr,...
    iterations_cminpack,...
    iterations_GpuFit)

figure('Name','GpuFit vs cminpack, snr presision','NumberTitle','off');
semilogx(...
    snr, iterations_GpuFit, 'red+', ...
    snr, iterations_cminpack, 'blueo', ...
    'LineWidth', 4, ...
    'LineStyle', '-', ...
    'MarkerSize', 20)
xlabel('SNR')
ylabel('number of iterations')
legend('GpuFit', 'Minpack')
xlim([-inf inf])
grid on;
box off;
current_figure = gca;
current_figure.FontSize = 25;
current_figure.LineWidth = 4;
end

function [] = Plot_GpuFit_cmpipack_precision(...
    snr,...
    results_gpufit_precision,...
    precision_cminpack)

figure('Name','GpuFit vs cminpack, snr presision','NumberTitle','off');
loglog(...
    snr, results_gpufit_precision, 'red+', ...
    snr, precision_cminpack, 'blues', ...
    'LineWidth', 4, ...
    'LineStyle', '-', ...
    'MarkerSize', 20)
xlabel('SNR')
ylabel('relative standard deviation')
legend(...
    'GpuFit center x',...
    'Minpack center x')
xlim([-inf inf])
grid on;
box off;
current_figure = gca;
current_figure.FontSize = 25;
current_figure.LineWidth = 4;

end