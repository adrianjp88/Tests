function [] = figure5_GpuFit_cminpack()

%% test parameters
n_graph_points = 20;
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


%% test setup

gauss_xpos_mean = (fit_size-1)/2;
gauss_ypos_mean = (fit_size-1)/2;

gauss_pos_x = gauss_xpos_mean + ( 2.0*gauss_pos_offset_max*(rand(n_fits, 1) - 0.5) );
gauss_pos_y = gauss_ypos_mean + ( 2.0*gauss_pos_offset_max*(rand(n_fits, 1) - 0.5) );

initial_guess_xpos_offset_max = initial_guess_offset_frac * gauss_width;
initial_guess_ypos_offset_max = initial_guess_offset_frac * gauss_width;
initial_guess_ampl_offset_max = initial_guess_offset_frac * gauss_amplitude;
initial_guess_width_offset_max = initial_guess_offset_frac * gauss_width;
initial_guess_baseline_offset_max = initial_guess_offset_frac * gauss_baseline;

initial_guess_xpos_offset = initial_guess_xpos_offset_max ...
                            * 2.0 * (rand(n_fits, 1) - 0.5);
                  
initial_guess_ypos_offset = initial_guess_ypos_offset_max ...
                            * 2.0 * (rand(n_fits, 1) - 0.5);

initial_guess_ampl_offset = initial_guess_ampl_offset_max ...
                            * 2.0 * (rand(n_fits, 1) - 0.5);
                        
initial_guess_width_offset = initial_guess_width_offset_max ...
                             * 2.0 * (rand(n_fits, 1) - 0.5);
                        
initial_guess_baseline_offset = initial_guess_baseline_offset_max ...
                                * 2.0 * (rand(n_fits, 1) - 0.5);
                        
initial_guess_xpos = gauss_pos_x + initial_guess_xpos_offset;
initial_guess_ypos = gauss_pos_y + initial_guess_ypos_offset;
initial_guess_ampl = gauss_amplitude + initial_guess_ampl_offset; 
initial_guess_width = gauss_width + initial_guess_width_offset; 
initial_guess_baseline = gauss_baseline + initial_guess_baseline_offset;

data_parameters.a = gauss_amplitude;
data_parameters.x0 = gauss_pos_x;
data_parameters.y0 = gauss_pos_y;
data_parameters.s = gauss_width;
data_parameters.b = gauss_baseline;

initial_guess_parameters = ones(1,n_parameters*n_fits);

for i = 1:n_fits
    tmp_index = (i-1) * n_parameters;
    initial_guess_parameters(1 + tmp_index) = initial_guess_ampl(i);
    initial_guess_parameters(2 + tmp_index) = initial_guess_xpos(i);
    initial_guess_parameters(3 + tmp_index) = initial_guess_ypos(i);
    initial_guess_parameters(4 + tmp_index) = initial_guess_width(i);
    initial_guess_parameters(5 + tmp_index) = initial_guess_baseline(i);
end

gauss_fwtm = 4.292 * gauss_width; % only valid for circular gaussian
fit_area = 3.1415 * (gauss_fwtm/2.0) * (gauss_fwtm/2.0);

mean_amplitude = 2*pi*gauss_amplitude*gauss_width*gauss_width/fit_area;

noise_std_dev = mean_amplitude ./ snr;

for i = 1:n_graph_points
    
    fprintf('SNR = %8.3f', snr(i));
    
    data = generate_2Dgaussians(data_parameters, n_fits, fit_size);
    data = data + noise_std_dev(i) * randn(fit_size,fit_size,n_fits);
    data = permute(data,[2,1,3]);
    
    %% run GpuFit
    [parameters_GpuFit, converged_GpuFit, chisquare_GpuFit, n_iterations_GpuFit, time_GpuFit]...
        = GpuFit(n_fits, data, model_id, n_parameters, initial_guess_parameters, ...
                 weights, tolerance, max_iterations, parameters_to_fit, estimator_id, user_info);
    
    gpufit_results.a  = parameters_GpuFit(1:n_parameters:end).';
    gpufit_results.x0 = parameters_GpuFit(2:n_parameters:end).';
    gpufit_results.y0 = parameters_GpuFit(3:n_parameters:end).';
    gpufit_results.s  = parameters_GpuFit(4:n_parameters:end).';
    gpufit_results.b  = parameters_GpuFit(5:n_parameters:end).';
    
    converged_GpuFit = converged_GpuFit + 1;
    
    valid_indices = get_valid_fit_results(converged_GpuFit, data_parameters, gpufit_results, chisquare_GpuFit);

    valid_n_iterations = n_iterations_GpuFit(valid_indices);
    
    valid_gpufit_results.a = gpufit_results.a(valid_indices);
    valid_gpufit_results.x0 = gpufit_results.x0(valid_indices);
    valid_gpufit_results.y0 = gpufit_results.y0(valid_indices);
    valid_gpufit_results.s = gpufit_results.s(valid_indices);
    valid_gpufit_results.b = gpufit_results.b(valid_indices);
    
    gpufit_abs_precision.a  = std(valid_gpufit_results.a  - data_parameters.a);
    gpufit_abs_precision.x0 = std(valid_gpufit_results.x0 - data_parameters.x0(valid_indices));
    gpufit_abs_precision.y0 = std(valid_gpufit_results.y0 - data_parameters.y0(valid_indices));
    gpufit_abs_precision.s  = std(valid_gpufit_results.s  - data_parameters.s);
    gpufit_abs_precision.b  = std(valid_gpufit_results.b  - data_parameters.b);
    
    %% save test results
    results_gpufit_mean_iterations(i) = mean(valid_n_iterations);
    results_gpufit_precision(i) = gpufit_abs_precision;
    
    print_fit_info(results_gpufit_precision(i), time_GpuFit, 'Gpufit', numel(valid_indices)/n_fits, results_gpufit_mean_iterations(i));

    %% run cminpack
    [parameters_cminpack, info_cminpack, n_iterations_cminpack, time_cminpack]...
        = cminpack(data, initial_guess_parameters, model_id, tolerance);
    
    cminpack_results.a = parameters_cminpack(1:n_parameters:end).';
    cminpack_results.x0= parameters_cminpack(2:n_parameters:end).';
    cminpack_results.y0 = parameters_cminpack(3:n_parameters:end).';
    cminpack_results.s = parameters_cminpack(4:n_parameters:end).';
    cminpack_results.b = parameters_cminpack(5:n_parameters:end).';
    
    converged_cminpack = (info_cminpack > 0) & (info_cminpack <= 3);
    
    valid_indices = get_valid_fit_results(converged_cminpack, data_parameters, cminpack_results, ones(1,n_fits));
    
    valid_n_iterations = n_iterations_cminpack(valid_indices);
    
    valid_cminpack_results.a = cminpack_results.a(valid_indices);
    valid_cminpack_results.x0 = cminpack_results.x0(valid_indices);
    valid_cminpack_results.y0 = cminpack_results.y0(valid_indices);
    valid_cminpack_results.s = cminpack_results.s(valid_indices);
    valid_cminpack_results.b = cminpack_results.b(valid_indices);

    cminpack_abs_precision.a  = std(valid_cminpack_results.a  - data_parameters.a);
    cminpack_abs_precision.x0 = std(valid_cminpack_results.x0 - data_parameters.x0(valid_indices));
    cminpack_abs_precision.y0 = std(valid_cminpack_results.y0 - data_parameters.y0(valid_indices));
    cminpack_abs_precision.s  = std(valid_cminpack_results.s  - data_parameters.s);
    cminpack_abs_precision.b  = std(valid_cminpack_results.b  - data_parameters.b);
    
    %% save test results
    results_cminpack_mean_iterations(i) = mean(valid_n_iterations);
    results_cminpack_precision(i) = cminpack_abs_precision;
    
    print_fit_info(results_cminpack_precision(i), time_cminpack, 'cminpack', numel(valid_indices)/n_fits, results_cminpack_mean_iterations(i));

end

%% output filename
filename = 'figure5_GpuFit_cminpack';

%% save data
save(filename);

%% write file precision
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