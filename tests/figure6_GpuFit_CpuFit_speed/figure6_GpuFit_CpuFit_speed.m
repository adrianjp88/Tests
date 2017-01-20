function [] = figure6_GpuFit_cpufit_speed()

%% test parameters
LogNFitsMin = 0;
LogNFitsMax = 6;
sampling_factor = 5;
skip_cpufit = 1;

%% set up n_fits parameter
ranges = logspace(LogNFitsMin,LogNFitsMax,LogNFitsMax-LogNFitsMin+1);
temp = zeros(LogNFitsMax-LogNFitsMin,10/sampling_factor);
stepslog = LogNFitsMin:LogNFitsMax;
for index = 1:length(stepslog)-1
    steps = 10^(stepslog(index)) * sampling_factor;
    temp(index,:) = steps:steps:ranges(index+1);
end
n_fits = reshape(temp', [1 10/sampling_factor*(LogNFitsMax-LogNFitsMin)]); 
n_fits = [10^LogNFitsMin n_fits];

%% parameters determining the data to be fit
fit_size = 5;
gauss_amplitude = 500;
gauss_width = 1.0;
gauss_baseline = 10;
noise = 'gauss';

%% parameters determining how the fit is carried out
weights = [];
sigma = ones(1,fit_size*fit_size);
max_iterations = 20;
model_id = 1; %GAUSS_2D
estimator_id = 0; %LSE
n_parameters = 5;
parameters_to_fit = ones(1,n_parameters);
user_info = 0;
tolerance = 0.0001;

%% parameters determining the randomness of the data
gauss_pos_offset_max = 0.5;
initial_guess_offset_frac = 0.1;
snr = 10;

%% test setup

n_fits_max = n_fits(end);

[data, data_parameters, initial_guess_parameters] = ...
    generate_gauss_fit_test_data(n_fits_max, fit_size, gauss_amplitude, gauss_width, ...
                                 gauss_baseline, noise, gauss_pos_offset_max, ...
                                 initial_guess_offset_frac, snr);

%% test loop
for i = 1:length(n_fits)

    tmp_n_fits = n_fits(i);
    
    fprintf('%d fits\n', tmp_n_fits);

    tmp_data = data(:,:,1:tmp_n_fits);
    tmp_initial_params = initial_guess_parameters(1:tmp_n_fits*n_parameters);
    
    tmp_data_params.a  = data_parameters.a;
    tmp_data_params.x0 = data_parameters.x0(1:tmp_n_fits);
    tmp_data_params.y0 = data_parameters.y0(1:tmp_n_fits);
    tmp_data_params.s  = data_parameters.s;
    tmp_data_params.b  = data_parameters.b;
    
    if skip_cpufit == 0 
    
        %% run cpufit
        [parameters_cpufit, converged_cpufit, chisquare_cpufit, n_iterations_cpufit, time_cpufit]...
            = CpuFit(tmp_data, sigma, max_iterations, tmp_initial_params, parameters_to_fit, model_id, estimator_id, tolerance, user_info);

        converged_cpufit = converged_cpufit + 1;

        chk_gpulmfit = 0;
        
        [valid_cpufit_results, cpufit_abs_precision, mean_n_iterations, valid_indices] = ...
            process_gaussian_fit_results(tmp_data_params, parameters_cpufit, converged_cpufit, ...
                                         chisquare_cpufit, n_iterations_cpufit, chk_gpulmfit);
        
        %% save test results
        speed_cpufit(i) = tmp_n_fits/time_cpufit;
        precision_cpufit(i) = cpufit_abs_precision;
        mean_iterations_cpufit(i) = mean_n_iterations;
        print_fit_info(precision_cpufit(i), time_cpufit, 'Cpufit', numel(valid_indices)/tmp_n_fits, mean_n_iterations);

        
    else
        
        %% save test results
        speed_cpufit(i) = 1.0;
        precision_cpufit(i) = 1.0;
        mean_iterations_cpufit(i) = 1.0;
        
    end
        
      
    %% run GpuFit
    [parameters_GpuFit, converged_GpuFit, chisquare_GpuFit, n_iterations_GpuFit, time_GpuFit]...
        = GpuFit(tmp_n_fits, tmp_data, model_id, tmp_initial_params, weights, tolerance, ...
                 max_iterations, parameters_to_fit, estimator_id, user_info);
    
    converged_GpuFit = converged_GpuFit + 1;

    chk_gpulmfit = 0;

    [valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(tmp_data_params, parameters_GpuFit, converged_GpuFit, ...
                                     chisquare_GpuFit, n_iterations_GpuFit, chk_gpulmfit);
    
    %% save test results
    speed_GpuFit(i) = tmp_n_fits/time_GpuFit;
    precision_GpuFit(i) = gpufit_abs_precision;
    mean_iterations_GpuFit(i) = mean_n_iterations;

    print_fit_info(precision_GpuFit(i), time_GpuFit, 'Gpufit', numel(valid_indices)/tmp_n_fits, mean_n_iterations);

    if skip_cpufit == 0 
        speed_increase_factor(i) = speed_GpuFit(i)/speed_cpufit(i);
        fprintf('Speedup factor = %f7.2\n', speed_increase_factor(i));
    else
        speed_increase_factor(i) = 1.0;
    end
    
end

%% output filename
filename = 'figure6_GpuFit_CpuFit_speed';

%% write file
xlsfilename = [filename '.xls'];

Raw(1:100, 1:100)=deal(NaN);
xlswrite(xlsfilename,Raw,1)

xlscolumns = {'number of fits' 'GpuFit' 'cpufit' 'speedup_factor'};
xlswrite(xlsfilename,xlscolumns,1,'A1')
xlsmat(:,1) = n_fits;
xlsmat(:,2) = speed_GpuFit;
xlsmat(:,3) = speed_cpufit;
xlsmat(:,4) = speed_increase_factor;
xlswrite(xlsfilename,xlsmat,1,'A2')

%% plot
Plot_GpuFit_cpufit_speed(n_fits, speed_GpuFit, speed_cpufit);
savefig(filename)

end

function [] = Plot_GpuFit_cpufit_speed(n_fits, speed_GpuFit, speed_cpufit)

while length(speed_GpuFit) > length(speed_cpufit)
    speed_cpufit = [speed_cpufit speed_cpufit(end)];
end

figure('Name', 'GpuFit vs cpufit', 'NumberTitle', 'off');
semilogx(...
    n_fits, speed_GpuFit, 'red.-', ...
    n_fits, speed_cpufit, 'blue.-', ...
    'LineWidth', 8)
xlabel('number of fits')
ylabel('fits per second')
legend('GpuFit', 'cpufit')

grid on;
box off;
current_figure = gca;
current_figure.FontSize = 37;
current_figure.LineWidth = 4;

end