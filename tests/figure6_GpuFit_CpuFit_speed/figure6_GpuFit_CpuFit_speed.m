function [] = figure6_gpufit_cpufit_speed()

%% test parameters
LogNFitsMin = 0;
LogNFitsMax = 8;
sampling_factor = 5;
n_timing_repetitions_cpufit = 3;
n_timing_repetitions_gpufit = 2;
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
max_iterations = 20;
model_id = 1; %GAUSS_2D
estimator_id = 0; %LSE
n_parameters = 5;
parameters_to_fit = ones(n_parameters,'int32');
user_info = [];
tolerance = 0.0001;

%% parameters determining the randomness of the data
gauss_pos_offset_max = 0.5;
initial_guess_offset_frac = 0.30;
snr = 10;

%% test setup

n_fits_max = n_fits(end);

[data, data_parameters, initial_guess_parameters] = ...
    generate_gauss_fit_test_data(n_fits_max, fit_size, gauss_amplitude, gauss_width, ...
                                 gauss_baseline, noise, gauss_pos_offset_max, ...
                                 initial_guess_offset_frac, snr);

%% test loop cpufit

for i = 1:length(n_fits)

    if skip_cpufit == 0 
    
        tmp_n_fits = n_fits(i);

        fprintf('%d fits\n', tmp_n_fits);

        tmp_data = data(:,1:tmp_n_fits);
        tmp_initial_params = initial_guess_parameters(:,1:tmp_n_fits);

        tmp_data_params.a  = data_parameters.a;
        tmp_data_params.x0 = data_parameters.x0(1:tmp_n_fits);
        tmp_data_params.y0 = data_parameters.y0(1:tmp_n_fits);
        tmp_data_params.s  = data_parameters.s;
        tmp_data_params.b  = data_parameters.b;
    
        tmp_timings = zeros(1,n_timing_repetitions_cpufit);
        
        for j = 1:n_timing_repetitions_cpufit
        
            %% run cpufit
            [parameters_cpufit, converged_cpufit, chisquare_cpufit, n_iterations_cpufit, time_cpufit]...
                = cpufit(tmp_data, weights, model_id, tmp_initial_params, tolerance, max_iterations, parameters_to_fit, estimator_id, user_info);
            
            tmp_timings(j) = time_cpufit;
            
        end
        
        time_cpufit = mean(tmp_timings);
        time_std_cpufit = std(tmp_timings);
        
        converged_cpufit = converged_cpufit + 1;

        chk_gpulmfit = 0;
        
        [valid_cpufit_results, cpufit_abs_precision, mean_n_iterations, valid_indices] = ...
            process_gaussian_fit_results(tmp_data_params, parameters_cpufit, converged_cpufit, ...
                                         chisquare_cpufit, n_iterations_cpufit, chk_gpulmfit);
        
        %% save test results
        speed_cpufit(i) = tmp_n_fits/time_cpufit;
        speed_std_cpufit(i) = (time_std_cpufit/time_cpufit) * speed_cpufit(i);
        precision_cpufit(i) = cpufit_abs_precision;
        mean_iterations_cpufit(i) = mean_n_iterations;
        print_fit_info(precision_cpufit(i), time_cpufit, 'Cpufit', numel(valid_indices)/tmp_n_fits, mean_n_iterations);

        
    else
        
        %% save test results
        speed_cpufit(i) = 1.0;
        speed_std_cpufit(i) = 1.0;
        precision_cpufit(i) = 1.0;
        mean_iterations_cpufit(i) = 1.0;
        
    end
end

%% test loop gpufit

for i = 1:length(n_fits)

    tmp_n_fits = n_fits(i);
    
    fprintf('%d fits\n', tmp_n_fits);

    tmp_data = data(:,1:tmp_n_fits);
    tmp_initial_params = initial_guess_parameters(:,1:tmp_n_fits);
    
    tmp_data_params.a  = data_parameters.a;
    tmp_data_params.x0 = data_parameters.x0(1:tmp_n_fits);
    tmp_data_params.y0 = data_parameters.y0(1:tmp_n_fits);
    tmp_data_params.s  = data_parameters.s;
    tmp_data_params.b  = data_parameters.b;
    

    tmp_timings = zeros(1,n_timing_repetitions_gpufit);
        
    for j = 1:n_timing_repetitions_gpufit
    
        %% run gpufit
        [parameters_gpufit, converged_gpufit, chisquare_gpufit, n_iterations_gpufit, time_gpufit]...
            = gpufit(tmp_data, weights, model_id, tmp_initial_params, tolerance, ...
                     max_iterations, parameters_to_fit, estimator_id, user_info);
               
        tmp_timings(j) = time_gpufit;
                 
    end
                 
    time_gpufit = mean(tmp_timings);
    time_std_gpufit = std(tmp_timings);
    
    converged_gpufit = converged_gpufit + 1;

    chk_gpulmfit = 0;

    [valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(tmp_data_params, parameters_gpufit, converged_gpufit, ...
                                     chisquare_gpufit, n_iterations_gpufit, chk_gpulmfit);
    
    %% save test results
    speed_gpufit(i) = tmp_n_fits/time_gpufit;
    speed_std_gpufit(i) = (time_std_gpufit/time_gpufit) * speed_gpufit(i);
    precision_gpufit(i) = gpufit_abs_precision;
    mean_iterations_gpufit(i) = mean_n_iterations;

    print_fit_info(precision_gpufit(i), time_gpufit, 'Gpufit', numel(valid_indices)/tmp_n_fits, mean_n_iterations);

    if skip_cpufit == 0 
        speed_increase_factor(i) = speed_gpufit(i)/speed_cpufit(i);
        fprintf('Speedup factor = %f7.2\n', speed_increase_factor(i));
    else
        speed_increase_factor(i) = 1.0;
    end
    
end


%% output filename
filename = 'figure6_gpufit_cpufit_speed';

%% write file
xlsfilename = [filename '.xls'];

Raw(1:100, 1:100)=deal(NaN);
xlswrite(xlsfilename,Raw,1)

xlscolumns = {'n_fits' 'speed_gpufit' 'std_gpufit' 'pct_gpufit' 'speed_cpufit' 'std_cpufit' 'pct_cpufit' 'speedup_factor'};
xlswrite(xlsfilename,xlscolumns,1,'A1')
xlsmat(:,1) = n_fits;
xlsmat(:,2) = speed_gpufit;
xlsmat(:,3) = speed_std_gpufit;
xlsmat(:,4) = speed_std_gpufit./speed_gpufit;
xlsmat(:,5) = speed_cpufit;
xlsmat(:,6) = speed_std_cpufit;
xlsmat(:,7) = speed_std_cpufit./speed_cpufit;
xlsmat(:,8) = speed_increase_factor;
xlswrite(xlsfilename,xlsmat,1,'A2')

%% plot
plot_gpufit_cpufit_speed(n_fits, speed_gpufit, speed_cpufit);

%savefig(filename)

end

function [] = plot_gpufit_cpufit_speed(n_fits, speed_gpufit, speed_cpufit)

while length(speed_gpufit) > length(speed_cpufit)
    speed_cpufit = [speed_cpufit speed_cpufit(end)];
end

figure('Name', 'gpufit vs cpufit', 'NumberTitle', 'off');
semilogx(...
    n_fits, speed_gpufit, 'red.-', ...
    n_fits, speed_cpufit, 'blue.-', ...
    'LineWidth', 8)
xlabel('number of fits')
ylabel('fits per second')
legend('gpufit', 'cpufit')

grid on;
box off;
current_figure = gca;
current_figure.FontSize = 37;
current_figure.LineWidth = 4;

end