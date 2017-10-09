function [] = Gpufit_speed_gpu_memory()

n_timing_repetitions_gpufit = 10;

%% gpu configuration

gpu_memory_fraction = 0.05:0.05:0.9;
n_fits_per_block = [];
min_threads = [];
max_threads = [];

%% parameters determining the data to be fit
n_fits = 10;
fit_size = 32;
gauss_amplitude = 500;
gauss_width = 1/7*fit_size;
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
gauss_pos_offset_max = 0.5;
initial_guess_offset_frac = 0.3;
snr = 60;

%% test setup

[data, data_parameters, initial_guess_parameters] = ...
    generate_gauss_fit_test_data(n_fits, fit_size, gauss_amplitude, gauss_width, ...
                                 gauss_baseline, noise, gauss_pos_offset_max, ...
                                 initial_guess_offset_frac, snr);

%% test loop
for i = 1:length(gpu_memory_fraction)

    fprintf('%d fits\n', n_fits);

    tmp_data = data(:,1:n_fits);
    tmp_initial_params = initial_guess_parameters(:,1:n_fits);
    
    tmp_data_params.a  = data_parameters.a;
    tmp_data_params.x0 = data_parameters.x0(1:n_fits);
    tmp_data_params.y0 = data_parameters.y0(1:n_fits);
    tmp_data_params.s  = data_parameters.s;
    tmp_data_params.b  = data_parameters.b;
        
    %% run gpufit
    for repetition = 1:n_timing_repetitions_gpufit
        [parameters_gpufit, converged_gpufit, chisquare_gpufit, n_iterations_gpufit, time_gpufit, info_gpufit(i,:)]...
            = gpufit(tmp_data, weights, model_id, tmp_initial_params, tolerance, ...
                     max_iterations, parameters_to_fit, estimator_id, user_info,...
                     gpu_memory_fraction(i),n_fits_per_block, min_threads, max_threads);
        tmp_timings(repetition) = time_gpufit;
    end
    
    time_gpufit = mean(tmp_timings);
    converged_gpufit = converged_gpufit + 1;
             
    [valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(tmp_data_params, parameters_gpufit, converged_gpufit, ...
                                     chisquare_gpufit, n_iterations_gpufit);

    %% save test results
    speed_gpufit(i) = n_fits/time_gpufit;  
    precision_gpufit(i) = gpufit_abs_precision.x0;

    print_fit_info(gpufit_abs_precision, time_gpufit, 'Gpufit', numel(valid_indices)/n_fits, mean_n_iterations);
end
%% test results
test.output = table(...
    gpu_memory_fraction',...
    speed_gpufit',...
    info_gpufit(:,1),...
    info_gpufit(:,2),...
    info_gpufit(:,3),...
    info_gpufit(:,4),...
    info_gpufit(:,5),...
    info_gpufit(:,6));
test.output.Properties.VariableNames = {...
    'gpu_memory_fraction',...
    'fits_per_second',...
    'n_chunks',...
    'min_chunk_size',...
    'max_chunk_size',...
    'min_fits_per_block',...
    'max_fits_per_block',...
    'blocks_per_fit'};
test.params.n_fits = n_fits;
test.params.fit_size = fit_size;
test.params.gpu_memory_fraction = gpu_memory_fraction;
test.params.min_threads = min_threads;
test.params.max_threads = max_threads;
test.params.n_fits_per_block = n_fits_per_block;
test1 = test;
if exist('test_memory.mat', 'file') == 2
    save('test_memory.mat', 'test1', '-append');
else
    save('test_memory.mat', 'test1');
end

%% plot
Plot_gpufit_speed(gpu_memory_fraction, speed_gpufit)

end

function [] = Plot_gpufit_speed(gpu_memory_fraction, speed_gpufit)

figure('Name','gpufit, variable used gpu memory fraction','NumberTitle','off');
plot(...
    gpu_memory_fraction, speed_gpufit, 'red.-', 'LineWidth', 8)
xlabel('used gpu memory fraction')
ylabel('fits per second')
legend('gpufit')

grid on;
box off;
current_figure = gca;
current_figure.FontSize = 37;
current_figure.LineWidth = 4;
end