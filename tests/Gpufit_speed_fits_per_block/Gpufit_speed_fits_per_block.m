function [] = Gpufit_speed_fits_per_block()

n_timing_repetitions_gpufit = 10;

%% gpu configuration

gpu_memory_fraction = 0.9;
min_threads = [];
max_threads = [];

%% parameters determining the data to be fit
n_fits = 300000;
fit_size = 5:1:22;
gauss_amplitude = 500;
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

info_gpufit = zeros(length(fit_size),32,6);
speed_gpufit = zeros(length(fit_size),32);
n_fits_per_block = 1:32;
hist = zeros(18,1024);

%% test loop
for j = 1:length(fit_size)
    
    n_points_pow2 = 2^(ceil(log(fit_size(j)^2)/log(2)));
    max_fits_per_block = 1024 / n_points_pow2;
    
    gauss_width = 1/7*fit_size(j);
    [data, data_parameters, initial_guess_parameters] = ...
    generate_gauss_fit_test_data(n_fits, fit_size(j), gauss_amplitude, gauss_width, ...
                                 gauss_baseline, noise, gauss_pos_offset_max, ...
                                 initial_guess_offset_frac, snr);
                             
    for i = 1:max_fits_per_block

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
            [parameters_gpufit, converged_gpufit, chisquare_gpufit, n_iterations_gpufit, time_gpufit, info_gpufit(j,i,:)]...
                = gpufit(tmp_data, weights, model_id, tmp_initial_params, tolerance, ...
                         max_iterations, parameters_to_fit, estimator_id, user_info,...
                         gpu_memory_fraction,n_fits_per_block(i), min_threads, max_threads);
            tmp_timings(repetition) = time_gpufit;
        end

        time_gpufit = mean(tmp_timings);
        converged_gpufit = converged_gpufit + 1;

        [valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
            process_gaussian_fit_results(tmp_data_params, parameters_gpufit, converged_gpufit, ...
                                         chisquare_gpufit, n_iterations_gpufit);

        %% save test results
        threads_per_block(j,i) = fit_size(j)^2 * n_fits_per_block(i);
        speed_gpufit(j,i) = n_fits * fit_size(j)^2 / time_gpufit;  
        precision_gpufit(i) = gpufit_abs_precision.x0;

        print_fit_info(gpufit_abs_precision, time_gpufit, 'Gpufit', numel(valid_indices)/n_fits, mean_n_iterations);
    end
    
    test(j).output = table(...
        repmat(fit_size(j), length(n_fits_per_block), 1),...
        n_fits_per_block',...
        n_fits_per_block' * n_points_pow2,...
        speed_gpufit(j,:)',...
        info_gpufit(j,:,1)',...
        info_gpufit(j,:,2)',...
        info_gpufit(j,:,3)',...
        info_gpufit(j,:,4)',...
        info_gpufit(j,:,5)',...
        info_gpufit(j,:,6)');
    test(j).output.Properties.VariableNames = {...
        'fit_size',...
        'n_fits_per_block',...
        'threads_per_block'...
        'data_points_per_second',...
        'n_chunks',...
        'min_chunk_size',...
        'max_chunk_size',...
        'min_fits_per_block',...
        'max_fits_per_block',...
        'blocks_per_fit'};
    test(j).params.n_fits = n_fits;
    test(j).params.fit_size = fit_size(j);
    test(j).params.gpu_memory_fraction = gpu_memory_fraction;
    test(j).params.min_threads = min_threads;
    test(j).params.max_threads = max_threads;
    test(j).params.n_fits_per_block = n_fits_per_block;
    
    hist(j,threads_per_block(j,1:max_fits_per_block)) = speed_gpufit(j,1:max_fits_per_block);
    
end
%% test results
save('test_fits_per_block.mat', 'test');

%% plot
Plot_gpufit_speed(hist);

end

function [] = Plot_gpufit_speed(threads_per_block)


figure('Name','gpufit, variable number of fits per block','NumberTitle','off');
imagesc(threads_per_block);
xlabel('number of threads per block')
ylabel('fit size')
colormap hot;
colorbar;

grid on;
box off;
current_figure = gca;
current_figure.FontSize = 30;
current_figure.LineWidth = 4;
end