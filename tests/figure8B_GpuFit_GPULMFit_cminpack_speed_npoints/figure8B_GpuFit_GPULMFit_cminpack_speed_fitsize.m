function [] = gpufit_GPULMFit_cminpack_speed_fitsize()

%% test parameters
fit_size = 5:1:25;
n_fits = 10000;
skip_cminpack = 0;
skip_gpulmfit = 0;

%% parameters determining the data to be fit
gauss_amplitude = 500;
gauss_baseline = 10;
noise = 'gauss';

%% parameters determining how the fit is carried out
weights = [];
max_iterations = 20;
model_id = 1; %GAUSS_2D
estimator_id = 0; %LSE
n_parameters = 5;
parameters_to_fit = ones(n_parameters,1);
user_info = [];
tolerance = 0.0001;

%% parameters determining the randomness of the data
gauss_pos_offset_max = 0.5;
initial_guess_offset_frac = 0.1;
snr = 10;


%% test loop
for i = 1:length(fit_size)

    fprintf('fit size: %d\n', fit_size(i));
    
    %%set weight
    sigma = ones(1,fit_size(i)*fit_size(i));
    
    %% set width
    gauss_width = 1/7*fit_size(i);
    
    [data, data_parameters, initial_guess_parameters] = ...
        generate_gauss_fit_test_data(n_fits, fit_size(i), gauss_amplitude, gauss_width, ...
                                     gauss_baseline, noise, gauss_pos_offset_max, ...
                                     initial_guess_offset_frac, snr);

    %% run gpufit
    [parameters_gpufit, converged_gpufit, chisquare_gpufit, n_iterations_gpufit, time_gpufit]...
        = gpufit(data, weights, model_id, initial_guess_parameters, tolerance, ...
                 max_iterations, parameters_to_fit, estimator_id, user_info);
             
    converged_gpufit = converged_gpufit + 1;
             
    [valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(data_parameters, parameters_gpufit, converged_gpufit, ...
                                     chisquare_gpufit, n_iterations_gpufit);
    
    speed_gpufit(i) = n_fits/time_gpufit;
    print_fit_info(gpufit_abs_precision, time_gpufit, 'Gpufit', numel(valid_indices)/n_fits, mean_n_iterations);
   
    
    %% run GPU-LMFit
   
    if skip_gpulmfit == 0
    
        [parameters_GPULMFit, info_GPULMFit, time_GPULMFit] = GPULMFit(...
            data,...
            estimator_id,...
            data_parameters.s,...
            fit_size(i));

        converged_GPULMFit = ones(1,n_fits);
        chisquare_GPULMFit = ones(1,n_fits);
        n_iterations_GPULMFit = ones(1,n_fits);
        chk_gpulmfit = 1;
        
        [valid_gpulmfit_results, gpulmfit_abs_precision, mean_n_iterations, valid_indices] = ...
            process_gaussian_fit_results(data_parameters, parameters_GPULMFit, converged_GPULMFit, ...
                                         chisquare_GPULMFit, n_iterations_GPULMFit, chk_gpulmfit);

        speed_GPULMFit(i) = n_fits/time_GPULMFit;
        print_fit_info(gpulmfit_abs_precision, time_GPULMFit, 'GPU-LMFit', numel(valid_indices)/n_fits, mean_n_iterations);

    else
        
        speed_GPULMFit(i) = 1.0;
        
    end
    
    %% run cminpack
    
    if skip_cminpack == 0
    
        [parameters_cminpack, info_cminpack, n_iterations_cminpack, time_cminpack]...
            = cminpack(data, initial_guess_parameters, model_id, tolerance);

        converged_cminpack = (info_cminpack > 0) & (info_cminpack <= 3);
        chisquare_cminpack = ones(1,n_fits);

        [valid_cminpack_results, cminpack_abs_precision, mean_n_iterations, valid_indices] = ...
            process_gaussian_fit_results(data_parameters, parameters_cminpack, converged_cminpack, ...
                                         chisquare_cminpack, n_iterations_cminpack);

        %% save test results
        speed_cminpack(i) = n_fits/time_cminpack;
        print_fit_info(cminpack_abs_precision, time_cminpack, 'C Minpack', numel(valid_indices)/n_fits, mean_n_iterations);

    else
        
        speed_cminpack(i) = 1.0;
        
    end
        
end

%% output filenames
filename = 'figure8B_gpufit_GPULMFit_cminpack_speed_fitsize';

%% write file
xlsfilename = [filename '.xls'];

Raw(1:100, 1:100)=deal(NaN);
xlswrite(xlsfilename,Raw,1)

xlscolumns = {'fit size' 'gpufit' 'GPU_LMFit' 'C Minpack'};
xlswrite(xlsfilename,xlscolumns,1,'A1')

xlsmat(:,1) = fit_size;
xlsmat(:,2) = speed_gpufit;
xlsmat(:,3) = speed_GPULMFit;
xlsmat(:,4) = speed_cminpack;
xlswrite(xlsfilename,xlsmat,1,'A2')

%% plot
Plot_gpufit_GPULMFit_Minpack_speed(...
    fit_size,...
    speed_gpufit,...
    speed_GPULMFit,...
    speed_cminpack)

savefig(filename)

end

function [] = Plot_gpufit_GPULMFit_Minpack_speed(...
    fit_size,...
    speed_gpufit,...
    speed_GPULMFit,...
    speed_cminpack)

figure('Name','gpufit vs GPU-LMFit vs Minpack, variable fit size','NumberTitle','off');
plot(...
    fit_size, speed_gpufit, 'red.-', ...
    fit_size, speed_GPULMFit, 'blue.-', ...
    fit_size, speed_cminpack, 'green.-', ...
    'LineWidth', 8)
xlabel('fit size')
ylabel('fits per second')
legend('gpufit', 'GPU-LMFit', 'Minpack')

grid on;
box off;
current_figure = gca;
current_figure.FontSize = 37;
current_figure.LineWidth = 4;
end