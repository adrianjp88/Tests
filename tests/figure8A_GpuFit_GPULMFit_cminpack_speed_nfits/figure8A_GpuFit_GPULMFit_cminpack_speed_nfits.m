function [] = figure8A_GpuFit_GPULMFit_cminpack_speed_nfits()

%% test parameters
LogNFitsMin = 0;
LogNFitsMax = 6;
sampling_factor = 5;
skip_cminpack = 0;
skip_gpulmfit = 0;

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
initial_guess_offset_frac = 0.4;
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

    %% run GpuFit
    [parameters_GpuFit, converged_GpuFit, chisquare_GpuFit, n_iterations_GpuFit, time_GpuFit]...
        = GpuFit(tmp_n_fits, tmp_data, model_id, tmp_initial_params, weights, tolerance, ...
                 max_iterations, parameters_to_fit, estimator_id, user_info);
    
    converged_GpuFit = converged_GpuFit + 1;
             
    [valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(tmp_data_params, parameters_GpuFit, converged_GpuFit, ...
                                     chisquare_GpuFit, n_iterations_GpuFit);

    %% save test results
    speed_GpuFit(i) = tmp_n_fits/time_GpuFit;
    precision_GpuFit(i) = gpufit_abs_precision.x0;

    print_fit_info(gpufit_abs_precision, time_GpuFit, 'Gpufit', numel(valid_indices)/tmp_n_fits, mean_n_iterations);
   
    %% run GPU-LMFit
    
    if skip_gpulmfit == 0
    
        [parameters_GPULMFit, info_GPULMFit, time_GPULMFit] = GPULMFit(...
            tmp_data,...
            estimator_id,...
            tmp_data_params.s,...
            fit_size);

        converged_GPULMFit = ones(1,tmp_n_fits);
        chisquare_GPULMFit = ones(1,tmp_n_fits);
        n_iterations_GPULMFit = ones(1,tmp_n_fits);
        chk_gpulmfit = 1;
        
        [valid_gpulmfit_results, gpulmfit_abs_precision, mean_n_iterations, valid_indices] = ...
            process_gaussian_fit_results(tmp_data_params, parameters_GPULMFit, converged_GPULMFit, ...
                                         chisquare_GPULMFit, n_iterations_GPULMFit, chk_gpulmfit);
        
        speed_GPULMFit(i) = n_fits(i)/time_GPULMFit;
        precision_GPULMFit(i) = gpulmfit_abs_precision.x0;

        print_fit_info(gpulmfit_abs_precision, time_GPULMFit, 'GPU-LMFit', numel(valid_indices)/tmp_n_fits, mean_n_iterations);
    
    else
       
        speed_GPULMFit(i) = 1.0;
        precision_GPULMFit(i) = 1.0;
        
    end
    
    
    %% run cminpack
    
    if skip_cminpack == 0 
    
        [parameters_cminpack, info_cminpack, n_iterations_cminpack, time_cminpack]...
            = cminpack(tmp_data, tmp_initial_params, model_id, tolerance);
        
        converged_cminpack = (info_cminpack > 0) & (info_cminpack <= 3);
        chisquare_cminpack = ones(1,tmp_n_fits);
        
        [valid_cminpack_results, cminpack_abs_precision, mean_n_iterations, valid_indices] = ...
            process_gaussian_fit_results(tmp_data_params, parameters_cminpack, converged_cminpack, ...
                                         chisquare_cminpack, n_iterations_cminpack);

        %% save test results
        speed_cminpack(i) = tmp_n_fits/time_cminpack;
        precision_cminpack(i) = cminpack_abs_precision.x0;

        print_fit_info(cminpack_abs_precision, time_cminpack, 'C Minpack', numel(valid_indices)/tmp_n_fits, mean_n_iterations);

    else
        
        %% save test results
        speed_cminpack(i) = 1.0;
        precision_cminpack(i) = 1.0;
        
    end

end

%% output filename
filename = 'figure8A_GpuFit_GPULMFit_cminpack_speed_nfits';

%% write file
xlsfilename = [filename '.xls'];

Raw(1:100, 1:100)=deal(NaN);
xlswrite(xlsfilename,Raw,1)

xlscolumns = {'fit count' 'GpuFit' 'GPU-LMFit' 'C Minpack'};
xlswrite(xlsfilename,xlscolumns,1,'A1')

xlsmat(:,1) = n_fits;
xlsmat(:,2) = speed_GpuFit;
xlsmat(:,3) = speed_GPULMFit;
xlsmat(:,4) = speed_cminpack;
xlswrite(xlsfilename,xlsmat,1,'A2')

%% plot
Plot_GpuFit_GPULMFit_cminpack_speed(...
    n_fits,...
    speed_GpuFit,...
    speed_GPULMFit,...
    speed_cminpack)

%% savefig(filename)

end

function [] = Plot_GpuFit_GPULMFit_cminpack_speed(...
    n_fits,...
    speed_GpuFit,...
    speed_GPULMFit,...
    speed_cminpack)

figure('Name','GpuFit vs GPU-LMFit vs Minpack, variable number of fits','NumberTitle','off');
semilogx(...
    n_fits, speed_GpuFit, 'red.-', ...
    n_fits, speed_GPULMFit, 'blue.-', ...
    n_fits, speed_cminpack, 'green.-', ...
    'LineWidth', 8)
xlabel('number of fits')
ylabel('fits per second')
legend('GpuFit', 'GPU-LMFit', 'Minpack')

grid on;
box off;
current_figure = gca;
current_figure.FontSize = 37;
current_figure.LineWidth = 4;
end