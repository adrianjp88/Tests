function [] = figure6_GpuFit_cpufit_speed()

%% test parameters
LogNFitsMin = 0;
LogNFitsMax = 6;
sampling_factor = 5;
skip_cpufit = 0;

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

%% generate the data for the largest number of fits, then use subsets
%% of this data in the test loop

n_fits_max = n_fits(end);

gauss_xpos_mean = (fit_size-1)/2;
gauss_ypos_mean = (fit_size-1)/2;

gauss_pos_x = gauss_xpos_mean + ( 2.0*gauss_pos_offset_max*(rand(n_fits_max, 1) - 0.5) );
gauss_pos_y = gauss_ypos_mean + ( 2.0*gauss_pos_offset_max*(rand(n_fits_max, 1) - 0.5) );

initial_guess_xpos_offset_max = initial_guess_offset_frac * gauss_width;
initial_guess_ypos_offset_max = initial_guess_offset_frac * gauss_width;
initial_guess_ampl_offset_max = initial_guess_offset_frac * gauss_amplitude;
initial_guess_width_offset_max = initial_guess_offset_frac * gauss_width;
initial_guess_baseline_offset_max = initial_guess_offset_frac * gauss_baseline;

initial_guess_xpos_offset = initial_guess_xpos_offset_max ...
                            * 2.0 * (rand(n_fits_max, 1) - 0.5);
                  
initial_guess_ypos_offset = initial_guess_ypos_offset_max ...
                            * 2.0 * (rand(n_fits_max, 1) - 0.5);

initial_guess_ampl_offset = initial_guess_ampl_offset_max ...
                            * 2.0 * (rand(n_fits_max, 1) - 0.5);
                        
initial_guess_width_offset = initial_guess_width_offset_max ...
                             * 2.0 * (rand(n_fits_max, 1) - 0.5);
                        
initial_guess_baseline_offset = initial_guess_baseline_offset_max ...
                                * 2.0 * (rand(n_fits_max, 1) - 0.5);
                        
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

initial_guess_parameters = ones(1,n_parameters*n_fits_max);

for i = 1:n_fits_max
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

%% generate the data

data = generate_2Dgaussians(data_parameters, n_fits_max, fit_size);
data = data + noise_std_dev * randn(fit_size,fit_size,n_fits_max);
data = permute(data,[2,1,3]);

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

        cpufit_results.a  = parameters_cpufit(1:n_parameters:end).';
        cpufit_results.x0 = parameters_cpufit(2:n_parameters:end).';
        cpufit_results.y0 = parameters_cpufit(3:n_parameters:end).';
        cpufit_results.s  = parameters_cpufit(4:n_parameters:end).';
        cpufit_results.b  = parameters_cpufit(5:n_parameters:end).';
    
        converged_cpufit = converged_cpufit + 1;

        valid_indices = get_valid_fit_results(converged_cpufit, tmp_data_params, cpufit_results, chisquare_cpufit);

        valid_n_iterations = n_iterations_cpufit(valid_indices);

        valid_cpufit_results.a = cpufit_results.a(valid_indices);
        valid_cpufit_results.x0 = cpufit_results.x0(valid_indices);
        valid_cpufit_results.y0 = cpufit_results.y0(valid_indices);
        valid_cpufit_results.s = cpufit_results.s(valid_indices);
        valid_cpufit_results.b = cpufit_results.b(valid_indices);

        cpufit_abs_precision.a  = std(valid_cpufit_results.a  - tmp_data_params.a);
        cpufit_abs_precision.x0 = std(valid_cpufit_results.x0 - tmp_data_params.x0(valid_indices));
        cpufit_abs_precision.y0 = std(valid_cpufit_results.y0 - tmp_data_params.y0(valid_indices));
        cpufit_abs_precision.s  = std(valid_cpufit_results.s  - tmp_data_params.s);
        cpufit_abs_precision.b  = std(valid_cpufit_results.b  - tmp_data_params.b);

        %% save test results
        speed_cpufit(i) = tmp_n_fits/time_cpufit;
        precision_cpufit(i) = cpufit_abs_precision;
        mean_iterations_cpufit(i) = mean(valid_n_iterations);
        print_fit_info(precision_cpufit(i), time_cpufit, 'Cpufit', numel(valid_indices)/tmp_n_fits, mean_iterations_cpufit(i));

        
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
    
    gpufit_results.a  = parameters_GpuFit(1:n_parameters:end).';
    gpufit_results.x0 = parameters_GpuFit(2:n_parameters:end).';
    gpufit_results.y0 = parameters_GpuFit(3:n_parameters:end).';
    gpufit_results.s  = parameters_GpuFit(4:n_parameters:end).';
    gpufit_results.b  = parameters_GpuFit(5:n_parameters:end).';
    
    converged_GpuFit = converged_GpuFit + 1;
    
    valid_indices = get_valid_fit_results(converged_GpuFit, tmp_data_params, gpufit_results, chisquare_GpuFit);

    valid_n_iterations = n_iterations_GpuFit(valid_indices);
    
    valid_gpufit_results.a = gpufit_results.a(valid_indices);
    valid_gpufit_results.x0 = gpufit_results.x0(valid_indices);
    valid_gpufit_results.y0 = gpufit_results.y0(valid_indices);
    valid_gpufit_results.s = gpufit_results.s(valid_indices);
    valid_gpufit_results.b = gpufit_results.b(valid_indices);
    
    gpufit_abs_precision.a  = std(valid_gpufit_results.a  - tmp_data_params.a);
    gpufit_abs_precision.x0 = std(valid_gpufit_results.x0 - tmp_data_params.x0(valid_indices));
    gpufit_abs_precision.y0 = std(valid_gpufit_results.y0 - tmp_data_params.y0(valid_indices));
    gpufit_abs_precision.s  = std(valid_gpufit_results.s  - tmp_data_params.s);
    gpufit_abs_precision.b  = std(valid_gpufit_results.b  - tmp_data_params.b);

    %% save test results
    speed_GpuFit(i) = tmp_n_fits/time_GpuFit;
    precision_GpuFit(i) = gpufit_abs_precision;
    mean_iterations_GpuFit(i) = mean(valid_n_iterations);

    print_fit_info(precision_GpuFit(i), time_GpuFit, 'Gpufit', numel(valid_indices)/tmp_n_fits, mean_iterations_GpuFit(i));

    if skip_cpufit == 0 
        speed_increase_factor(i) = speed_GpuFit(i)/speed_cpufit(i);
        fprintf('Speedup factor = %f7.2\n', speed_increase_factor(i));
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