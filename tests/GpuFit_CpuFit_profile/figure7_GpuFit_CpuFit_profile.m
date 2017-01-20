function [] = figure7_GpuFit_CpuFit_profile()


%% number of fits per test point
n_fits = 500000;

%% parameters determining the data to be fit
fit_size = 15;
gauss_amplitude = 500;
gauss_width = 1.0;
gauss_baseline = 10;

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
gauss_pos_offset_max = 2.0;
initial_guess_offset_frac = 0.5;
snr = 10; 
noise = 'gauss';

[data, data_parameters, initial_guess_parameters] = ...
    generate_gauss_fit_test_data(n_fits, fit_size, gauss_amplitude, gauss_width, ...
                                 gauss_baseline, noise, gauss_pos_offset_max, ...
                                 initial_guess_offset_frac, snr);

%% run CpuFit
[parameters_CpuFit,...
 converged_CpuFit,...
 chisquare_CpuFit,...
 n_iterations_CpuFit,...
 times_CpuFit,...
 time_CpuFit]...
    = CpuFit_profile(...
    data,...
    sigma,...
    max_iterations,...
    initial_guess_parameters,...
    parameters_to_fit,...
    model_id,...
    estimator_id,...
    tolerance);

cpufit_results.a  = parameters_CpuFit(1:n_parameters:end).';
cpufit_results.x0 = parameters_CpuFit(2:n_parameters:end).';
cpufit_results.y0 = parameters_CpuFit(3:n_parameters:end).';
cpufit_results.s  = parameters_CpuFit(4:n_parameters:end).';
cpufit_results.b  = parameters_CpuFit(5:n_parameters:end).';

valid_indices = get_valid_fit_results(converged_CpuFit, data_parameters, cpufit_results, chisquare_CpuFit);

valid_n_iterations = n_iterations_CpuFit(valid_indices);

valid_cpufit_results.a  = cpufit_results.a(valid_indices);
valid_cpufit_results.x0 = cpufit_results.x0(valid_indices);
valid_cpufit_results.y0 = cpufit_results.y0(valid_indices);
valid_cpufit_results.s  = cpufit_results.s(valid_indices);
valid_cpufit_results.b  = cpufit_results.b(valid_indices);

cpufit_abs_precision.a  = std(valid_cpufit_results.a  - data_parameters.a);
cpufit_abs_precision.x0 = std(valid_cpufit_results.x0 - data_parameters.x0(valid_indices));
cpufit_abs_precision.y0 = std(valid_cpufit_results.y0 - data_parameters.y0(valid_indices));
cpufit_abs_precision.s  = std(valid_cpufit_results.s  - data_parameters.s);
cpufit_abs_precision.b  = std(valid_cpufit_results.b  - data_parameters.b);

precision_CpuFit = cpufit_abs_precision;
mean_iterations_CpuFit = mean(valid_n_iterations);
print_fit_info(precision_CpuFit, time_CpuFit, 'Cpufit', numel(valid_indices)/n_fits, mean_iterations_CpuFit);

%% run GpuFit
[parameters_GpuFit,...
 converged_GpuFit,...
 chisquare_GpuFit,...
 n_iterations_GpuFit,...
 times_GpuFit,...
 time_GpuFit]...
    = GpuFit_profile(...
    data,...
    sigma,...
    max_iterations,...
    initial_guess_parameters,...
    parameters_to_fit,...
    model_id,...
    estimator_id,...
    tolerance);

gpufit_results.a  = parameters_GpuFit(1:n_parameters:end).';
gpufit_results.x0 = parameters_GpuFit(2:n_parameters:end).';
gpufit_results.y0 = parameters_GpuFit(3:n_parameters:end).';
gpufit_results.s  = parameters_GpuFit(4:n_parameters:end).';
gpufit_results.b  = parameters_GpuFit(5:n_parameters:end).';


valid_indices = get_valid_fit_results(converged_GpuFit, data_parameters, gpufit_results, chisquare_GpuFit);

valid_n_iterations = n_iterations_GpuFit(valid_indices);

valid_gpufit_results.a  = gpufit_results.a(valid_indices);
valid_gpufit_results.x0 = gpufit_results.x0(valid_indices);
valid_gpufit_results.y0 = gpufit_results.y0(valid_indices);
valid_gpufit_results.s  = gpufit_results.s(valid_indices);
valid_gpufit_results.b  = gpufit_results.b(valid_indices);

gpufit_abs_precision.a  = std(valid_gpufit_results.a  - data_parameters.a);
gpufit_abs_precision.x0 = std(valid_gpufit_results.x0 - data_parameters.x0(valid_indices));
gpufit_abs_precision.y0 = std(valid_gpufit_results.y0 - data_parameters.y0(valid_indices));
gpufit_abs_precision.s  = std(valid_gpufit_results.s  - data_parameters.s);
gpufit_abs_precision.b  = std(valid_gpufit_results.b  - data_parameters.b);

precision_GpuFit = gpufit_abs_precision;
mean_iterations_GpuFit = mean(valid_n_iterations);
print_fit_info(precision_GpuFit, time_GpuFit, 'Gpufit', numel(valid_indices)/n_fits, mean_iterations_GpuFit);



%% output filenames
filename = 'figure7_GpuFit_CpuFit_speed';


%% step identifier
step_identifier = {...
    'initialize LM',...
    'allocate GPU memory',...
    'copy data to GPU',...
    'calculate fittting curve',...
    'calculate Chi-square',...
    'calculate gradient & hessian',...
    'Gauss-Jordan elimination',...
    'evaluate fit iteration',...
    'read results from GPU',...
    'get results',...
    'free GPU memory'};

%% write file
xlsfilename = [filename '.xls'];

Raw(1:100, 1:100)=deal(NaN);
xlswrite(xlsfilename,Raw,1)

xlscolumns = {'GpuFit' 'CpuFit'};
xlswrite(xlsfilename,xlscolumns,1,'B1')

xlsrows = step_identifier.';
xlswrite(xlsfilename,xlsrows,1,'A2')

xlsmat(:,1) = times_GpuFit;
xlsmat(:,2) = times_CpuFit;
xlswrite(xlsfilename,xlsmat,1,'B2')

%% plot
Plot_GpuFit_CpuFit_times(times_GpuFit, times_CpuFit, step_identifier);
savefig(filename)

end

function [] = Plot_GpuFit_CpuFit_times(times_GpuFit, times_CpuFit, identifier)

times = [log10(single(times_GpuFit(1:10))+1) ; log10(single(times_CpuFit(1:10))+1)];
b = barh(times.', 1);
set(b(1),'FaceColor','red');
set(b(2),'FaceColor','blue');
legend('GpuFit', 'CpuFit')

set(gca,'yticklabel',identifier)

end