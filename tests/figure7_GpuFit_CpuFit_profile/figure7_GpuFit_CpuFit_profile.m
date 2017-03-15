function [] = figure7_GpuFit_CpuFit_profile()

%% number of fits per test point
n_fits = 5000000;

%% parameters determining the data to be fit
fit_size = 15;
gauss_amplitude = 500;
gauss_width = 1.0;
gauss_baseline = 10;

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
gauss_pos_offset_max = 1.0;
initial_guess_offset_frac = 0.3;
snr = 60; 
noise = 'gauss';

[data, data_parameters, initial_guess_parameters] = ...
    generate_gauss_fit_test_data(n_fits, fit_size, gauss_amplitude, gauss_width, ...
                                 gauss_baseline, noise, gauss_pos_offset_max, ...
                                 initial_guess_offset_frac, snr);

%% run CpuFit
[parameters_cpufit,...
 converged_cpufit,...
 chisquare_cpufit,...
 n_iterations_cpufit,...
 times_cpufit,...
 time_cpufit]...
    = cpufit_profile(...
    data,...
    weights,...
    model_id,...
    initial_guess_parameters,...
    tolerance,...
    max_iterations,...
    parameters_to_fit,...
    estimator_id,...
    user_info);

converged_cpufit = converged_cpufit + 1;

times_cpufit = double(times_cpufit);
cpu_times_frac = 100*times_cpufit/sum(times_cpufit);

chk_gpulmfit = 0;
[valid_cpufit_results, cpufit_abs_precision, mean_n_iterations, valid_indices] = ...
    process_gaussian_fit_results(data_parameters, parameters_cpufit, converged_cpufit, ...
                                 chisquare_cpufit, n_iterations_cpufit, chk_gpulmfit);

precision_CpuFit = cpufit_abs_precision;
mean_iterations_CpuFit = mean_n_iterations;
print_fit_info(precision_CpuFit, time_cpufit, 'Cpufit', numel(valid_indices)/n_fits, mean_n_iterations);

%% run GpuFit
[parameters_gpufit,...
 converged_gpufit,...
 chisquare_gpufit,...
 n_iterations_gpufit,...
 times_gpufit,...
 time_gpufit]...
    = gpufit_profile(...
    data,...
    weights,...
    model_id,...
    initial_guess_parameters,...
    tolerance,...
    max_iterations,...
    parameters_to_fit,...
    estimator_id,...
    user_info);

converged_gpufit = converged_gpufit + 1;

times_gpufit = double(times_gpufit);
gpu_times_frac = 100*times_gpufit/sum(times_gpufit);

[valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
    process_gaussian_fit_results(data_parameters, parameters_gpufit, converged_gpufit, ...
                                 chisquare_gpufit, n_iterations_gpufit);


precision_GpuFit = gpufit_abs_precision;
mean_iterations_GpuFit = mean_n_iterations;
print_fit_info(precision_GpuFit, time_gpufit, 'Gpufit', numel(valid_indices)/n_fits, mean_n_iterations);



%% output filenames
filename = 'figure7_GpuFit_CpuFit_speed';


%% step identifier
step_identifier = {...
    'initialize LM',...
    'allocate GPU memory',...
    'copy data to GPU',...
    'fittting curve',...
    'Chi-square',...
    'gradient & hessian',...
    'Gauss-Jordan elimination',...
    'evaluate iteration',...
    'read results from GPU',...
    'free GPU memory'};

% write file
xlsfilename = [filename '.xls'];

Raw(1:100, 1:100)=deal(NaN);
xlswrite(xlsfilename,Raw,1)

xlscolumns = {'GpuFit' '%' 'CpuFit' '%'};
xlswrite(xlsfilename,xlscolumns,1,'B1')

xlsrows = step_identifier.';
xlswrite(xlsfilename,xlsrows,1,'A2')

xlsmat(:,1) = times_gpufit;
xlsmat(:,2) = gpu_times_frac;
xlsmat(:,3) = times_cpufit;
xlsmat(:,4) = cpu_times_frac;
xlswrite(xlsfilename,xlsmat,1,'B2')

%% plot
Plot_GpuFit_CpuFit_times(times_gpufit, times_cpufit, step_identifier);

end

function [] = Plot_GpuFit_CpuFit_times(times_gpufit, times_cpufit, identifier)
%% Cpufit & GPUfit
times = [times_gpufit; times_cpufit];
b = barh(times.');
set(b(1),'FaceColor','red');
set(b(2),'FaceColor','blue');
legend('GpuFit', 'CpuFit')
set(gca,'yticklabel',identifier)
current_figure = gca;
current_figure.FontSize = 20;
current_figure.LineWidth = 4;
current_figure.XScale = 'log';
savefig('figure_7')

figure;
%% CPUfit
times = [times_cpufit; zeros(1,10, 'double')];
cpu_times_frac = 100*times_cpufit/sum(times_cpufit);
cpufit_legend{1} = ['initialize LM: ', num2str(cpu_times_frac(1)), ' %'];
cpufit_legend{2} = '--';
cpufit_legend{3} = '--';
cpufit_legend{4} = ['fittting curve : ', num2str(cpu_times_frac(4)), ' %'];
cpufit_legend{5} = ['Chi-square : ', num2str(cpu_times_frac(5)), ' %'];
cpufit_legend{6} = ['gradient & hessian : ', num2str(cpu_times_frac(6)), ' %'];
cpufit_legend{7} = ['Gauss-Jordan elimination : ', num2str(cpu_times_frac(7)), ' %'];
cpufit_legend{8} = ['evaluate iteration : ', num2str(cpu_times_frac(8)), ' %'];
cpufit_legend{9} = '--';
cpufit_legend{10} = '--';

subplot(2,1,1)
barh(times, 'stacked');
set(gca,'yticklabel',{'CPUfit',''})
legend(cpufit_legend)
current_figure = gca;
current_figure.FontSize = 20;
current_figure.LineWidth = 4;

%% GPUfit
times = [times_gpufit; zeros(1,10, 'double')];
gpu_times_frac = 100*times_gpufit/sum(times_gpufit);
gpufit_legend{1} = ['initialize LM: ', num2str(gpu_times_frac(1)), ' %'];
gpufit_legend{2} = ['allocate GPU memory :', num2str(gpu_times_frac(2)), ' %'];
gpufit_legend{3} = ['copy data to GPU : ', num2str(gpu_times_frac(3)), ' %'];
gpufit_legend{4} = ['fittting curve : ', num2str(gpu_times_frac(4)), ' %'];
gpufit_legend{5} = ['Chi-square : ', num2str(gpu_times_frac(5)), ' %'];
gpufit_legend{6} = ['gradient & hessian : ', num2str(gpu_times_frac(6)), ' %'];
gpufit_legend{7} = ['Gauss-Jordan elimination : ', num2str(gpu_times_frac(7)), ' %'];
gpufit_legend{8} = ['evaluate iteration : ', num2str(gpu_times_frac(8)), ' %'];
gpufit_legend{9} = ['read results from GPU : ', num2str(gpu_times_frac(9)), ' %'];
gpufit_legend{10} = ['free GPU memory : ', num2str(gpu_times_frac(10)), ' %'];

subplot(2,1,2)
barh(times, 'stacked');
set(gca,'yticklabel',{'GPUfit',''})
legend(gpufit_legend)
current_figure = gca;
current_figure.FontSize = 20;
current_figure.LineWidth = 4;
savefig('figure_1_and_3')

end