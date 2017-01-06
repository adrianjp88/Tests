function [] = GpuFit_CpuFit_profile()

fit_size = 5;
n_fits = 5000000;

sigma = ones(1,fit_size*fit_size);
model_id = 1; %GAUSS_2D
estimator_id = 0; %LSE
n_parameters = 5;
max_iterations = 10;
parameters_to_fit = ones(1,n_parameters);
user_info = 0;
tolerance = 0.001;
snr = 10;
noise = 'gauss';

%% set start values
initial_parameters = ones(1,n_parameters);
initial_parameters(1) = 500;
initial_parameters(2) = (fit_size-1)/2;
initial_parameters(3) = (fit_size-1)/2;
initial_parameters(4) = 1;
initial_parameters(5) = 10;

%% noise
gauss_fwtm = 4.292 * initial_parameters(4); % only valid for circular gaussian
fit_area = gauss_fwtm*gauss_fwtm;
mean_amplitude...
    =2*pi*initial_parameters(1)*initial_parameters(4)*initial_parameters(4)/fit_area;
noise_std_dev = mean_amplitude ./ snr;


fprintf('%d fits\n', n_fits);

%% generate data
xpos_mean = initial_parameters(2);
xpos_offset = rand(n_fits, 1) - 0.5;
input_xpos = xpos_mean + xpos_offset;

ypos_mean = initial_parameters(3);
ypos_offset = rand(n_fits, 1) - 0.5;
input_ypos = ypos_mean + ypos_offset;

parameters.a = initial_parameters(1);
parameters.x0 = input_xpos;
parameters.y0 = input_ypos;
parameters.s = initial_parameters(4);
parameters.b = initial_parameters(5);
data = generate_2Dgaussians(parameters, n_fits, fit_size);
data = data + noise_std_dev * randn(fit_size,fit_size,n_fits);
data = permute(data,[2,1,3]);

initial_parameter_set = repmat(initial_parameters, 1, n_fits);

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
    initial_parameter_set,...
    parameters_to_fit,...
    model_id,...
    estimator_id,...
    tolerance);

converged_CpuFit = converged_CpuFit;
calculated.a = parameters_CpuFit(1:n_parameters:end).';
calculated.x0 = parameters_CpuFit(2:n_parameters:end).';
calculated.y0 = parameters_CpuFit(3:n_parameters:end).';
calculated.s = parameters_CpuFit(4:n_parameters:end).';
calculated.b = parameters_CpuFit(5:n_parameters:end).';

valid_indices = get_valid_fit_indices(converged_CpuFit, chisquare_CpuFit);
valid_n_iterations = n_iterations_CpuFit(valid_indices);
precision_CpuFit = calculate_precision(calculated, parameters, valid_indices);
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
    initial_parameter_set,...
    parameters_to_fit,...
    model_id,...
    estimator_id,...
    tolerance);

converged_GpuFit = converged_GpuFit;
calculated.a = parameters_GpuFit(1:n_parameters:end).';
calculated.x0 = parameters_GpuFit(2:n_parameters:end).';
calculated.y0 = parameters_GpuFit(3:n_parameters:end).';
calculated.s = parameters_GpuFit(4:n_parameters:end).';
calculated.b = parameters_GpuFit(5:n_parameters:end).';

valid_indices = get_valid_fit_indices(converged_GpuFit, chisquare_GpuFit);
valid_n_iterations = n_iterations_GpuFit(valid_indices);
precision_GpuFit = calculate_precision(calculated, parameters, valid_indices);
mean_iterations_GpuFit = mean(valid_n_iterations);
print_fit_info(precision_GpuFit, time_GpuFit, 'Gpufit', numel(valid_indices)/n_fits, mean_iterations_GpuFit);

%% save test info
info.parameters = parameters;
info.initial_parameters = initial_parameter_set;
info.noise = noise;
info.snr = snr;
info.fit_size = fit_size;
info.n_fits = n_fits;
info.model_id = model_id;

%% output filenames
filename = 'GpuFit_CpuFit_speed';

%% save data
save(filename);

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

write_test_info(xlsfilename, info);

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