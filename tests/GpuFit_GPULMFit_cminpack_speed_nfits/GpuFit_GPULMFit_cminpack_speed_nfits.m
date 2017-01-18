function [] = GpuFit_GPULMFit_cminpack_speed_nfits()

fit_size = 16;
sigma = [];
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

%% test case
LogNFitsMin = 0;
LogNFitsMax = 5;
sampling_factor = 5;
ranges = logspace(LogNFitsMin,LogNFitsMax,LogNFitsMax-LogNFitsMin+1);
temp = zeros(LogNFitsMax-LogNFitsMin,10/sampling_factor);
stepslog = LogNFitsMin:LogNFitsMax;
for index = 1:length(stepslog)-1
    steps = 10^(stepslog(index)) * sampling_factor;
    temp(index,:) = steps:steps:ranges(index+1);
end
n_fits = reshape(temp', [1 10/sampling_factor*(LogNFitsMax-LogNFitsMin)]); 
n_fits = [10^LogNFitsMin n_fits];

%% noise
mean_amplitude = 2 * pi * initial_parameters(1) * initial_parameters(4) * initial_parameters(4)/fit_size/fit_size;
noise_std_dev = mean_amplitude / snr;

%% test loop
for i = 1:length(n_fits)

    fprintf('%d fits\n', n_fits(i));

    %% generate data
    xpos_mean = initial_parameters(2);
    xpos_offset = rand(n_fits(i), 1) - 0.5;
    input_xpos = xpos_mean + xpos_offset;

    ypos_mean = initial_parameters(3);
    ypos_offset = rand(n_fits(i), 1) - 0.5;
    input_ypos = ypos_mean + ypos_offset;

    parameters.a = initial_parameters(1);
    parameters.x0 = input_xpos;
    parameters.y0 = input_ypos;
    parameters.s = initial_parameters(4);
    parameters.b = initial_parameters(5);
    data = generate_2Dgaussians(parameters, n_fits(i), fit_size);
    data = data + noise_std_dev * randn(fit_size,fit_size,n_fits(i));
    data = permute(data,[2,1,3]);

    initial_parameter_set = repmat(initial_parameters, 1, n_fits(i));

    %% run GpuFit
    [parameters_GpuFit, converged_GpuFit, chisquare_GpuFit, n_iterations_GpuFit, time_GpuFit]...
        = GpuFit(data, sigma, fit_size*fit_size, max_iterations, initial_parameter_set, parameters_to_fit, model_id, estimator_id, tolerance, user_info);
    converged_GpuFit = converged_GpuFit + 1;
    calculated.a = parameters_GpuFit(1:n_parameters:end).';
    calculated.x0 = parameters_GpuFit(2:n_parameters:end).';
    calculated.y0 = parameters_GpuFit(3:n_parameters:end).';
    calculated.s = parameters_GpuFit(4:n_parameters:end).';
    calculated.b = parameters_GpuFit(5:n_parameters:end).';
    
    speed_GpuFit(i) = n_fits(i)/time_GpuFit;
    valid_indices = get_valid_fit_indices(converged_GpuFit, chisquare_GpuFit);
    valid_n_iterations = n_iterations_GpuFit(valid_indices);
    precision_GpuFit(i) = calculate_precision(calculated, parameters, valid_indices);
    mean_iterations_GpuFit(i) = mean(valid_n_iterations);
    print_fit_info(precision_GpuFit(i), time_GpuFit, 'Gpufit', numel(valid_indices)/n_fits(i), mean_iterations_GpuFit(i));
    
    %% run GPU-LMFit
    
    [parameters_GPULMFit, info_GPULMFit, time_GPULMFit] = GPULMFit(...
        data,...
        estimator_id,...
        parameters.s,...
        fit_size);
    calculated.a = parameters_GPULMFit(2:n_parameters:end).';
    calculated.x0 = parameters_GPULMFit(3:n_parameters:end).';
    calculated.y0 = parameters_GPULMFit(4:n_parameters:end).';
    calculated.s = parameters_GPULMFit(5:n_parameters:end).';
    calculated.b = parameters_GPULMFit(1:n_parameters:end).';
    
    speed_GPULMFit(i) = n_fits(i)/time_GPULMFit;
    valid_indices = get_valid_fit_indices(ones(1,n_fits(i)), ones(1,n_fits(i)));
    precision_GPULMFit(i) = calculate_precision(calculated, parameters, valid_indices);
    print_fit_info(precision_GPULMFit(i), time_GPULMFit, 'GPU-LMFit', 0, 0);
    
    
    %% run cminpack
    [parameters_cminpack, info_cminpack, n_iterations_cminpack, time_cminpack]...
        = cminpack(data, initial_parameter_set, model_id, tolerance);
    converged_cminpack = (info_cminpack > 0) & (info_cminpack <= 3);
    calculated.a = parameters_cminpack(1:n_parameters:end).';
    calculated.x0= parameters_cminpack(2:n_parameters:end).';
    calculated.y0 = parameters_cminpack(3:n_parameters:end).';
    calculated.s = parameters_cminpack(4:n_parameters:end).';
    calculated.b = parameters_cminpack(5:n_parameters:end).';
    
    speed_cminpack(i) = n_fits(i)/time_cminpack;
    valid_indices = get_valid_fit_indices(converged_cminpack, ones(1,n_fits(i)));
    valid_n_iterations = n_iterations_cminpack(valid_indices);
    mean_iterations_cminpack(i) = mean(valid_n_iterations);
    precision_cminpack(i) = calculate_precision(calculated, parameters, valid_indices);
    print_fit_info(precision_cminpack(i), time_cminpack, 'cminpack', numel(valid_indices)/n_fits(i), mean_iterations_cminpack(i));
end
%% save test info
info.parameters = parameters;
info.initial_parameters = initial_parameter_set;
info.noise = noise;
info.snr = snr;
info.fit_size = fit_size;
info.n_fits = n_fits;
info.model_id = model_id;

%% output filename
filename = 'GpuFit_GPULMFit_cminpack_speed_nfits';
    
%% save data
save(filename);

%% write file
xlsfilename = [filename '.xls'];

Raw(1:100, 1:100)=deal(NaN);
xlswrite(xlsfilename,Raw,1)

xlscolumns = {'fit count' 'GpuFit' 'GPU-LMFit' 'Minpack'};
xlswrite(xlsfilename,xlscolumns,1,'A1')

xlsmat(:,1) = n_fits;
xlsmat(:,2) = speed_GpuFit;
xlsmat(:,3) = speed_GPULMFit;
xlsmat(:,4) = speed_cminpack;
xlswrite(xlsfilename,xlsmat,1,'A2')

write_test_info(xlsfilename, info);

%% plot
Plot_GpuFit_GPULMFit_cminpack_speed(...
    n_fits,...
    speed_GpuFit,...
    speed_GPULMFit,...
    speed_cminpack)
savefig(filename)
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