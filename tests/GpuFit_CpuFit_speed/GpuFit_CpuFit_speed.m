function [] = GpuFit_CpuFit_speed()
fit_size = 5;
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

%% test case
LogNFitsMin = 0;
LogNFitsMax = 6;
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
gauss_fwtm = 4.292 * initial_parameters(4); % only valid for circular gaussian
fit_area = gauss_fwtm*gauss_fwtm;
mean_amplitude...
    =2*pi*initial_parameters(1)*initial_parameters(4)*initial_parameters(4)/fit_area;
noise_std_dev = mean_amplitude ./ snr;

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

    %% run CpuFit
    [parameters_CpuFit, converged_CpuFit, chisquare_CpuFit, n_iterations_CpuFit, time_CpuFit]...
        = CpuFit(data, sigma, max_iterations, initial_parameter_set, parameters_to_fit, model_id, estimator_id, tolerance, user_info);
    converged_CpuFit = converged_CpuFit + 1;
    calculated.a = parameters_CpuFit(1:n_parameters:end).';
    calculated.x0 = parameters_CpuFit(2:n_parameters:end).';
    calculated.y0 = parameters_CpuFit(3:n_parameters:end).';
    calculated.s = parameters_CpuFit(4:n_parameters:end).';
    calculated.b = parameters_CpuFit(5:n_parameters:end).';
    
    speed_CpuFit(i) = n_fits(i)/time_CpuFit;
    valid_indices = get_valid_fit_indices(converged_CpuFit, chisquare_CpuFit);
    valid_n_iterations = n_iterations_CpuFit(valid_indices);
    precision_CpuFit(i) = calculate_precision(calculated, parameters, valid_indices);
    mean_iterations_CpuFit(i) = mean(valid_n_iterations);
    print_fit_info(precision_CpuFit(i), time_CpuFit, 'Gpufit', numel(valid_indices)/n_fits(i), mean_iterations_CpuFit(i));
      
    %% run GpuFit
    [parameters_GpuFit, converged_GpuFit, chisquare_GpuFit, n_iterations_GpuFit, time_GpuFit]...
        = GpuFit(data, sigma, max_iterations, initial_parameter_set, parameters_to_fit, model_id, estimator_id, tolerance, user_info);
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
    
end
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

%% write file
xlsfilename = [filename '.xls'];

Raw(1:100, 1:100)=deal(NaN);
xlswrite(xlsfilename,Raw,1)

xlscolumns = {'number of fits' 'GpuFit' 'CpuFit'};
xlswrite(xlsfilename,xlscolumns,1,'A1')
xlsmat(:,1) = n_fits;
xlsmat(:,2) = speed_GpuFit;
xlsmat(:,3) = speed_CpuFit;
xlswrite(xlsfilename,xlsmat,1,'A2')

write_test_info(xlsfilename, info);

%% plot
Plot_GpuFit_CpuFit_speed(n_fits, speed_GpuFit, speed_CpuFit);
savefig(filename)

end

function [] = Plot_GpuFit_CpuFit_speed(n_fits, speed_GpuFit, speed_CpuFit)

while length(speed_GpuFit) > length(speed_CpuFit)
    speed_CpuFit = [speed_CpuFit speed_CpuFit(end)];
end

figure('Name', 'GpuFit vs CpuFit', 'NumberTitle', 'off');
semilogx(...
    n_fits, speed_GpuFit, 'red.-', ...
    n_fits, speed_CpuFit, 'blue.-', ...
    'LineWidth', 8)
xlabel('number of fits')
ylabel('fits per second')
legend('GpuFit', 'CpuFit')

grid on;
box off;
current_figure = gca;
current_figure.FontSize = 37;
current_figure.LineWidth = 4;

end