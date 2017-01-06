function [] = GpuFit_cminpack()

n_fits = 100;
fit_size = 15;

sigma = ones(1,fit_size*fit_size);

model_id = 1; %GAUSS_2D
estimator_id = 0; %LSE
n_parameters = 5;
max_iterations = 20;
parameters_to_fit = ones(1,n_parameters);
user_info = 0;
tolerance = 0.0001;
noise = 'gauss';

pos_offset_dist = 2;

xpos_mean = (fit_size-1)/2;
xpos_offset = pos_offset_dist*rand(n_fits, 1) - (pos_offset_dist/2.0);
input_xpos = xpos_mean + xpos_offset;

ypos_mean = (fit_size-1)/2;
ypos_offset = pos_offset_dist*rand(n_fits, 1) - (pos_offset_dist/2.0);
input_ypos = ypos_mean + ypos_offset;

parameters.a = 500;
parameters.x0 = input_xpos;
parameters.y0 = input_ypos;
parameters.s = 1;
parameters.b = 10;

initial_parameters = ones(1,n_parameters);
initial_parameters(1) = parameters.a * 1.1;
initial_parameters(2) = (fit_size-1)/2;
initial_parameters(3) = (fit_size-1)/2;
initial_parameters(4) = parameters.s * 1.1;
initial_parameters(5) = parameters.b * 1.1;
initial_parameters = repmat(initial_parameters, 1, n_fits);

gauss_fwtm = 4.292 * parameters.s % only valid for circular gaussian
fit_area = gauss_fwtm*gauss_fwtm;

mean_amplitude = 2*pi*parameters.a(1)*parameters.s(1)*parameters.s(1)/fit_area;

n_points = 20;

snr_min = 2;
snr_max = 100000000;
log_min = log10(snr_min);
log_max = log10(snr_max);
log_snr = linspace(log_min, log_max, n_points);
snr = 10.^log_snr;

noise_std_dev = mean_amplitude ./ snr;

for i = 1:n_points
    
    data = generate_2Dgaussians(parameters, n_fits, fit_size);
    data = data + noise_std_dev(i) * randn(fit_size,fit_size,n_fits);
    data = permute(data,[2,1,3]);
    
    %% run GpuFit
    [parameters_GpuFit, converged_GpuFit, chisquare_GpuFit, n_iterations_GpuFit, time_GpuFit]...
        = GpuFit(data, sigma, max_iterations, initial_parameters, parameters_to_fit, model_id, estimator_id, tolerance, user_info);
    converged_GpuFit = converged_GpuFit + 1;
    calculated.a = parameters_GpuFit(1:n_parameters:end).';
    calculated.x0 = parameters_GpuFit(2:n_parameters:end).';
    calculated.y0 = parameters_GpuFit(3:n_parameters:end).';
    calculated.s = parameters_GpuFit(4:n_parameters:end).';
    calculated.b = parameters_GpuFit(5:n_parameters:end).';
    
    valid_indices = get_valid_fit_indices(converged_GpuFit, chisquare_GpuFit);
    valid_n_iterations = n_iterations_GpuFit(valid_indices);
    precision_GpuFit(i) = calculate_precision(calculated, parameters, valid_indices);
    mean_iterations_GpuFit(i) = mean(valid_n_iterations);
    print_fit_info(precision_GpuFit(i), time_GpuFit, 'Gpufit', numel(valid_indices)/n_fits, mean_iterations_GpuFit(i));

    %% run cminpack
    [parameters_cminpack, info_cminpack, n_iterations_cminpack, time_cminpack]...
        = cminpack(data, initial_parameters, model_id, tolerance);
    converged_cminpack = (info_cminpack > 0) & (info_cminpack <= 3);
    calculated.a = parameters_cminpack(1:n_parameters:end).';
    calculated.x0= parameters_cminpack(2:n_parameters:end).';
    calculated.y0 = parameters_cminpack(3:n_parameters:end).';
    calculated.s = parameters_cminpack(4:n_parameters:end).';
    calculated.b = parameters_cminpack(5:n_parameters:end).';
    
    valid_indices = get_valid_fit_indices(converged_cminpack, ones(1,n_fits));
    valid_n_iterations = n_iterations_cminpack(valid_indices);
    mean_iterations_cminpack(i) = mean(valid_n_iterations);
    precision_cminpack(i) = calculate_precision(calculated, parameters, valid_indices);
    print_fit_info(precision_cminpack(i), time_cminpack, 'cminpack', numel(valid_indices)/n_fits, mean_iterations_cminpack(i));
end
%% save test info
info.parameters = parameters;
info.initial_parameters = initial_parameters;
info.noise = noise;
info.snr = snr;
info.fit_size = fit_size;
info.n_fits = n_fits;
info.model_id = model_id;

%% output filename
filename = 'GpuFit_cminpack';

%% save data
save(filename);

%% write file precision
xlsfilename = [filename '.xls'];

Raw(1:100, 1:100)=deal(NaN);
xlswrite(xlsfilename,Raw,1)

xlscolumns = {'SNR' 'GpuFit' 'Minpack'};
xlswrite(xlsfilename,xlscolumns,1,'A1')

xlsmat(:,1) = snr;
xlsmat(:,2) = [precision_GpuFit.x0];
xlsmat(:,3) = [precision_cminpack.x0];
xlswrite(xlsfilename,xlsmat,1,'A2')
clear xlsmat

write_test_info(xlsfilename, info);

%% write file iterations
xlsfilename = [filename '_iterations.xls'];
xlscolumns = {'SNR' 'GpuFit' 'Minpack'};
xlswrite(xlsfilename,xlscolumns,1,'A1')

xlsmat(:,1) = snr;
xlsmat(:,2) = mean_iterations_GpuFit;
xlsmat(:,3) = mean_iterations_cminpack;

xlswrite(xlsfilename,xlsmat,1,'A2')
clear xlsmat

write_test_info(xlsfilename, info);

%% plot
Plot_GpuFit_cmpipack_precision(snr, [precision_GpuFit.x0], [precision_cminpack.x0])
savefig(filename)
Plot_GpuFit_cmpipack_iterations(snr, mean_iterations_cminpack, mean_iterations_GpuFit)
savefig([filename '_iterations'])

end

function [] = Plot_GpuFit_cmpipack_iterations(...
    snr,...
    iterations_cminpack,...
    iterations_GpuFit)

figure('Name','GpuFit vs cminpack, snr presision','NumberTitle','off');
semilogx(...
    snr, iterations_GpuFit, 'red+', ...
    snr, iterations_cminpack, 'blueo', ...
    'LineWidth', 4, ...
    'LineStyle', '-', ...
    'MarkerSize', 20)
xlabel('SNR')
ylabel('number of iterations')
legend('GpuFit', 'Minpack')
xlim([-inf inf])
grid on;
box off;
current_figure = gca;
current_figure.FontSize = 25;
current_figure.LineWidth = 4;
end

function [] = Plot_GpuFit_cmpipack_precision(...
    snr,...
    precision_GpuFit,...
    precision_cminpack)

figure('Name','GpuFit vs cminpack, snr presision','NumberTitle','off');
loglog(...
    snr, precision_GpuFit, 'red+', ...
    snr, precision_cminpack, 'blues', ...
    'LineWidth', 4, ...
    'LineStyle', '-', ...
    'MarkerSize', 20)
xlabel('SNR')
ylabel('relative standard deviation')
legend(...
    'GpuFit center x',...
    'Minpack center x')
xlim([-inf inf])
grid on;
box off;
current_figure = gca;
current_figure.FontSize = 25;
current_figure.LineWidth = 4;

end