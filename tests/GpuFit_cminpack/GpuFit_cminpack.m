function [] = GpuFit_cminpack()

n_graph_points = 20;
snr_min = 2;
snr_max = 100000000;
log_min = log10(snr_min);
log_max = log10(snr_max);
log_snr = linspace(log_min, log_max, n_graph_points);
snr = 10.^log_snr;

n_fits = 1000;
fit_size = 15;

sigma = [];
max_iterations = 20;
model_id = 1; %GAUSS_2D
estimator_id = 0; %LSE
n_parameters = 5;
parameters_to_fit = ones(1,n_parameters);
user_info = 0;
tolerance = 0.0001;
noise = 'gauss';

gauss_xpos_mean = (fit_size-1)/2;
gauss_ypos_mean = (fit_size-1)/2;

gauss_pos_offset_max = 1.0;

gauss_pos_x = gauss_xpos_mean + ( 2.0*gauss_pos_offset_max*(rand(n_fits, 1) - 0.5) );
gauss_pos_y = gauss_ypos_mean + ( 2.0*gauss_pos_offset_max*(rand(n_fits, 1) - 0.5) );

gauss_amplitude = 500;
gauss_width = 1.0;
gauss_baseline = 10;

initial_guess_offset_frac = 0.25;

initial_guess_xpos_offset_max = initial_guess_offset_frac * gauss_xpos_mean;
initial_guess_ypos_offset_max = initial_guess_offset_frac * gauss_ypos_mean;
initial_guess_ampl_offset_max = initial_guess_offset_frac * gauss_amplitude;
initial_guess_width_offset_max = initial_guess_offset_frac * gauss_width;
initial_guess_baseline_offset_max = initial_guess_offset_frac * gauss_baseline;

initial_guess_xpos_offset = initial_guess_xpos_offset_max ...
                            * 2.0 * (rand(n_fits, 1) - 0.5);
                  
initial_guess_ypos_offset = initial_guess_ypos_offset_max ...
                            * 2.0 * (rand(n_fits, 1) - 0.5);

initial_guess_ampl_offset = initial_guess_ampl_offset_max ...
                            * 2.0 * (rand(n_fits, 1) - 0.5);
                        
initial_guess_width_offset = initial_guess_width_offset_max ...
                             * 2.0 * (rand(n_fits, 1) - 0.5);
                        
initial_guess_baseline_offset = initial_guess_baseline_offset_max ...
                                * 2.0 * (rand(n_fits, 1) - 0.5);
                        
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

initial_guess_parameters = ones(1,n_parameters*n_fits);

for i = 1:n_fits
    tmp_index = (i-1) * n_fits
    initial_guess_parameters(1 + tmp_index) = initial_guess_ampl(i);
    initial_guess_parameters(2 + tmp_index) = initial_guess_xpos(i);
    initial_guess_parameters(3 + tmp_index) = initial_guess_ypos(i);
    initial_guess_parameters(4 + tmp_index) = initial_guess_width(i);
    initial_guess_parameters(5 + tmp_index) = initial_guess_baseline(i);
end

gauss_fwtm = 4.292 * gauss_width % only valid for circular gaussian
fit_area = 3.1415 * (gauss_fwtm/2.0) * (gauss_fwtm/2.0);

mean_amplitude = 2*pi*gauss_amplitude*gauss_width*gauss_width/fit_area;

noise_std_dev = mean_amplitude / snr;

for i = 1:n_graph_points
    
    data = generate_2Dgaussians(data_parameters, n_fits, fit_size);
    data = data + noise_std_dev * randn(fit_size,fit_size,n_fits);
    data = permute(data,[2,1,3]);
    
    %% run GpuFit
    [parameters_GpuFit, converged_GpuFit, chisquare_GpuFit, n_iterations_GpuFit, time_GpuFit]...
        = GpuFit(data, sigma, fit_size*fit_size, max_iterations, initial_guess_parameters, parameters_to_fit, model_id, estimator_id, tolerance, user_info);

    converged_GpuFit = converged_GpuFit + 1;
    
    valid_indices = get_valid_fit_results(converged_GpuFit, chisquare_GpuFit);
    valid_n_iterations = n_iterations_GpuFit(valid_indices);
    
    gpufit_results.a = parameters_GpuFit(1:n_parameters:end).';
    gpufit_results.x0 = parameters_GpuFit(2:n_parameters:end).';
    gpufit_results.y0 = parameters_GpuFit(3:n_parameters:end).';
    gpufit_results.s = parameters_GpuFit(4:n_parameters:end).';
    gpufit_results.b = parameters_GpuFit(5:n_parameters:end).';
    
    valid_gpufit_results = gpufit_results(valid_indices);
    
    precision_GpuFit(i) = calculate_precision(valid_gpufit_results, data_parameters);
    
    mean_iterations_GpuFit(i) = mean(valid_n_iterations);

    print_fit_info(precision_GpuFit(i), time_GpuFit, 'Gpufit', numel(valid_indices)/n_fits, mean_iterations_GpuFit(i));

    %% run cminpack
    [parameters_cminpack, info_cminpack, n_iterations_cminpack, time_cminpack]...
        = cminpack(data, initial_guess_parameters, model_id, tolerance);
    converged_cminpack = (info_cminpack > 0) & (info_cminpack <= 3);
    calculated.a = parameters_cminpack(1:n_parameters:end).';
    calculated.x0= parameters_cminpack(2:n_parameters:end).';
    calculated.y0 = parameters_cminpack(3:n_parameters:end).';
    calculated.s = parameters_cminpack(4:n_parameters:end).';
    calculated.b = parameters_cminpack(5:n_parameters:end).';
    
    valid_indices = get_valid_fit_results(converged_cminpack, ones(1,n_fits));
    valid_n_iterations = n_iterations_cminpack(valid_indices);
    mean_iterations_cminpack(i) = mean(valid_n_iterations);
    precision_cminpack(i) = calculate_precision(calculated, data_parameters, valid_indices);
    print_fit_info(precision_cminpack(i), time_cminpack, 'cminpack', numel(valid_indices)/n_fits, mean_iterations_cminpack(i));

end

%% save test info
info.parameters = data_parameters;
info.initial_guess_parameters = initial_guess_parameters;
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