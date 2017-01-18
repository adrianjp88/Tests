function [] = MLE_LSE()

n_fits = 1000;
fit_size = 15;

sigma = [];

model_id = 1; %GAUSS_2D
LSE = 0;
MLE = 1;
n_curve_parameters = 5;
max_iterations = 100;
parameters_to_fit = ones(1,n_curve_parameters);
user_info = 0;
tolerance = 0.0001;
noise = 'poisson';

pos_offset_dist = 1.0;

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
parameters.b = 20;

start_values = ones(1,n_curve_parameters);

start_values(2) = (fit_size-1)/2;
start_values(3) = (fit_size-1)/2;
start_values(4) = parameters.s * 1.1;
start_values(5) = parameters.b * 1.0;
start_values = repmat(start_values, 1, n_fits);

n_points = 5;

amp_min = 100;
amp_max = 10000;
log_min = log10(amp_min);
log_max = log10(amp_max);
log_amp = linspace(log_min, log_max, n_points);
amp = 10.^log_amp;

for i = 1:n_points
    parameters.a = amp(i);
    start_values(1:n_curve_parameters:end) = parameters.a * 1.1;
    data = generate_2Dgaussians(parameters, n_fits, fit_size);
    data = poissrnd(data,fit_size,fit_size,n_fits);
    data = permute(data,[2,1,3]);
    
    %% run GpuFit MLE
    [curve_parameters_MLE, converged_MLE, chisquare_MLE, n_iterations_MLE, time_MLE]...
        = GpuFit(data, sigma, fit_size*fit_size, max_iterations, start_values, parameters_to_fit, model_id, MLE, tolerance, user_info);
    converged_MLE = converged_MLE + 1;
    calculated.a = curve_parameters_MLE(1:n_curve_parameters:end).';
    calculated.x0 = curve_parameters_MLE(2:n_curve_parameters:end).';
    calculated.y0 = curve_parameters_MLE(3:n_curve_parameters:end).';
    calculated.s = curve_parameters_MLE(4:n_curve_parameters:end).';
    calculated.b = curve_parameters_MLE(5:n_curve_parameters:end).';
    
    valid_indices = get_valid_fit_indices(converged_MLE, chisquare_MLE);

    valid_n_iterations = n_iterations_MLE(valid_indices);

    mean_iterations_MLE(i) = mean(valid_n_iterations);

    precision_MLE(i) = calculate_precision(calculated, parameters, valid_indices);

    print_fit_info(precision_MLE(i), time_MLE, 'MLE', numel(valid_indices)/n_fits, mean_iterations_MLE(i));
    
    %% run GpuFit LSE
    reshaped_data = reshape(data,1,fit_size*fit_size*n_fits);
    [curve_parameters_LSE, converged_LSE, chisquare_LSE, n_iterations_LSE, time_LSE]...
        = GpuFit(reshaped_data, sigma, fit_size*fit_size, max_iterations, start_values, parameters_to_fit, model_id, LSE, tolerance, user_info);
    converged_LSE = converged_LSE + 1;
    calculated.a = curve_parameters_LSE(1:n_curve_parameters:end).';
    calculated.x0 = curve_parameters_LSE(2:n_curve_parameters:end).';
    calculated.y0 = curve_parameters_LSE(3:n_curve_parameters:end).';
    calculated.s = curve_parameters_LSE(4:n_curve_parameters:end).';
    calculated.b = curve_parameters_LSE(5:n_curve_parameters:end).';
    
    valid_indices = get_valid_fit_indices(converged_LSE, chisquare_LSE);
    
    valid_n_iterations = n_iterations_LSE(valid_indices);

    mean_iterations_LSE(i) = mean(valid_n_iterations);

    precision_LSE(i) = calculate_precision(calculated, parameters, valid_indices);

    print_fit_info(precision_LSE(i), time_LSE, 'LSE', numel(valid_indices)/n_fits, mean_iterations_LSE(i));
    
end
%% save test info
info.parameters = parameters;
info.initial_parameters = start_values;
info.noise = noise;
info.fit_size = fit_size;
info.n_fits = n_fits;
info.model_id = model_id;


%% output filename
filename = 'MLE_LSE';

%% save data
save(filename);

%% write .xls file precision
xlsfilename = [filename '.xls'];

Raw(1:100, 1:100)=deal(NaN);
xlswrite(xlsfilename,Raw,1)

xlscolumns(1,2) = {'MLE'};
xlscolumns(1,7) = {'LSE'};
xlscolumns(2,1) = {'SNR'};
xlswrite(xlsfilename,xlscolumns,1,'A1')

xlscolumns = {'amplitude' 'center x' 'center y' 'width' 'background'};
xlscolumns = [xlscolumns xlscolumns];
xlswrite(xlsfilename,xlscolumns,1,'B2')

xlsmat(:,1) = amp;

xlsmat(:,2) = [precision_MLE.a];
xlsmat(:,3) = [precision_MLE.x0];
xlsmat(:,4) = [precision_MLE.y0];
xlsmat(:,5) = [precision_MLE.s];
xlsmat(:,6) = [precision_MLE.b];

xlsmat(:,7) = [precision_LSE.a];
xlsmat(:,8) = [precision_LSE.x0];
xlsmat(:,9) = [precision_LSE.y0];
xlsmat(:,10) = [precision_LSE.s];
xlsmat(:,11) = [precision_LSE.b];

xlswrite(xlsfilename,xlsmat,1,'A3')
clear xlsmat

write_test_info(xlsfilename, info);

%% write file precision x0
xlsfilename = [filename '_x0.xls'];

xlscolumns = {'SNR','MLE x0','LSE x0'};
xlswrite(xlsfilename,xlscolumns,1,'A1')

xlsmat(:,1) = amp;
xlsmat(:,2) = [precision_MLE.x0];
xlsmat(:,3) = [precision_LSE.x0];

xlswrite(xlsfilename,xlsmat,1,'A2')
clear xlsmat

write_test_info(xlsfilename, info);

%% write file iterations
xlsfilename = [filename '_iterations.xls'];
xlscolumns = {'SNR' 'MLE' 'LSE'};
xlswrite(xlsfilename,xlscolumns,1,'A1')

xlsmat(:,1) = amp;
xlsmat(:,2) = mean_iterations_MLE;
xlsmat(:,3) = mean_iterations_LSE;

xlswrite(xlsfilename,xlsmat,1,'A2')
clear xlsmat

write_test_info(xlsfilename, info);

%% plot
 Plot_MLE_LSE_precision(amp, precision_MLE, precision_LSE)
 savefig(filename)
 Plot_MLE_LSE_precision_x0(amp, [precision_MLE.x0], [precision_LSE.x0])
 savefig([filename '_x0'])
 Plot_MLE_LSE_iterations(amp, mean_iterations_MLE, mean_iterations_LSE)
 savefig([filename '_iterations'])

end

function [] = Plot_MLE_LSE_precision(amp, precision_MLE, precision_LSE)

figure('Name','MLE vs LSE, snr presision','NumberTitle','off');
loglog(...
    amp, [precision_MLE.a], 'red+', ...
    amp, [precision_MLE.x0], 'blue+', ...
    amp, [precision_MLE.s], 'green+', ...
    amp, [precision_MLE.b], 'black+', ...
    amp, [precision_LSE.a], 'redo', ...
    amp, [precision_LSE.x0], 'blueo', ...
    amp, [precision_LSE.s], 'greeno', ...
    amp, [precision_LSE.b], 'blacko', ...
    'LineWidth', 4, ...
    'LineStyle', '-', ...
    'MarkerSize', 20)
xlabel('amplitude')
ylabel('relative standard deviation')
legend(...
    'MLE amplitude',...
    'MLE center',...
    'MLE width',...
    'MLE background',...
    'LSE amplitude',...
    'LSE center',...
    'LSE width',...
    'LSE background')
xlim([-inf inf])
grid on;
box off;
current_figure = gca;
current_figure.FontSize = 25;
current_figure.LineWidth = 4;

end

function [] = Plot_MLE_LSE_precision_x0(amp, precision_MLE_x0, precision_LSE_x0)

figure('Name','MLE vs LSE, snr x center presision','NumberTitle','off');
loglog(...
    amp, precision_MLE_x0, 'red+', ...
    amp, precision_LSE_x0, 'blues', ...
    'LineWidth', 4, ...
    'LineStyle', 'none', ...
    'MarkerSize', 20)
xlabel('amplitude')
ylabel('relative standard deviation')
legend('MLE', 'LSE')
xlim([-inf inf])
grid on;
box off;
current_figure = gca;
current_figure.FontSize = 25;
current_figure.LineWidth = 4;
end

function [] = Plot_MLE_LSE_iterations(amp, iterations_MLE, iterations_LSE)

figure('Name','MLE vs LSE, snr presision','NumberTitle','off');
semilogx(...
    amp, iterations_MLE, 'red+', ...
    amp, iterations_LSE, 'blueo', ...
    'LineWidth', 4, ...
    'LineStyle', '-', ...
    'MarkerSize', 20)
xlabel('amplitude')
ylabel('number of iterations')
legend('MLE', 'LSE')
xlim([-inf inf])
grid on;
box off;
current_figure = gca;
current_figure.FontSize = 25;
current_figure.LineWidth = 4;
end