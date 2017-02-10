function [] = figure11_MLE_LSE()

%% test data

n_points = 20;
amp_min = 100;
amp_max = 10000;
log_min = log10(amp_min);
log_max = log10(amp_max);
log_amp = linspace(log_min, log_max, n_points);
amp = 10.^log_amp;

%% number of fits per test point
n_fits = 10000;

%% parameters determining the data to be fit
fit_size = 15;
gauss_width = 1.0;
gauss_baseline = 10;
noise = 'poisson';

%% parameters determining how the fit is carried out
weights = [];
max_iterations = 20;
model_id = 1; %GAUSS_2D
n_parameters = 5;
parameters_to_fit = ones(1,n_parameters,'int32')';
user_info = [];
tolerance = 0.0001;

%% parameters determining the randomness of the data
gauss_pos_offset_max = 1.0;
initial_guess_offset_frac = 0.3;
snr = 1.0;

for i = 1:numel(amp)

    tmp_ampl = amp(i);
    
    [data, data_parameters, initial_guess_parameters] = ...
        generate_gauss_fit_test_data(n_fits, fit_size, tmp_ampl, gauss_width, ...
                                     gauss_baseline, noise, gauss_pos_offset_max, ...
                                     initial_guess_offset_frac, snr);
    

    %% run gpufit MLE
    
    tmp_estimator = 1; 
             
    [parameters_MLE, converged_MLE, chisquare_MLE, n_iterations_MLE, time_MLE]...
        = gpufit(data, weights, model_id, initial_guess_parameters, tolerance, ...
                 max_iterations, parameters_to_fit, tmp_estimator, user_info);
               
             
    converged_MLE = converged_MLE + 1;
             
    chk_gpulmfit = 0;

    [valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(data_parameters, parameters_MLE, converged_MLE, ...
                                     chisquare_MLE, n_iterations_MLE, chk_gpulmfit);
                       
    speed_MLE(i) = n_fits/time_MLE;
    precision_MLE(i) = gpufit_abs_precision;
    mean_iterations_MLE(i) = mean_n_iterations;

    print_fit_info(gpufit_abs_precision, time_MLE, 'MLE', numel(valid_indices)/n_fits, mean_n_iterations);

    %% run gpufit LSE

    tmp_estimator = 0; 
             
    [parameters_LSE, converged_LSE, chisquare_LSE, n_iterations_LSE, time_LSE]...
        = gpufit(data, weights, model_id, initial_guess_parameters, tolerance, ...
                 max_iterations, parameters_to_fit, tmp_estimator, user_info);
    
    converged_LSE = converged_LSE + 1;
             
    chk_gpulmfit = 0;

    [valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
        process_gaussian_fit_results(data_parameters, parameters_LSE, converged_LSE, ...
                                     chisquare_LSE, n_iterations_LSE, chk_gpulmfit);
                       
    speed_LSE(i) = n_fits/time_LSE;
    precision_LSE(i) = gpufit_abs_precision;
    mean_iterations_LSE(i) = mean_n_iterations;

    print_fit_info(gpufit_abs_precision, time_LSE, 'LSE', numel(valid_indices)/n_fits, mean_n_iterations);

end

%% output filename
filename = 'figure11_MLE_LSE';

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

%% write file precision x0
xlsfilename = [filename '_x0.xls'];

xlscolumns = {'SNR','MLE x0','LSE x0'};
xlswrite(xlsfilename,xlscolumns,1,'A1')

xlsmat(:,1) = amp;
xlsmat(:,2) = [precision_MLE.x0];
xlsmat(:,3) = [precision_LSE.x0];

xlswrite(xlsfilename,xlsmat,1,'A2')
clear xlsmat

%% write file iterations
xlsfilename = [filename '_iterations.xls'];
xlscolumns = {'SNR' 'MLE' 'LSE'};
xlswrite(xlsfilename,xlscolumns,1,'A1')

xlsmat(:,1) = amp;
xlsmat(:,2) = mean_iterations_MLE;
xlsmat(:,3) = mean_iterations_LSE;

xlswrite(xlsfilename,xlsmat,1,'A2')
clear xlsmat

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