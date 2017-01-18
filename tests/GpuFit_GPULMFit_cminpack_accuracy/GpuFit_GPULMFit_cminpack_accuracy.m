function [] = GpuFit_GPULMFit_cminpack_accuracy()
%% set info
n_fits = 100000;
fit_size = 15;

sigma = [];

model_id = 1; %GAUSS_2D
estimator_id = 0; %LSE
n_curve_parameters = 5;
max_iterations = 20;
parameters_to_fit = ones(1,n_curve_parameters);
user_info = 0;
tolerance = 0.0001;
noise = 'gauss';

pos_offset_dist = 1.0;

%% set parameters
x0 = (fit_size-1)/2;
x0 = repmat(x0, n_fits, 1);
y0 = (fit_size-1)/2;
y0 = repmat(y0, n_fits,1);

xpos_offset = 0.0;
ypos_offset = 0.0;

% xpos_offset = pos_offset_dist*rand(n_fits, 1) - (pos_offset_dist/2.0);
% ypos_offset = pos_offset_dist*rand(n_fits, 1) - (pos_offset_dist/2.0);

x0 = x0 + xpos_offset;
y0 = y0 + ypos_offset;

parameters.a = 500;
parameters.x0 = x0;
parameters.y0 = y0;
parameters.s = 1;
parameters.b = 10;

%% set start values

initial_parameters = ones(1,n_curve_parameters);
initial_parameters(1) = parameters.a * 1.1;
initial_parameters(2) = (fit_size-1)/2 * 1.1;
initial_parameters(3) = (fit_size-1)/2 * 1.1;
initial_parameters(4) = parameters.s * 1.1;
initial_parameters(5) = parameters.b * 1.1;
initial_parameters = repmat(initial_parameters, 1, n_fits);

%% calculate noise standard deviations
snr = 100;

gauss_fwtm = 4.292 * parameters.s; % only valid for circular gaussian
fit_area = gauss_fwtm*gauss_fwtm;

mean_amplitude = 2*pi*parameters.a(1)*parameters.s(1)*parameters.s(1)/fit_area;

noise_std_dev = mean_amplitude / snr;

%% generate data
data = generate_2Dgaussians(parameters, n_fits, fit_size);
data = data + noise_std_dev * randn(fit_size,fit_size,n_fits);
data = permute(data,[2,1,3]);

%% run GpuFit
[parameters_GpuFit, converged_GpuFit, chisquare_GpuFit, n_iterations_GpuFit, time_GpuFit]...
    = GpuFit(data, sigma, fit_size*fit_size, max_iterations, initial_parameters, parameters_to_fit, model_id, estimator_id, tolerance, user_info);
converged_GpuFit = converged_GpuFit + 1;
calculated.a = parameters_GpuFit(1:n_curve_parameters:end).';
calculated.x0 = parameters_GpuFit(2:n_curve_parameters:end).';
calculated.y0 = parameters_GpuFit(3:n_curve_parameters:end).';
calculated.s = parameters_GpuFit(4:n_curve_parameters:end).';
calculated.b = parameters_GpuFit(5:n_curve_parameters:end).';

errors_gpufit.a = calculated.a - parameters.a;
errors_gpufit.x0 = calculated.x0 - parameters.x0;
errors_gpufit.y0 = calculated.y0 - parameters.y0;
errors_gpufit.s = calculated.s - parameters.s;
errors_gpufit.b = calculated.b - parameters.b;

valid_indices = get_valid_fit_indices(converged_GpuFit, chisquare_GpuFit);
valid_n_iterations = n_iterations_GpuFit(valid_indices);
precision_GpuFit = calculate_precision(calculated, parameters, valid_indices);
mean_iterations_GpuFit = mean(valid_n_iterations);
print_fit_info(precision_GpuFit, time_GpuFit, 'Gpufit', numel(valid_indices)/n_fits, mean_iterations_GpuFit);

%% run GPU-LMFit
[parameters_GPULMFit, info_GPULMFit, time_GPULMFit] = GPULMFit(...
        data,...
        estimator_id,...
        parameters.s,...
        fit_size);
calculated.a = parameters_GPULMFit(2:n_curve_parameters:end).';
calculated.x0 = parameters_GPULMFit(3:n_curve_parameters:end).';
calculated.y0 = parameters_GPULMFit(4:n_curve_parameters:end).';
calculated.s = parameters_GPULMFit(5:n_curve_parameters:end).';
calculated.b = parameters_GPULMFit(1:n_curve_parameters:end).';

valid_indices = get_valid_fit_indices(ones(1,n_fits), ones(1,n_fits));
precision_GPULMFit = calculate_precision(calculated, parameters, valid_indices);
print_fit_info(precision_GPULMFit, time_GPULMFit, 'GPU-LMFit', 0, 0);


%% run cminpack
[parameters_cminpack, info_cminpack, n_iterations_cminpack, time_cminpack]...
    = cminpack(data, initial_parameters, model_id, tolerance);
converged_cminpack = (info_cminpack > 0) & (info_cminpack <= 3);
calculated.a = parameters_cminpack(1:n_curve_parameters:end).';
calculated.x0= parameters_cminpack(2:n_curve_parameters:end).';
calculated.y0 = parameters_cminpack(3:n_curve_parameters:end).';
calculated.s = parameters_cminpack(4:n_curve_parameters:end).';
calculated.b = parameters_cminpack(5:n_curve_parameters:end).';

valid_indices = get_valid_fit_indices(converged_cminpack, ones(1,n_fits));
valid_n_iterations = n_iterations_cminpack(valid_indices);
mean_iterations_cminpack = mean(valid_n_iterations);
precision_cminpack = calculate_precision(calculated, parameters, valid_indices);
print_fit_info(precision_cminpack, time_cminpack, 'cminpack', numel(valid_indices)/n_fits, mean_iterations_cminpack);

%% save test info
info.parameters = parameters;
info.initial_parameters = initial_parameters;
info.noise = noise;
info.snr = snr;
info.fit_size = fit_size;
info.n_fits = n_fits;
info.model_id = model_id;

%% filename
filename = 'GpuFit_GPULMFit_cminpack_accuracy';

%% save data
save(filename);

%% write file
xlsfilename = [filename '.xlsx'];

xlscolumns(1,1) = {'GpuFit'};
xlscolumns(1,6) = {'GPULMFit'};
xlscolumns(1,11) = {'Minpack'};
xlswrite(xlsfilename, xlscolumns, 1, 'A1')

xlscolumns = {'amplitude' 'center x' 'center y' 'width' 'background'};
xlscolumns = [xlscolumns xlscolumns xlscolumns];
xlswrite(xlsfilename,xlscolumns,1,'A2')

xlsmat(:,1) = parameters_GpuFit(1:5:end);
xlsmat(:,2) = parameters_GpuFit(2:5:end);
xlsmat(:,3) = parameters_GpuFit(3:5:end);
xlsmat(:,4) = parameters_GpuFit(4:5:end);
xlsmat(:,5) = parameters_GpuFit(5:5:end);

xlsmat(:,6) = parameters_GPULMFit(2:5:end);
xlsmat(:,7) = parameters_GPULMFit(2:5:end);
xlsmat(:,8) = parameters_GPULMFit(3:5:end);
xlsmat(:,9) = parameters_GPULMFit(4:5:end);
xlsmat(:,10) = parameters_GPULMFit(1:5:end);

xlsmat(:,11) = parameters_cminpack(1:5:end);
xlsmat(:,12) = parameters_cminpack(2:5:end);
xlsmat(:,13) = parameters_cminpack(3:5:end);
xlsmat(:,14) = parameters_cminpack(4:5:end);
xlsmat(:,15) = parameters_cminpack(5:5:end);

xlswrite(xlsfilename,xlsmat,1,'A3')

%% save info
write_test_info(xlsfilename, info);
            
%% plot
Plot_GpuFit_GPULMFit_cmpipack_accuracy(...
    parameters,...
    parameters_GpuFit,...
    parameters_cminpack,...
    parameters_GPULMFit)
savefig(filename)
end

function [] = Plot_GpuFit_GPULMFit_cmpipack_accuracy(...
    parameters,...
    parameters_GpuFit,...
    parameters_cminpack,...
    parameters_GPULMFit)
%%
temp_max(2) = max(abs(parameters_GpuFit(1:5:end)));
temp_max(3) = max(abs(parameters_cminpack(1:5:end)));
temp_max(4) = max(abs(parameters_GPULMFit(2:5:end)));
amp_max = max(temp_max);

temp_max(2) = max(abs(parameters_GpuFit(2:5:end)));
temp_max(3) = max(abs(parameters_cminpack(2:5:end)));
temp_max(4) = max(abs(parameters_GPULMFit(3:5:end)));
x_max = max(temp_max);

temp_max(2) = max(abs(parameters_GpuFit(3:5:end)));
temp_max(3) = max(abs(parameters_cminpack(3:5:end)));
temp_max(4) = max(abs(parameters_GPULMFit(4:5:end)));
y_max = max(temp_max);

temp_max(2) = max(abs(parameters_GpuFit(4:5:end)));
temp_max(3) = max(abs(parameters_cminpack(4:5:end)));
temp_max(4) = max(abs(parameters_GPULMFit(5:5:end)));
width_max = max(temp_max);

temp_max(2) = max(abs(parameters_GpuFit(5:5:end)));
temp_max(3) = max(abs(parameters_cminpack(5:5:end)));
temp_max(4) = max(abs(parameters_GPULMFit(1:5:end)));
bg_max = max(temp_max);

subplot(3,5,0*5+1);
plotSingleHist(parameters_GpuFit(1:5:end), parameters.a, amp_max);
ylabel({'\fontsize{40}GpuFit'; ['{' ' ' '}']}, 'FontWeight', 'bold');
title({'\fontsize{40}amplitude'; ['\fontsize{30}' '{\color{red}' num2str(parameters.a) '}']}, 'FontWeight', 'bold')

subplot(3,5,0*5+2);
plotSingleHist(parameters_GpuFit(2:5:end), parameters.x0(1), x_max);
title({'\fontsize{40}center x'; ['\fontsize{30}' '{\color{red}' num2str(parameters.x0(1)) '}']}, 'FontWeight', 'bold')

subplot(3,5,0*5+3);
plotSingleHist(parameters_GpuFit(3:5:end), parameters.y0(1), y_max);
title({'\fontsize{40}center y'; ['\fontsize{30}' '{\color{red}' num2str(parameters.y0(1)) '}']}, 'FontWeight', 'bold')

subplot(3,5,0*5+4);
plotSingleHist(parameters_GpuFit(4:5:end), parameters.s, width_max);
title({'\fontsize{40}width'; ['\fontsize{30}' '{\color{red}' num2str(parameters.s) '}']}, 'FontWeight', 'bold')

subplot(3,5,0*5+5);
plotSingleHist(parameters_GpuFit(5:5:end), parameters.b, bg_max);
title({'\fontsize{40}background'; ['\fontsize{30}' '{\color{red}' num2str(parameters.b) '}']}, 'FontWeight', 'bold')
%% ----------------
subplot(3,5,1*5+1);
plotSingleHist(parameters_GPULMFit(2:5:end), parameters.a, amp_max);
ylabel({'\fontsize{40}GPU-LMFit'; ['{' ' ' '}']}, 'FontWeight', 'bold');

subplot(3,5,1*5+2);
plotSingleHist(parameters_GPULMFit(3:5:end), parameters.x0(1), x_max);

subplot(3,5,1*5+3);
plotSingleHist(parameters_GPULMFit(4:5:end), parameters.y0(1), y_max);

subplot(3,5,1*5+4);
plotSingleHist(parameters_GPULMFit(5:5:end), parameters.s, width_max);

subplot(3,5,1*5+5);
plotSingleHist(parameters_GPULMFit(1:5:end), parameters.b, bg_max);
%% ----------------
subplot(3,5,2*5+1);
plotSingleHist(parameters_cminpack(1:5:end), parameters.a, amp_max);
ylabel({'\fontsize{40}Minpack'; ['{' ' ' '}']}, 'FontWeight', 'bold');
 
subplot(3,5,2*5+2);
plotSingleHist(parameters_cminpack(2:5:end), parameters.x0(1), x_max);

subplot(3,5,2*5+3);
plotSingleHist(parameters_cminpack(3:5:end), parameters.y0(1), y_max);

subplot(3,5,2*5+4);
plotSingleHist(parameters_cminpack(4:5:end), parameters.s, width_max);

subplot(3,5,2*5+5);
plotSingleHist(parameters_cminpack(5:5:end), parameters.b, bg_max);
end

%%
function[] = plotSingleHist(fitted_parameters, actually_parameters, max_value)

    histogram(fitted_parameters, 50);
    ylim([0 8500])
    yL = get(gca,'YLim');
    line([actually_parameters actually_parameters],yL,'Color','r', 'LineWidth', 1);
    m = mean(fitted_parameters);
    s = std(fitted_parameters);
    xlabel({['\fontsize{22}µ: ' num2str(m)] ; ['\fontsize{22}' '{' '  \fontsize{30}\sigma: ' '\fontsize{22}' num2str(s) '}']}, 'FontWeight', 'bold')
    xlim([-max_value+2*actually_parameters max_value])
    current_hist=gca;
    current_hist.FontSize = 16;
    current_hist.LineWidth = 4;
    current_hist.XTick = [];
    current_hist.YTick = [];
end
