function [] = figure10_gpufit_gpulmfit_cminpack_accuracy()

%% number of fits per test point
n_fits = 100000;

%% parameters determining the data to be fit
fit_size = 15;
gauss_amplitude = 500;
gauss_width = 1.0;
gauss_baseline = 10;
noise = 'gauss';

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
gauss_pos_offset_max = 0.0;
initial_guess_offset_frac = 0.3;
snr = 60;

[data, data_parameters, initial_guess_parameters] = ...
    generate_gauss_fit_test_data(n_fits, fit_size, gauss_amplitude, gauss_width, ...
                                 gauss_baseline, noise, gauss_pos_offset_max, ...
                                 initial_guess_offset_frac, snr);


                             
%% run gpufit
[parameters_gpufit, converged_gpufit, chisquare_gpufit, n_iterations_gpufit, time_gpufit]...
    = gpufit(data, weights, model_id, initial_guess_parameters, tolerance, ...
             max_iterations, parameters_to_fit, estimator_id, user_info);

converged_gpufit = converged_gpufit + 1;

chk_gpulmfit = 0;

[valid_gpufit_results, gpufit_abs_precision, mean_n_iterations, valid_indices] = ...
    process_gaussian_fit_results(data_parameters, parameters_gpufit, converged_gpufit, ...
                                 chisquare_gpufit, n_iterations_gpufit, chk_gpulmfit);
                       
print_fit_info(gpufit_abs_precision, time_gpufit, 'Gpufit', numel(valid_indices)/n_fits, mean_n_iterations);

%% run GPU-LMFit
[parameters_gpulmfit, info_gpulmfit, chisquare_gpulmfit, time_gpulmfit] = GPULMFit(...
        data,...
        estimator_id,...
        initial_guess_parameters,...
        fit_size);
    
chisquare_gpulmfit = chisquare_gpulmfit .* chisquare_gpulmfit;
converged_gpulmfit = (info_gpulmfit > 0) & (info_gpulmfit <= 4);
n_iterations_gpulmfit = ones(1,n_fits);
chk_gpulmfit = 0;

[valid_gpulmfit_results, gpulmfit_abs_precision, mean_n_iterations, valid_indices] = ...
    process_gaussian_fit_results(data_parameters, parameters_gpulmfit, converged_gpulmfit, ...
                                 chisquare_gpulmfit, n_iterations_gpulmfit, chk_gpulmfit);

print_fit_info(gpulmfit_abs_precision, time_gpulmfit, 'GPU-LMFit', numel(valid_indices)/n_fits, mean_n_iterations);


%% run cminpack
[parameters_cminpack, info_cminpack, n_iterations_cminpack, time_cminpack]...
    = cminpack(data, initial_guess_parameters, model_id, tolerance);

converged_cminpack = (info_cminpack > 0) & (info_cminpack <= 3);
chisquare_cminpack = ones(1,n_fits);

[valid_cminpack_results, cminpack_abs_precision, mean_n_iterations, valid_indices] = ...
    process_gaussian_fit_results(data_parameters, parameters_cminpack, converged_cminpack, ...
                                 chisquare_cminpack, n_iterations_cminpack);

print_fit_info(cminpack_abs_precision, time_cminpack, 'C Minpack', numel(valid_indices)/n_fits, mean_n_iterations);


%% filename
filename = 'figure10_gpufit_gpulmfit_cminpack_accuracy';

%% write file
xlsfilename = [filename '.xlsx'];

xlscolumns(1,1) = {'Gpufit'};
xlscolumns(1,6) = {'GpuLMFit'};
xlscolumns(1,11) = {'Minpack'};
xlswrite(xlsfilename, xlscolumns, 1, 'A1')

xlscolumns = {'amplitude' 'center x' 'center y' 'width' 'background'};
xlscolumns = [xlscolumns xlscolumns xlscolumns];
xlswrite(xlsfilename,xlscolumns,1,'A2')

xlsmat(:,1) = valid_gpufit_results.a;
xlsmat(:,2) = valid_gpufit_results.x0;
xlsmat(:,3) = valid_gpufit_results.y0;
xlsmat(:,4) = valid_gpufit_results.s;
xlsmat(:,5) = valid_gpufit_results.b;

xlsmat(:,6) = valid_gpulmfit_results.a;
xlsmat(:,7) = valid_gpulmfit_results.x0;
xlsmat(:,8) = valid_gpulmfit_results.y0;
xlsmat(:,9) = valid_gpulmfit_results.s;
xlsmat(:,10) = valid_gpulmfit_results.b;

xlsmat(:,11) = valid_cminpack_results.a;
xlsmat(:,12) = valid_cminpack_results.x0;
xlsmat(:,13) = valid_cminpack_results.y0;
xlsmat(:,14) = valid_cminpack_results.s;
xlsmat(:,15) = valid_cminpack_results.b;

xlswrite(xlsfilename,xlsmat,1,'A3')

            
%% plot
Plot_gpufit_gpulmfit_cmpipack_accuracy(...
    data_parameters,...
    valid_gpufit_results,...
    valid_cminpack_results,...
    valid_gpulmfit_results)
savefig(filename)

end

function [] = Plot_gpufit_gpulmfit_cmpipack_accuracy(...
    parameters,...
    parameters_gpufit,...
    parameters_cminpack,...
    parameters_gpulmfit)
%%
temp_max(2) = max(abs(parameters_gpufit.a));
temp_max(3) = max(abs(parameters_cminpack.a));
temp_max(4) = max(abs(parameters_gpulmfit.a));
amp_max = max(temp_max);

temp_max(2) = max(abs(parameters_gpufit.x0));
temp_max(3) = max(abs(parameters_cminpack.x0));
temp_max(4) = max(abs(parameters_gpulmfit.x0));
x_max = max(temp_max);

temp_max(2) = max(abs(parameters_gpufit.y0));
temp_max(3) = max(abs(parameters_cminpack.y0));
temp_max(4) = max(abs(parameters_gpulmfit.y0));
y_max = max(temp_max);

temp_max(2) = max(abs(parameters_gpufit.s));
temp_max(3) = max(abs(parameters_cminpack.s));
temp_max(4) = max(abs(parameters_gpulmfit.s));
width_max = max(temp_max);

temp_max(2) = max(abs(parameters_gpufit.b));
temp_max(3) = max(abs(parameters_cminpack.b));
temp_max(4) = max(abs(parameters_gpulmfit.b));
bg_max = max(temp_max);

subplot(3,5,0*5+1);
plotSingleHist(parameters_gpufit.a, parameters.a, amp_max);
ylabel({'\fontsize{40}gpufit'; ['{' ' ' '}']}, 'FontWeight', 'bold');
title({'\fontsize{40}amplitude'; ['\fontsize{30}' '{\color{red}' num2str(parameters.a) '}']}, 'FontWeight', 'bold')

subplot(3,5,0*5+2);
plotSingleHist(parameters_gpufit.x0, parameters.x0(1), x_max);
title({'\fontsize{40}center x'; ['\fontsize{30}' '{\color{red}' num2str(parameters.x0(1)) '}']}, 'FontWeight', 'bold')

subplot(3,5,0*5+3);
plotSingleHist(parameters_gpufit.y0, parameters.y0(1), y_max);
title({'\fontsize{40}center y'; ['\fontsize{30}' '{\color{red}' num2str(parameters.y0(1)) '}']}, 'FontWeight', 'bold')

subplot(3,5,0*5+4);
plotSingleHist(parameters_gpufit.s, parameters.s, width_max);
title({'\fontsize{40}width'; ['\fontsize{30}' '{\color{red}' num2str(parameters.s) '}']}, 'FontWeight', 'bold')

subplot(3,5,0*5+5);
plotSingleHist(parameters_gpufit.b, parameters.b, bg_max);
title({'\fontsize{40}background'; ['\fontsize{30}' '{\color{red}' num2str(parameters.b) '}']}, 'FontWeight', 'bold')
%% ----------------
subplot(3,5,1*5+1);
plotSingleHist(parameters_gpulmfit.a, parameters.a, amp_max);
ylabel({'\fontsize{40}GPU-LMFit'; ['{' ' ' '}']}, 'FontWeight', 'bold');

subplot(3,5,1*5+2);
plotSingleHist(parameters_gpulmfit.x0, parameters.x0(1), x_max);

subplot(3,5,1*5+3);
plotSingleHist(parameters_gpulmfit.y0, parameters.y0(1), y_max);

subplot(3,5,1*5+4);
plotSingleHist(parameters_gpulmfit.s, parameters.s, width_max);

subplot(3,5,1*5+5);
plotSingleHist(parameters_gpulmfit.b, parameters.b, bg_max);
%% ----------------
subplot(3,5,2*5+1);
plotSingleHist(parameters_cminpack.a, parameters.a, amp_max);
ylabel({'\fontsize{40}Minpack'; ['{' ' ' '}']}, 'FontWeight', 'bold');
 
subplot(3,5,2*5+2);
plotSingleHist(parameters_cminpack.x0, parameters.x0(1), x_max);

subplot(3,5,2*5+3);
plotSingleHist(parameters_cminpack.y0, parameters.y0(1), y_max);

subplot(3,5,2*5+4);
plotSingleHist(parameters_cminpack.s, parameters.s, width_max);

subplot(3,5,2*5+5);
plotSingleHist(parameters_cminpack.b, parameters.b, bg_max);
end

%%
function[] = plotSingleHist(fitted_parameters, actual_parameters, max_value)

    histogram(fitted_parameters, 50);
    ylim([0 850])
    yL = get(gca,'YLim');
    line([actual_parameters actual_parameters],yL,'Color','r', 'LineWidth', 1);
    m = mean(fitted_parameters);
    s = std(fitted_parameters);
    xlabel({['\fontsize{22}µ: ' num2str(m)] ; ['\fontsize{22}' '{' '  \fontsize{30}\sigma: ' '\fontsize{22}' num2str(s) '}']}, 'FontWeight', 'bold')
    xlim([-max_value+2*actual_parameters max_value])
    current_hist=gca;
    current_hist.FontSize = 16;
    current_hist.LineWidth = 4;
    current_hist.XTick = [];
    current_hist.YTick = [];
end
