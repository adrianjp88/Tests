%% test parameters
snr = 1;

%% number of fits per test point
n_fits = 1;

%% parameters determining the data to be fit
fit_size = 20;
%% parameters determining how the fit is carried out
weights = [];
max_iterations = 1000;
model_id = ModelID.BROWN_DENNIS; %GAUSS_2D
estimator_id = 0; %LSE
n_parameters = 4;
parameters_to_fit = ones(n_parameters,1);
user_info = [];
tolerance = 1e-8;

data = zeros(fit_size, 1, 'single');
    
initial_guess_parameters = single([25; 5; -5; 1] * 10);
    
    %% run gpufit
    [parameters_gpufit, converged_gpufit, chisquare_gpufit, n_iterations_gpufit, time_gpufit, lambda_info_gpufit]...
        = cpufit(data, weights, model_id, initial_guess_parameters, tolerance, ...
                 max_iterations, parameters_to_fit, estimator_id, user_info);
             
             lambda_info_gpufit = reshape(lambda_info_gpufit, 100, 10);
             lambda_info_gpufit = table(...
                 lambda_info_gpufit(:,1),...
                 lambda_info_gpufit(:,2),...
                 lambda_info_gpufit(:,3),...
                 lambda_info_gpufit(:,4),...
                 lambda_info_gpufit(:,5),...
                 lambda_info_gpufit(:,6),...
                 lambda_info_gpufit(:,7),...
                 lambda_info_gpufit(:,8),...
                 lambda_info_gpufit(:,9),...
                 lambda_info_gpufit(:,10));
             lambda_info_gpufit.Properties.VariableNames...
                 = {'lambda'...
                    'lower_bound'...
                    'upper_bound'...
                    'step_bound'...
                    'predicted_reduction'...
                    'actual_reduction'...
                    'directive_derivative'...
                    'phi'...
                    'chi'...
                    'previous_chi'};
             
    converged_gpufit = converged_gpufit + 1;

    %% run cminpack
    [parameters_cminpack, info_cminpack, n_iterations_cminpack, time_cminpack, lambda_info_cminpack]...
        = cminpack(data, double(initial_guess_parameters), 11, tolerance);
    
    lambda_info_cminpack = reshape(lambda_info_cminpack, 100, 10);
             lambda_info_cminpack = table(...
                 lambda_info_cminpack(:,1),...
                 lambda_info_cminpack(:,2),...
                 lambda_info_cminpack(:,3),...
                 lambda_info_cminpack(:,4),...
                 lambda_info_cminpack(:,5),...
                 lambda_info_cminpack(:,6),...
                 lambda_info_cminpack(:,7),...
                 lambda_info_cminpack(:,8),...
                 lambda_info_cminpack(:,9),...
                 lambda_info_cminpack(:,10));
             lambda_info_cminpack.Properties.VariableNames...
                 = {'lambda'...
                    'lower_bound'...
                    'upper_bound'...
                    'step_bound'...
                    'predicted_reduction'...
                    'actual_reduction'...
                    'directive_derivative'...
                    'phi'...
                    'chi'...
                    'previous_chi'};
    
    converged_cminpack = (info_cminpack > 0) & (info_cminpack <= 3);