function Gauss_1D()

    n_points = 5;
    
	true_parameters = [4; 2; 0.5; 1];
    xdata = 0:n_points-1;
	ydata = generate_gauss1d(true_parameters, xdata);

    model_id = ModelID.GAUSS_1D;
    estimator_id = EstimatorID.LSE;

	initial_parameters = [ 2; 1.5; 0.3; 1 ];
	max_iterations = 1000;

for i = 1:3
    
    %% cpufit
    [parameters_cpufit, states_cpufit(i), chi_square_cpufit, iterations_cpufit(i), time_cpufit(i)]...
        = cpufit(ydata,[], model_id, initial_parameters, 1e-8, max_iterations,[], estimator_id);

    %%gpufit
    [parameters_gpufit, states_gpufit(i), chi_square_gpufit, iterations_gpufit(i), time_gpufit(i)]...
        = gpufit(ydata,[], model_id, initial_parameters, 1e-8, max_iterations,[], estimator_id);

    %% cminpack
    [parameters_minpack, states_minpack(i), iterations_minpack(i), time_minpack(i)]...
        = cminpack(ydata, initial_parameters, model_id, 1e-8, xdata);

    %% calculate resudials
    % cpufit
    f_estimated_cpufit = generate_gauss1d(parameters_cpufit, xdata);
    resudial_cpufit(i) = norm(f_estimated_cpufit-ydata);
    
    % gpufit
    f_estimated_gpufit = generate_gauss1d(parameters_gpufit, xdata);
    resudial_gpufit(i) = norm(f_estimated_gpufit-ydata);

    % minpack
    f_estimated_minpack = generate_gauss1d(parameters_minpack, xdata);
    resudial_minpack(i) = norm(f_estimated_minpack-ydata);

    solution = true_parameters;

    f_solution = generate_gauss1d(solution, xdata);
    solution_resudial = norm(f_solution);

    f_initial = generate_gauss1d(initial_parameters, xdata);

    %% plot
    figure(i);
    plot(f_initial);
    hold on;
    plot(f_solution);
    plot(f_estimated_cpufit);
    plot(f_estimated_gpufit);
    plot(f_estimated_minpack);
    legend(...
        'initial',...
        'solution',...
        'estimated cpufit',...
        'estimated gpufit',...
        'estimated minpack');

    %% scale initial parameters
    initial_parameters(1) = initial_parameters(1) * 10;
end

names = {'residuals__1x0_10x0_100x0', 'iterations__1x0_10x0_100x0', 'states__1x0_10x0_100x0'};
methods = {'cpufit'; 'gpufit'; 'minpack'};
resudials = [resudial_cpufit; resudial_gpufit; resudial_minpack];
iterations = [iterations_cpufit; iterations_gpufit; iterations_minpack];
states = [states_cpufit; states_gpufit; states_minpack];

T_gauss_1d = table(resudials , iterations, states, 'RowNames', methods, 'VariableNames', names);
save('T_gauss_1d','T_gauss_1d');
    
end

%%
function f = generate_gauss1d(parameters, x)
    
    a = parameters(1);
    x0 = parameters(2);
    s = parameters(3);
    b = parameters(4);
    
    f(:,:) = a*exp(-1/2*((x-x0)/s).^2) + b;
    f = f';
end
