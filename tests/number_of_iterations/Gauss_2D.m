function Gauss_2D()

n_fits = 1;
n_points = 15;

xdata = 0:n_points-1;

gauss_amplitude = 500;
gauss_width = 1.0;
gauss_baseline = 10;
noise = 'gauss';
gauss_pos_offset_max = 0.5;
initial_guess_offset_frac = 0.5;
snr = 10;
rng(1);

[ydata, true_parameters, initial_parameters] = ...
generate_gauss_fit_test_data(n_fits, n_points, gauss_amplitude, gauss_width, ...
                             gauss_baseline, noise, gauss_pos_offset_max, ...
                             initial_guess_offset_frac, snr);
ydata = double(ydata);
initial_parameters = double(initial_parameters);

model_id = ModelID.GAUSS_2D;
estimator_id = EstimatorID.LSE;

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
        = cminpack(ydata, initial_parameters, model_id, 1e-8);

    %% calculate resudials
    % cpufit
    f_estimated_cpufit = reshape(...
        generate_gauss2d(parameters_cpufit, n_points),...
        n_points * n_points, 1);
    resudial_cpufit(i) = norm(f_estimated_cpufit-ydata);
    
    % gpufit
    f_estimated_gpufit = reshape(...
        generate_gauss2d(parameters_gpufit, n_points),...
        n_points * n_points, 1);
    resudial_gpufit(i) = norm(f_estimated_gpufit-ydata);

    % minpack
    f_estimated_minpack = reshape(...
        generate_gauss2d(parameters_minpack, n_points),...
        n_points * n_points, 1);
    resudial_minpack(i) = norm(f_estimated_minpack-ydata);

    solution(1) = true_parameters.a;
    solution(2) = true_parameters.x0;
    solution(3) = true_parameters.y0;
    solution(4) = true_parameters.s;
    solution(5) = true_parameters.b;

    f_solution = reshape(...
        generate_gauss2d(solution, n_points),...
        n_points * n_points, 1);
    resudial_solution = norm(f_solution-ydata);

    f_initial = generate_gauss2d(initial_parameters, n_points);

    %% scale initial parameters
    initial_parameters(1) = initial_parameters(1) * 10;
    initial_parameters(2) = initial_parameters(2) * 1.2;
    initial_parameters(3) = initial_parameters(3) * 1.2;
    initial_parameters(4) = initial_parameters(4) * 5;
    initial_parameters(5) = initial_parameters(5) * 10;
end

names = {'residuals__1x0_10x0_100x0', 'iterations__1x0_10x0_100x0', 'states__1x0_10x0_100x0'};
methods = {'cpufit'; 'gpufit'; 'minpack'};
resudials = [resudial_cpufit; resudial_gpufit; resudial_minpack];
iterations = [iterations_cpufit; iterations_gpufit; iterations_minpack];
states = [states_cpufit; states_gpufit; states_minpack];

T_gauss_2d = table(resudials , iterations, states, 'RowNames', methods, 'VariableNames', names);
save('T_gauss_2d','T_gauss_2d');
    
end

%%
function f = generate_gauss2d(parameters, size)
    
[x,y] = meshgrid(0:size-1,0:size-1);

x = repmat(x, [1 1 1])';
y = repmat(y, [1 1 1])';

a = parameters(1);
x0 = parameters(2);
y0 = parameters(3);
s = parameters(4);
b = parameters(5);

f(:,:) = a*exp(-1/2*((x-x0)/s).^2)...
       .*exp(-1/2*((y-y0)/s).^2)...
       + b;
end
