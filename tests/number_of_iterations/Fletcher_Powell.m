function Fletcher_Powell()

n_fits_gpu = 1000000;
n_fits_cpu = 1;

initial_parameters_gpu = repmat([-1; 0; 0], 1, n_fits_gpu);
initial_parameters_cpu = repmat([-1; 0; 0], 1, n_fits_cpu);

ydata_gpu = repmat([0; 0; 0], 1, n_fits_gpu);
ydata_cpu = repmat([0; 0; 0], 1, n_fits_cpu);

for i = 1:1
    
    %% cpufit
    [parameters_cpufit, states_cpufit(:,i), chi_square_cpufit, iterations_cpufit(:,i), time_cpufit(i)]...
        = cpufit(ydata_cpu,[], ModelID.FLETCHER_POWELL, initial_parameters_cpu, 1e-8, 10000,[1; 1; 1], EstimatorID.LSE);

    %%gpufit
    [parameters_gpufit, states_gpufit(:,i), chi_square_gpufit, iterations_gpufit(:,i), time_gpufit(i)]...
        = gpufit(ydata_gpu,[], ModelID.FLETCHER_POWELL, initial_parameters_gpu, 1e-8, 10000,[1; 1; 1], EstimatorID.LSE);

    %% cminpack
    xdata=[0; 1; 2];
    [parameters_minpack, states_minpack(:,i), iterations_minpack(:,i), time_minpack(i)]...
        = cminpack(ydata_cpu, initial_parameters_cpu, ModelID.FLETCHER_POWELL, 1e-8, xdata);

    %% calculate resudials
    % cpufit
    f_estimated_cpufit = calculate_fletcher_powell(parameters_cpufit(:,1));
    resudial_cpufit(i) = norm(f_estimated_cpufit);
    
    % gpufit
    f_estimated_gpufit = calculate_fletcher_powell(parameters_gpufit(:,1));
    resudial_gpufit(i) = norm(f_estimated_gpufit);

    % minpack
    f_estimated_minpack = calculate_fletcher_powell(parameters_minpack(:,1));
    resudial_minpack(i) = norm(f_estimated_minpack);

    solution = [1, 0, 0];

    f_solution = calculate_fletcher_powell(solution);
    solution_resudial = norm(f_solution);

    f_initial = calculate_fletcher_powell(initial_parameters_gpu);
    
    %% calculate speed
    speed_cpufit = n_fits_cpu ./ time_cpufit;
    speed_gpufit = n_fits_gpu ./ time_gpufit;
    speed_minpack = n_fits_cpu ./ time_minpack;

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
    initial_parameters_gpu = initial_parameters_gpu * 10;
    initial_parameters_cpu = initial_parameters_cpu * 10;
end

names = {'residuals__1x0_10x0_100x0', 'iterations__1x0_10x0_100x0', 'states__1x0_10x0_100x0', 'speed__1x0_10x0_100x0'};
methods = {'cpufit'; 'gpufit'; 'minpack'};
resudials = [resudial_cpufit; resudial_gpufit; resudial_minpack];
iterations = [iterations_cpufit(1,:); iterations_gpufit(1,:); iterations_minpack(1,:)];
states = [states_cpufit(1,:); states_gpufit(1,:); states_minpack(1,:)];
speed = [speed_cpufit; speed_gpufit; speed_minpack];

T_fletcher_powell = table(resudials , iterations, states, speed, 'RowNames', methods, 'VariableNames', names);
save('T_fletcher_powell','T_fletcher_powell');
    
end

%%
function f = calculate_fletcher_powell(p)

f(1) = 10*(p(3) - 10*theta(p(1),p(2)));
f(2) = 10*(sqrt(p(1)^2+p(2)^2) - 1);
f(3) = p(3);

end

function t = theta(p1,p2)

t = 1/(2*pi) * atan(p2/p1);

if p1 < 0
    t = t + 0.5;
end

end
