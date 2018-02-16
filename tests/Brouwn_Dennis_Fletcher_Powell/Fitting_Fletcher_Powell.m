function Fitting_Fletcher_Powell()

init_factors = [1;10;100];
initial_parameters = single([-1; 0; 0]);
ydata = single([0; 0; 0]);

for i = 1:numel(init_factors)

[parameters_gpufit, states_gpufit(i), chi_square_gpufit, iterations_gpufit(i), time_gpufit(i)]...
    = gpufit(ydata,[], ModelID.FLETCHER_POWELL, initial_parameters, 1e-8, 10000,[1; 1; 1], EstimatorID.LSE);

xdata=[0; 1; 2];
[parameters_minpack, states_minpack(i), iterations_minpack(i), time_minpack(i)]...
    = cminpack(double(ydata), double(initial_parameters), 8, 1e-8, xdata);

f_estimated_gpufit = calculate_fletcher_powell(parameters_gpufit);
resudial_gpufit(i) = norm(f_estimated_gpufit);

f_estimated_minpack = calculate_fletcher_powell(parameters_minpack);
resudial_minpack(i) = norm(f_estimated_minpack);

solution = [1, 0, 0];

f_solution = calculate_fletcher_powell(solution);
solution_resudial = norm(f_solution);

f_initial = calculate_fletcher_powell(initial_parameters);

figure(i);
plot(f_initial);
hold on;
plot(f_solution);
plot(f_estimated_gpufit);
plot(f_estimated_minpack);
legend(...
    'initial',...
    'solution',...
    'estimated gpufit',...
    'estimated minpack');

initial_parameters = initial_parameters * 10;
end

names = {'residuals__1x0_10x0_100x0', 'iterations__1x0_10x0_100x0', 'states__1x0_10x0_100x0'};
methods = {'gpufit'; 'minpack'};
resudials = [resudial_gpufit; resudial_minpack];
iterations = [iterations_gpufit; iterations_minpack];
states = [states_gpufit; states_minpack];

T_fletcher_powell = table(resudials , iterations, states, 'RowNames', methods, 'VariableNames', names);
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
