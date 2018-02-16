function Fitting_Brown_Dennis()

init_factors = [1;10;100];
initial_parameters = single([25; 5; -5; 1]);
ydata = zeros(20,1, 'single');

for i = 1:numel(init_factors)

[parameters_gpufit, states_gpufit(i), chi_square_gpufit, iterations_gpufit(i), time_gpufit]...
    = gpufit(ydata,[], ModelID.BROWN_DENNIS, initial_parameters, 1e-8, 10000,[1; 1; 1; 1], EstimatorID.LSE);

xdata = 0:19;
[parameters_minpack, states_minpack(i), iterations_minpack(i), time_minpack]...
    = cminpack(double(ydata), double(initial_parameters), 11, 1e-30, xdata);

f_estimated_gpufit = calculate_brown_dennis(parameters_gpufit);
resudial_gpufit(i) = norm(f_estimated_gpufit);

f_estimated_minpack = calculate_brown_dennis(parameters_minpack);
resudial_minpack(i) = norm(f_estimated_minpack);

solution = [ -11.5844, 13.1999, -0.406200, 0.240998 ];

f_solution = calculate_brown_dennis(solution);
solution_resudial = norm(f_solution);

f_initial = calculate_brown_dennis(initial_parameters);

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

T_brown_dennis = table(resudials , iterations, states, 'RowNames', methods, 'VariableNames', names);
save('T_brown_dennis', 'T_brown_dennis');

end

function f = calculate_brown_dennis(p)

x = 0:19;

t = x / 5;

arg1 = p(1) + p(2) * t - exp(t);
arg2 = p(3) + p(4) * sin(t) - cos(t);

f = arg1.*arg1 + arg2.*arg2;

end