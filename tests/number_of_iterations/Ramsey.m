function Ramsey()

n_fits = 1;
data = importdata('ramsey_data.mat');

xdata  = data(:,1);
ydata = data(:,2);
amp = (max(ydata) - min(ydata)) / 5;
offset = (max(ydata) + min(ydata)) / 2;
samplerate = length(xdata) / max(xdata);
nsamples = length(xdata);
yfft = abs(fft(ydata));
yfft = yfft(1:floor(nsamples / 2));
frange = samplerate / nsamples * (0:nsamples / 2 - 1);
yfft = yfft(frange > 0.00);
frange = frange(frange > 0.00);
minpeakheigth = (max(yfft) - min(yfft)) / 5;
[peaks, locs] = findpeaks(yfft,'SORTSTR','descend','MINPEAKHEIGHT',minpeakheigth);   
freq = frange(locs(1:2));
t2star1 = 0.75e-06;
phase = 2.8e-08 ./ (2.*pi) .*10.^6;

initial_parameters = [amp amp offset freq(1) freq(2) 1 t2star1 0 0]';

model_id = ModelID.RAMSEY_VAR_P;
estimator_id = EstimatorID.LSE;
tolerance = 1e-8;
max_iterations = 100;

for i = 1:3
    
    %% cpufit
    [parameters_cpufit, states_cpufit(i), chi_square_cpufit, iterations_cpufit(i), time_cpufit(i)]...
        = cpufit(ydata,[], model_id, initial_parameters, tolerance, max_iterations,[], estimator_id, xdata);

    %%gpufit
    [parameters_gpufit, states_gpufit(i), chi_square_gpufit, iterations_gpufit(i), time_gpufit(i)]...
        = gpufit(ydata,[], model_id, initial_parameters, tolerance, max_iterations,[], estimator_id, xdata);

    %% cminpack
    [parameters_minpack, states_minpack(i), iterations_minpack(i), time_minpack(i)]...
        = cminpack(ydata, initial_parameters, model_id, tolerance, xdata);

    %% calculate resudials
    % cpufit
    f_estimated_cpufit = generate_ramsey(parameters_cpufit, xdata);
    resudial_cpufit(i) = norm(f_estimated_cpufit-ydata);
    
    % gpufit
    f_estimated_gpufit = generate_ramsey(parameters_gpufit, xdata);
    resudial_gpufit(i) = norm(f_estimated_gpufit-ydata);

    % minpack
    f_estimated_minpack = generate_ramsey(parameters_minpack, xdata);
    resudial_minpack(i) = norm(f_estimated_minpack-ydata);

    f_initial = generate_ramsey(initial_parameters, xdata);
    
    %% plot
    figure(i);
    plot(f_initial);
    hold on;
    plot(ydata);
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
    initial_parameters(8) = initial_parameters(8) + 1e-4;
    initial_parameters(9) = initial_parameters(9) + 1e-4;
end

names = {'residuals__1x0_10x0_100x0', 'iterations__1x0_10x0_100x0', 'states__1x0_10x0_100x0'};
methods = {'cpufit'; 'gpufit'; 'minpack'};
resudials = [resudial_cpufit; resudial_gpufit; resudial_minpack];
iterations = [iterations_cpufit; iterations_gpufit; iterations_minpack];
states = [states_cpufit; states_gpufit; states_minpack];

T_ramsey = table(resudials , iterations, states, 'RowNames', methods, 'VariableNames', names);
save('T_ramsey','T_ramsey');
    
end

%%
function f = generate_ramsey(parameters, x)

a1 = parameters(1);
a2 = parameters(2);
b = parameters(3);
f1 = parameters(4);
f2 = parameters(5);
p = parameters(6);
t = parameters(7);
phase1 = parameters(8);
phase2 = parameters(9);

f(:,:) =  exp(-(x./t).^p)...
  .* ((a1.*cos(2.*pi.*f1.*(x - phase1)) + a2.*cos(2.*pi.*f2.*(x-phase2))))...
  +  b;

end
