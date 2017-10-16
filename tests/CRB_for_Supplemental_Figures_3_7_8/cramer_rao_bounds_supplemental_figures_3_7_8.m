function cramer_rao_bounds_supplemental_figures_3_7_8()
% Cramer-Rao Bounds for fit precision (standard deviation of x center
% coordinate) for some supplemental figures

%% supplemental figure 3
f3.likelihood = 'gaussian';
f3.A = 500;
f3.L = 15;
f3.c = (f3.L-1)/2 * [1,1];
f3.s = 1;
f3.b = 10;


% taken from figure5_GpuFit_cminpack.m
n_graph_points = 30;
snr_min = 5.0;
snr_max = 100000000;
log_min = log10(snr_min);
log_max = log10(snr_max);
log_snr = linspace(log_min, log_max, n_graph_points);
f3.snr = 10.^log_snr;

f3.noise_std = f3.A ./ (f3.snr * log(10.0));

f3_crb = compute_crb(f3);

%% supplemental figure 7
f7.likelihood = 'gaussian';
f7.A = 500;
f7.L = 15;
f7.c = (f7.L-1)/2 * [1,1];
f7.s = 1;
f7.b = 10;

% taken from figure9_GpuFit_GPULMFIT_cminpack_precision.m
n_graph_points = 20;
snr_min = 2;
snr_max = 1000;
log_min = log10(snr_min);
log_max = log10(snr_max);
log_snr = linspace(log_min, log_max, n_graph_points);
f7.snr = 10.^log_snr;

f7.noise_std = f7.A ./ (f7.snr * log(10.0));

f7_crb = compute_crb(f7);

%% supplemental figure 8
f8.likelihood = 'poisson';

n_points = 20;
amp_min = 5;
amp_max = 200;
log_min = log10(amp_min);
log_max = log10(amp_max);
log_amp = linspace(log_min, log_max, n_points);
f8.A = 10.^log_amp;

f8.L = 15;
f8.c = (f8.L-1)/2 * [1,1];
f8.s = 2;
f8.b = 5;

f8_crb = compute_crb(f8);

% write to xls
filename = 'crb_supplemental_figures_3_7_8.xls';

xlscolumns = {'SNR', 'CRB-SupFig3', 'SNR', 'CRB-SupFig7', 'Gauss-Peak-Amplitude', 'CRB-SupFig8'};
xlswrite(filename, xlscolumns, 1, 'A1');

xlswrite(filename,[f3.snr', f3_crb'],1,'A2');
xlswrite(filename,[f7.snr', f7_crb'],1,'C2');
xlswrite(filename,[f8.A', f8_crb'],1,'E2');
end

function crb = compute_crb(f)

% grid
[x,y] = ndgrid(0:f.L-1,0:f.L-1);

% compute fisher information

switch f.likelihood
    case 'gaussian'
        % compute model and model derivatives
        exp_factor = exp(-((x-f.c(1)).^2+(y-f.c(2)).^2)/(2*f.s^2));
        % mi = f.A * exp_factor + f.b;
        dmidxc = f.A * exp_factor .* (x - f.c(1)) /f.s^2;        
        I = 1 ./ (f.noise_std.^2) * sum(dmidxc(:).^2);
    case 'poisson'
        n = length(f.A);
        I = zeros(1, n);
        for i = 1 : n
            Ai = f.A(i);
            exp_factor = exp(-((x-f.c(1)).^2+(y-f.c(2)).^2)/(2*f.s^2));
            mi = Ai * exp_factor + f.b;
            dmidxc = Ai * exp_factor .* (x - f.c(1)) /f.s^2;                    
            I(i) = sum(dmidxc(:).^2 ./ mi(:));
        end
    otherwise
        error('unknown likelihood');
end

crb = 1 ./ sqrt(I);

end