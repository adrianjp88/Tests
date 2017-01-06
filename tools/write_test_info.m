function [] = write_test_info(filename, info)

architecture = computer('arch');

parameters_names = fieldnames(info.parameters).';
parameters_cell = struct2cell(info.parameters).';

n_parameters = length(parameters_cell);

for i = 1:n_parameters
    parameters_cell{i} = get_range(cell2mat(parameters_cell(i)));
    initial_parameters{i} = get_range(info.initial_parameters(i));
end

if strcmp(info.noise, 'gauss')
    snr = info.snr;
end

fit_size_range = get_range(info.fit_size);
n_fits_range = get_range(info.n_fits);

if strcmp(info.noise, 'gauss')
    snr_range = get_range(snr);
end

switch info.model_id
    case 0
        model = 'Gauss 1D';
    case 1
        model = 'Gauss 2D';
    case 2
        model = 'Gauss 2D elliptic';
    case 3
        model = 'Gauss 2D rotated';
    case 4
        model = 'Cauchy 2D elliptic';
end

warning('off','MATLAB:xlswrite:AddSheet');

xlsrows = {...
    'fit size'... 
    'number of fits'...
    'architecture'... 
    'model'... 
    'parameters'... 
    'test parameter values '... 
    'initial parameter values '... 
    'noise'}.';
if strcmp(info.noise, 'gauss')
    xlsrows = [xlsrows ; 'SNR'];
end

xlsmat(:,1) = xlsrows;

xlsmat(1,2) = {fit_size_range};
xlsmat(2,2) = {n_fits_range};
xlsmat(3,2) = {architecture};
xlsmat(4,2) = {model};
xlsmat(5,2:n_parameters+1) = parameters_names;
xlsmat(6,2:n_parameters+1) = parameters_cell;
xlsmat(7,2:n_parameters+1) = initial_parameters;
xlsmat(8,2) = {info.noise};
if strcmp(info.noise, 'gauss')
    xlsmat(9,2) = {snr_range};
end

Raw(1:100, 1:100)=deal(NaN);
xlswrite(filename,Raw,2)
xlswrite(filename, xlsmat, 2, 'A1')

end

function [range] = get_range(array)

if length(array) < 1
    range = 'not defined';
elseif length(array) == 1
    range = num2str(array);
elseif abs(min(array) - max(array)) < 0.0000000001
    range = num2str(array(1));
elseif length(array) > 1
    range = ['[' num2str(min(array)) '...' num2str(max(array)) ']'];
end

end