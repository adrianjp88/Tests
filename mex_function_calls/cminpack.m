function [out,info, iteration_count, time]...
    = cminpack(data, initial_parameters, functionID, tolerance)

versionID = 6;

%% conficure parameters
if functionID == 1
    n_curve_parameters = 5;
elseif functionID == 2
    n_curve_parameters = 6;
end

data_size = size(data);

fit_size = data_size(1) * data_size(2);
if (ndims(data) == 3)
    n_fits = data_size(3);
else 
    n_fits = 1;
end

if length(initial_parameters) < n_fits * n_curve_parameters
    initial_parameters = repmat(initial_parameters, 1, n_fits);
end

if ~isa(data, 'double')
    data = double(data);
end

                
%% run cminpack
tic;
[out, info, iteration_count]...
    = cminpackMex(versionID, data, fit_size, n_fits, n_curve_parameters, initial_parameters, functionID, tolerance);
time = toc;

end
