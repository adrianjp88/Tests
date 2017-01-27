function [out,info, iteration_count, time]...
    = cminpack(data, initial_parameters, functionID, tolerance)

versionID = 6;

%% configure parameters
if functionID == 1
    n_parameters = 5;
elseif functionID == 2
    n_parameters = 6;
end

data_size = size(data);

fit_size = data_size(1);
if (ndims(data) == 2)
    n_fits = data_size(2);
else 
    n_fits = 1;
end

if ~isa(data, 'double')
    data = cast(data,'double');
end

%% run cminpack
tic;
[out, info, iteration_count]...
    = cminpackMex(versionID, tmp_data, fit_size, n_fits, n_parameters, tmp_init_params, functionID, tolerance);
time = toc;

out = reshape(out, n_parameters, n_fits);

end
