data = importdata('data.mat');

 scaling = 62500000;
 scaling = 1;
        xdata  = data(:,1)*scaling;
        ydata = data(:,2);
 
    
%% Multiplies data by #fits, needed for performance testing
    number_fits = 1;
    ydata_mat = repmat(ydata, [1, number_fits]);

%% calculate fft (need to feed frequencies as initial parameters)
    samplerate = length(xdata) / max(xdata);                      % sample rate (Hz)
    nsamples = length(xdata);
    yfft = abs(fft(ydata));
    yfft = yfft(1:floor(nsamples / 2));                           % discard half of points due to nyquist
    frange = samplerate / nsamples * (0:nsamples / 2 - 1);        % generate frequency range in MHz

    yfft = yfft(frange > 0.00);                                   % discard frequencies below 1 MHz
    frange = frange(frange > 0.00);

%% estimate initial values    
    
    % get doublet hyperfine frequencies from fft
    minpeakheigth = (max(yfft) - min(yfft)) / 5;
    [peaks, locs] = findpeaks(yfft,'SORTSTR','descend','MINPEAKHEIGHT',minpeakheigth);   
    
    freq = frange(locs(1:2));
    freqamp = peaks(1:2);   
     
    % Estimates for amplitudes, offset and phases; t2* currently hardcoded
    amp = (max(ydata) - min(ydata)) / 5;
    offset = (max(ydata) + min(ydata)) / 2; 

    phase = 2.8e-08 ./ (2.*pi) .*10.^6;  % Chinese Paper relation, crudest approximation: phi = detuning * pulse length
    
    t2star1 =0.75e-06;
%% fitting

    % parameters: [A1 A2 c f1 f2 t2star x1 x2]
    %startpoints = [amp amp offset freq(1) freq(2) 1 t2star1 phase phase];
    startpoints = [amp amp offset freq(1) freq(2) 1 t2star1*scaling 0 0];
    
% startpoints = [...
%     0.0217773902098990,...
%     -0.0137695590306338,...
%     -0.00472719405093488,...
%     4987771.34558005,...
%     8015555.12534668,...
%     1.09058846045393,...
%     7.50364681222176e-07,...
%     0.00441853099480146,...
%     0.00439933503950098];

    parameters_to_fit = ones(numel(startpoints),1);

    startpoints_mat = repmat(double(startpoints'), [1, number_fits]);
    
    single_data = single(ydata);
    double_data = double(single_data);
    
    single_startpoints = single(startpoints_mat);
    double_startpoints = double(single_startpoints);
    
    single_tol = single(1e-19);
    double_tol = double(single_tol);
    
    single_xdata = single(xdata);
    double_xdata = double(single_xdata);
    
    [parameters_cpufit, states_cpufit, chi, n_iterations_cpufit, time_cpufit]...
        = cpufit(single_data,[], ModelID.RAMSEY_REDUCED, single_startpoints, single_tol, 1000, parameters_to_fit, EstimatorID.LSE, single_xdata);
    
%     lambda_info_cpufit = reshape(lambda_info_cpufit, 1000, 10);
%              lambda_info_cpufit = table(...
%                  lambda_info_cpufit(:,1),...
%                  lambda_info_cpufit(:,2),...
%                  lambda_info_cpufit(:,3),...
%                  lambda_info_cpufit(:,4),...
%                  lambda_info_cpufit(:,5),...
%                  lambda_info_cpufit(:,6),...
%                  lambda_info_cpufit(:,7),...
%                  lambda_info_cpufit(:,8),...
%                  lambda_info_cpufit(:,9),...
%                  lambda_info_cpufit(:,10));
%              lambda_info_cpufit.Properties.VariableNames...
%                  = {'lambda'...
%                     'lower_bound'...
%                     'upper_bound'...
%                     'step_bound'...
%                     'predicted_reduction'...
%                     'actual_reduction'...
%                     'directive_derivative'...
%                     'phi'...
%                     'chi'...
%                     'previous_chi'};
             
    converged_cpufit = states_cpufit + 1;
    
    [parameters_cminpack, states_cminpack, n_iterations_cminpack, time_cminpack, lambda_info_cminpack]...
                = cminpack(ydata, startpoints, 12, 1e-30, xdata);
            
            lambda_info_cminpack = reshape(lambda_info_cminpack, 1000, 10);
             lambda_info_cminpack = table(...
                 lambda_info_cminpack(:,1),...
                 lambda_info_cminpack(:,2),...
                 lambda_info_cminpack(:,3),...
                 lambda_info_cminpack(:,4),...
                 lambda_info_cminpack(:,5),...
                 lambda_info_cminpack(:,6),...
                 lambda_info_cminpack(:,7),...
                 lambda_info_cminpack(:,8),...
                 lambda_info_cminpack(:,9),...
                 lambda_info_cminpack(:,10));
             lambda_info_cminpack.Properties.VariableNames...
                 = {'lambda'...
                    'lower_bound'...
                    'upper_bound'...
                    'step_bound'...
                    'predicted_reduction'...
                    'actual_reduction'...
                    'directive_derivative'...
                    'phi'...
                    'chi'...
                    'previous_chi'};