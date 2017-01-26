function [data, true_parameters, initial_guess_parameters] = ...
    generate_gauss_fit_test_data(n_fits, fit_size, mean_amplitude, mean_width, ...
                                 mean_baseline, noise_type, center_pos_offset_max, ...
                                 initial_guess_max_offset_frac, snr)

    n_parameters = 5;

    gauss_xpos_mean = (fit_size-1)/2;
    gauss_ypos_mean = (fit_size-1)/2;

    gauss_pos_x = gauss_xpos_mean + ( 2.0*center_pos_offset_max*(rand(n_fits, 1) - 0.5) );
    gauss_pos_y = gauss_ypos_mean + ( 2.0*center_pos_offset_max*(rand(n_fits, 1) - 0.5) );

    initial_guess_xpos_offset_max = initial_guess_max_offset_frac * mean_width;
    initial_guess_ypos_offset_max = initial_guess_max_offset_frac * mean_width;
    initial_guess_ampl_offset_max = initial_guess_max_offset_frac * mean_amplitude;
    initial_guess_width_offset_max = initial_guess_max_offset_frac * mean_width;
    initial_guess_baseline_offset_max = initial_guess_max_offset_frac * mean_baseline;

    initial_guess_xpos_offset = initial_guess_xpos_offset_max ...
                                * 2.0 * (rand(n_fits, 1) - 0.5);

    initial_guess_ypos_offset = initial_guess_ypos_offset_max ...
                                * 2.0 * (rand(n_fits, 1) - 0.5);

    initial_guess_ampl_offset = initial_guess_ampl_offset_max ...
                                * 2.0 * (rand(n_fits, 1) - 0.5);

    initial_guess_width_offset = initial_guess_width_offset_max ...
                                 * 2.0 * (rand(n_fits, 1) - 0.5);

    initial_guess_baseline_offset = initial_guess_baseline_offset_max ...
                                    * 2.0 * (rand(n_fits, 1) - 0.5);

    initial_guess_xpos = gauss_pos_x + initial_guess_xpos_offset;
    initial_guess_ypos = gauss_pos_y + initial_guess_ypos_offset;
    initial_guess_ampl = mean_amplitude + initial_guess_ampl_offset; 
    initial_guess_width = mean_width + initial_guess_width_offset; 
    initial_guess_baseline = mean_baseline + initial_guess_baseline_offset;

    initial_guess_parameters = ones(n_parameters,n_fits,'single');

    initial_guess_parameters(1,:) = initial_guess_ampl;
    initial_guess_parameters(2,:) = initial_guess_xpos;
    initial_guess_parameters(3,:) = initial_guess_ypos;
    initial_guess_parameters(4,:) = initial_guess_width;
    initial_guess_parameters(5,:) = initial_guess_baseline;

    true_parameters.a  = mean_amplitude;
    true_parameters.x0 = gauss_pos_x;
    true_parameters.y0 = gauss_pos_y;
    true_parameters.s  = mean_width;
    true_parameters.b  = mean_baseline;

    noise_std_dev = mean_amplitude ./ (snr * log(10.0));

    %% generate the data

    data = generate_2Dgaussians(true_parameters, n_fits, fit_size, noise_type, noise_std_dev);

end
