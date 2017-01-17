function [precision] = calculate_abs_precision(fit_results, data_parameters)

    precision.a = std(calc_a-data_parameters.a) / data_parameters.a;
    precision.x0 = std(calc_x0-data_parameters.x0(fit_valid_index));
    precision.y0 = std(calc_y0-data_parameters.y0(fit_valid_index));
    precision.s = std(calc_s-data_parameters.s) / data_parameters.s;
    precision.b = std(calc_b-data_parameters.b) / data_parameters.b;

end
