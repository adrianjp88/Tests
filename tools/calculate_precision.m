function [precision] = calculate_precision(calculated, parameters, fit_valid_index)

    calc_a  = calculated.a(fit_valid_index);
    calc_x0 = calculated.x0(fit_valid_index);
    calc_y0 = calculated.y0(fit_valid_index);
    calc_s  = calculated.s(fit_valid_index);
    calc_b  = calculated.b(fit_valid_index);

    precision.a = std(calc_a-parameters.a) / parameters.a;
    precision.x0 = std(calc_x0-parameters.x0(fit_valid_index));
    precision.y0 = std(calc_y0-parameters.y0(fit_valid_index));
    precision.s = std(calc_s-parameters.s) / parameters.s;
    precision.b = std(calc_b-parameters.b) / parameters.b;

end
