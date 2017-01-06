function [] = print_fit_info(precision, time, version, frac_converged, mean_iterations)

    fprintf('***%s***\t', version);
    fprintf('x precision: %g\t|\t', precision.x0);
    fprintf('time: %g\t|\t', time);
    fprintf('converged: %f \t|\t', frac_converged);
    fprintf('iterations: %f\n', mean_iterations);

end