function Nonlinear()
    clc; clear;

    rng(0); % for reproducibility

    algorithms = {@newton_raphson, @hestenes_stiefel, @polak_ribiere, @fleetcher_reeves};
    algo_names = ["Newton-Raphson", "Hestenes-Stiefel", "Polak-Ribiere", "Fletcher-Reeves"];
    
    % 3 initial guesses from uniform distribution in [-5,5]
    initial_points = -5 + 10 * rand(4, 3); 

    epsilon = 1e-4;

    for a = 1:length(algorithms)
        fprintf("\n=== %s ===\n", algo_names(a));
        for i = 1:3
            x0 = initial_points(:, i);
            fprintf("\nInitial Point %d:\n", i);
            tic;
            [xmin, fmin, steps] = algorithms{a}(x0, epsilon);
            t = toc;
            fprintf("Minimizer x = [%s]\n", num2str(xmin', '%.4f '));
            fprintf("f(x) = %.4f\n", fmin);
            fprintf("Steps: %d, Time: %.4f sec\n", steps, t);
        end
    end
end
