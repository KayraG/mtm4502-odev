% NONLINEAR.M - Main script for optimizing the Kowalik Problem using
%               various algorithms.

clear; clc; close all;

% --- Problem Parameters ---
epsilon = 1e-4; % Absolute error bound for stopping criteria
max_iterations = 5000; % Maximum iterations to prevent infinite loops

% Domain bounds for Kowalik Problem (0 <= xi <= 0.42 for i=1,2,3,4)
x_min_bound = 0;
x_max_bound = 0.42;
n = 4; % Dimension of the problem

% --- Initial Guesses ---
% Generate 3 initial guesses from a uniform distribution within the bounds
% x0 = U_n(x0,min, x0,max)
num_initial_guesses = 3;
initial_guesses = zeros(n, num_initial_guesses);
rng('default'); % For reproducibility

for k = 1:num_initial_guesses
    % Generate random initial points uniformly within the specified bounds
    initial_guesses(:, k) = x_min_bound + (x_max_bound - x_min_bound) * rand(n, 1);
end

% --- Optimization Algorithms ---

% List of algorithms to run
algorithms = {'Newton-Raphson', 'Hestenes-Stiefel', 'Polak-Ribiere', 'Fletcher-Reeves'};

% Store results for plotting and benchmarking
results = struct();

for alg_idx = 1:length(algorithms)
    current_algorithm = algorithms{alg_idx};
    fprintf('--- Running %s ---\n', current_algorithm);

    results.(current_algorithm).iterations = [];
    results.(current_algorithm).exec_time = [];
    results.(current_algorithm).final_fval = [];
    results.(current_algorithm).final_x = [];
    results.(current_algorithm).fval_history = cell(num_initial_guesses, 1);
    results.(current_algorithm).grad_norm_history = cell(num_initial_guesses, 1);

    for guess_idx = 1:num_initial_guesses
        x0 = initial_guesses(:, guess_idx);
        fprintf('  Initial Guess %d: [%s]\n', guess_idx, sprintf('%.4f ', x0));

        x_k = x0;
        f_k = func(x_k);
        grad_k = gradfunc(x_k);
        norm_grad_k = norm(grad_k);

        fval_history = [f_k];
        grad_norm_history = [norm_grad_k];
        path_history = x_k'; % Store the path for visualization if needed

        tic; % Start timer

        for iter = 1:max_iterations
            % Store previous values for stopping criteria
            x_prev = x_k;
            f_prev = f_k;

            % --- Determine Search Direction and Step Size based on Algorithm ---
            if strcmp(current_algorithm, 'Newton-Raphson')
                H_k = hessianfunc(x_k);

                % Add a small multiple of identity to Hessian for robustness if it's singular or ill-conditioned
                % This helps prevent issues with non-positive definite Hessians far from the minimum
                if cond(H_k) > 1e10 || any(eig(H_k) <= 0)
                    H_k = H_k + eye(n) * 1e-6;
                end

                p_k = -H_k \ grad_k; % Newton direction
                alpha_k = 1; % Full Newton step (often line search is used in practice)
                             % For this problem, we'll start with full step and adjust if needed

                % Simple backtracking line search for robustness
                % This ensures that the function value decreases
                c1 = 0.0001; % Armijo condition constant
                rho = 0.5;   % Step size reduction factor
                while func(x_k + alpha_k * p_k) > f_k + c1 * alpha_k * grad_k' * p_k
                    alpha_k = rho * alpha_k;
                    if alpha_k < 1e-10 % Prevent too small step size
                        break;
                    end
                end

            elseif strcmp(current_algorithm, 'Hestenes-Stiefel')
                if iter == 1
                    p_k = -grad_k; % Steepest descent for first iteration
                else
                    y_prev = grad_k - grad_prev;
                    beta_k = (grad_k' * y_prev) / (y_prev' * p_prev);
                    p_k = -grad_k + beta_k * p_prev;
                end

                % Line search (e.g., Armijo backtracking)
                alpha_k = 1;
                c1 = 0.0001;
                rho = 0.5;
                while func(x_k + alpha_k * p_k) > f_k + c1 * alpha_k * grad_k' * p_k
                    alpha_k = rho * alpha_k;
                    if alpha_k < 1e-10
                        break;
                    end
                end

            elseif strcmp(current_algorithm, 'Polak-Ribiere')
                if iter == 1
                    p_k = -grad_k; % Steepest descent for first iteration
                else
                    y_prev = grad_k - grad_prev;
                    beta_k = (grad_k' * y_prev) / (norm(grad_prev)^2);
                    beta_k = max(0, beta_k); % Ensure non-negative beta (restart if needed)
                    p_k = -grad_k + beta_k * p_prev;
                end

                % Line search
                alpha_k = 1;
                c1 = 0.0001;
                rho = 0.5;
                while func(x_k + alpha_k * p_k) > f_k + c1 * alpha_k * grad_k' * p_k
                    alpha_k = rho * alpha_k;
                    if alpha_k < 1e-10
                        break;
                    end
                end

            elseif strcmp(current_algorithm, 'Fletcher-Reeves')
                if iter == 1
                    p_k = -grad_k; % Steepest descent for first iteration
                else
                    beta_k = (norm(grad_k)^2) / (norm(grad_prev)^2);
                    p_k = -grad_k + beta_k * p_prev;
                end

                % Line search
                alpha_k = 1;
                c1 = 0.0001;
                rho = 0.5;
                while func(x_k + alpha_k * p_k) > f_k + c1 * alpha_k * grad_k' * p_k
                    alpha_k = rho * alpha_k;
                    if alpha_k < 1e-10
                        break;
                    end
                end
            end

            % Store previous gradient and search direction for next iteration
            grad_prev = grad_k;
            p_prev = p_k;

            % Update x_k
            x_k = x_k + alpha_k * p_k;

            % Ensure x_k stays within bounds (projection)
            x_k = max(x_min_bound, min(x_max_bound, x_k));

            % Recalculate function value and gradient at new point
            f_k = func(x_k);
            grad_k = gradfunc(x_k);
            norm_grad_k = norm(grad_k);

            % Store history
            fval_history = [fval_history; f_k];
            grad_norm_history = [grad_norm_history; norm_grad_k];
            path_history = [path_history; x_k'];

            % --- Check Stopping Criteria ---
            % ||gradient|| <= epsilon AND |f(xk+1) - f(xk)| <= epsilon
            if norm_grad_k <= epsilon && abs(f_k - f_prev) <= epsilon
                fprintf('    Converged in %d iterations.\n', iter);
                break;
            end

            if iter == max_iterations
                fprintf('    Max iterations reached (%d) without full convergence.\n', max_iterations);
            end
        end

        exec_time = toc; % End timer

        % Store results for this initial guess
        results.(current_algorithm).iterations = [results.(current_algorithm).iterations; iter];
        results.(current_algorithm).exec_time = [results.(current_algorithm).exec_time; exec_time];
        results.(current_algorithm).final_fval = [results.(current_algorithm).final_fval; f_k];
        results.(current_algorithm).final_x = [results.(current_algorithm).final_x; x_k'];
        results.(current_algorithm).fval_history{guess_idx} = fval_history;
        results.(current_algorithm).grad_norm_history{guess_idx} = grad_norm_history;

        fprintf('    Final f(x): %.4e, Final x: [%s], Iterations: %d, Time: %.4f s\n', ...
            f_k, sprintf('%.4f ', x_k), iter, exec_time);
    end
    fprintf('\n');
end

% --- Plotting Results ---
colors = {'b', 'r', 'g', 'm'}; % Colors for algorithms
line_styles = {'-', '--', ':', '-.'}; % Line styles for initial guesses

% Figure 1: Function Value Convergence
figure;
hold on;
title('Function Value Convergence for Kowalik Problem');
xlabel('Iteration');
ylabel('Function Value');
set(gca, 'YScale', 'log'); % Log scale for better visualization of convergence
grid on;

legend_entries_fval = {};
for alg_idx = 1:length(algorithms)
    current_algorithm = algorithms{alg_idx};
    for guess_idx = 1:num_initial_guesses
        plot(1:length(results.(current_algorithm).fval_history{guess_idx}), ...
             results.(current_algorithm).fval_history{guess_idx}, ...
             'Color', colors{alg_idx}, 'LineStyle', line_styles{guess_idx}, 'LineWidth', 1.5);
        legend_entries_fval{end+1} = sprintf('%s - Guess %d', current_algorithm, guess_idx);
    end
end
legend(legend_entries_fval, 'Location', 'best');
hold off;

% Figure 2: Gradient Norm Convergence
figure;
hold on;
title('Gradient Norm Convergence for Kowalik Problem');
xlabel('Iteration');
ylabel('||∇f(x)||');
set(gca, 'YScale', 'log'); % Log scale for better visualization of convergence
grid on;

legend_entries_grad = {};
for alg_idx = 1:length(algorithms)
    current_algorithm = algorithms{alg_idx};
    for guess_idx = 1:num_initial_guesses
        plot(1:length(results.(current_algorithm).grad_norm_history{guess_idx}), ...
             results.(current_algorithm).grad_norm_history{guess_idx}, ...
             'Color', colors{alg_idx}, 'LineStyle', line_styles{guess_idx}, 'LineWidth', 1.5);
        legend_entries_grad{end+1} = sprintf('%s - Guess %d', current_algorithm, guess_idx);
    end
end
legend(legend_entries_grad, 'Location', 'best');
hold off;

% --- Benchmarking Table (for Report) ---
fprintf('\n--- Benchmarking Results ---\n');
fprintf('%-20s | %-15s | %-15s | %-20s | %-25s\n', 'Algorithm', 'Initial Guess', 'Iterations', 'Execution Time (s)', 'Final f(x) (4 sig figs)');
fprintf(repmat('-', 100, 1));
fprintf('\n');

for alg_idx = 1:length(algorithms)
    current_algorithm = algorithms{alg_idx};
    for guess_idx = 1:num_initial_guesses
        fprintf('%-20s | %-15d | %-15d | %-15.4f | %-25.4e\n', ...
            current_algorithm, ...
            guess_idx, ...
            results.(current_algorithm).iterations(guess_idx), ...
            results.(current_algorithm).exec_time(guess_idx), ...
            results.(current_algorithm).final_fval(guess_idx));
    end
end

