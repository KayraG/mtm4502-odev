function [x, fval, steps] = hestenes_stiefel(x0, epsilon)
    max_iter = 1000;
    x = x0;
    g = gradfunc(x);
    d = -g;
    steps = 0;

    while norm(g) > epsilon && steps < max_iter
        alpha = line_search(x, d); % optional line search
        x_new = x + alpha * d;
        g_new = gradfunc(x_new);
        
        beta = (g_new' * (g_new - g)) / (d' * (g_new - g));
        d = -g_new + beta * d;

        x = x_new;
        g = g_new;
        steps = steps + 1;

        if abs(func(x + d) - func(x)) < epsilon
            break;
        end
    end

    fval = func(x);
end
