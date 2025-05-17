function [x, fval, steps] = newton_raphson(x0, epsilon)
    max_iter = 1000;
    x = x0;
    steps = 0;

    while steps < max_iter
        g = gradfunc(x);
        H = hessianfunc(x);

        if norm(g) <= epsilon
            break;
        end

        % Solve H * p = -g
        p = -H \ g;

        % Line search or fixed step size
        alpha = 1.0;

        x = x + alpha * p;
        steps = steps + 1;

        % Optional: check f difference too
        if abs(func(x + p) - func(x)) <= epsilon
            break;
        end
    end

    fval = func(x);
end
