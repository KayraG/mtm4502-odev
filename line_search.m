function alpha = line_search(x, d)
    alpha = 1.0;
    rho = 0.5;  % reduction factor
    c = 1e-4;

    while func(x + alpha * d) > func(x) + c * alpha * (gradfunc(x)' * d)
        alpha = rho * alpha;
        if alpha < 1e-6
            break;
        end
    end
end
