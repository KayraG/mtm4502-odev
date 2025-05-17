function g = gradfunc(x)
    % Numerical gradient using central difference method
    n = length(x);
    g = zeros(n,1);
    h = 1e-6;  % small perturbation

    for i = 1:n
        dx = zeros(n,1);
        dx(i) = h;
        g(i) = (func(x + dx) - func(x - dx)) / (2*h);
    end
end
