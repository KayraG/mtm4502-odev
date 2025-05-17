function H = hessianfunc(x)
    % Numerical Hessian using central difference method
    n = length(x);
    H = zeros(n,n);
    h = 1e-5;

    for i = 1:n
        for j = 1:n
            dx_i = zeros(n,1);
            dx_j = zeros(n,1);
            dx_i(i) = h;
            dx_j(j) = h;

            H(i,j) = (func(x + dx_i + dx_j) - func(x + dx_i - dx_j) - ...
                      func(x - dx_i + dx_j) + func(x - dx_i - dx_j)) / (4*h^2);
        end
    end
end
