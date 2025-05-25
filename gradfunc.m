function g = gradfunc(x)
% GRADFUNC calculates the gradient vector for the Kowalik Problem.
%
%   g = gradfunc(x)
%
%   Input:
%       x - A 4x1 column vector representing the current point (x1, x2, x3, x4).
%
%   Output:
%       g - A 4x1 column vector representing the gradient [df/dx1; df/dx2; df/dx3; df/dx4].

% Define the constants ai and bi for the Kowalik Problem
a_vals = [0.1957, 0.1947, 0.1735, 0.16, 0.0844, 0.0627, 0.0456, 0.0342, 0.0323, 0.0235, 0.0246];
b_vals = [0.25, 0.50, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0];

% Initialize gradient vector
g = zeros(4, 1);

% Loop through i from 1 to 11 to sum up contributions to the gradient
for i = 1:11
    ai = a_vals(i);
    bi = b_vals(i);

    % Define D and N for convenience
    D = 1 + x(3) * bi + x(4) * bi^2;
    N = x(1) * (1 + x(2) * bi);

    % Calculate the error term e_i
    e_i = ai - (N / D);

    % Calculate partial derivatives of e_i with respect to x1, x2, x3, x4
    de_i_dx1 = -(1 + x(2) * bi) / D;
    de_i_dx2 = -(x(1) * bi) / D;
    de_i_dx3 = (N * bi) / D^2;
    de_i_dx4 = (N * bi^2) / D^2;

    % Add contribution to the gradient vector (2 * e_i * de_i/dx_j)
    g(1) = g(1) + 2 * e_i * de_i_dx1;
    g(2) = g(2) + 2 * e_i * de_i_dx2;
    g(3) = g(3) + 2 * e_i * de_i_dx3;
    g(4) = g(4) + 2 * e_i * de_i_dx4;
end

end
