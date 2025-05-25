function H = hessianfunc(x)
% HESSIANFUNC calculates the Hessian matrix for the Kowalik Problem.
%
%   H = hessianfunc(x)
%
%   Input:
%       x - A 4x1 column vector representing the current point (x1, x2, x3, x4).
%
%   Output:
%       H - A 4x4 symmetric matrix representing the Hessian.

% Define the constants ai and bi for the Kowalik Problem
a_vals = [0.1957, 0.1947, 0.1735, 0.16, 0.0844, 0.0627, 0.0456, 0.0342, 0.0323, 0.0235, 0.0246];
b_vals = [0.25, 0.50, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0];

% Initialize Hessian matrix
H = zeros(4, 4);

% Loop through i from 1 to 11 to sum up contributions to the Hessian
for i = 1:11
    ai = a_vals(i);
    bi = b_vals(i);

    % Define D and N for convenience
    D = 1 + x(3) * bi + x(4) * bi^2;
    N = x(1) * (1 + x(2) * bi);

    % Calculate the error term e_i
    e_i = ai - (N / D);

    % First partial derivatives of e_i
    de_i_dx1 = -(1 + x(2) * bi) / D;
    de_i_dx2 = -(x(1) * bi) / D;
    de_i_dx3 = (N * bi) / D^2;
    de_i_dx4 = (N * bi^2) / D^2;

    % Second partial derivatives of e_i
    d2e_i_dx1_dx1 = 0;
    d2e_i_dx2_dx2 = 0;
    d2e_i_dx3_dx3 = -2 * N * bi^2 / D^3;
    d2e_i_dx4_dx4 = -2 * N * bi^4 / D^3;

    d2e_i_dx1_dx2 = -bi / D;
    d2e_i_dx1_dx3 = bi * (1 + x(2) * bi) / D^2;
    d2e_i_dx1_dx4 = bi^2 * (1 + x(2) * bi) / D^2;

    d2e_i_dx2_dx3 = x(1) * bi^2 / D^2;
    d2e_i_dx2_dx4 = x(1) * bi^3 / D^2;

    d2e_i_dx3_dx4 = -2 * N * bi^3 / D^3;

    % Construct the Jacobian matrix of e_i (4x1 vector of first derivatives)
    J_e_i = [de_i_dx1; de_i_dx2; de_i_dx3; de_i_dx4];

    % Construct the Hessian of e_i (4x4 matrix of second derivatives)
    H_e_i = zeros(4,4);
    H_e_i(1,1) = d2e_i_dx1_dx1;
    H_e_i(2,2) = d2e_i_dx2_dx2;
    H_e_i(3,3) = d2e_i_dx3_dx3;
    H_e_i(4,4) = d2e_i_dx4_dx4;

    H_e_i(1,2) = d2e_i_dx1_dx2; H_e_i(2,1) = d2e_i_dx1_dx2;
    H_e_i(1,3) = d2e_i_dx1_dx3; H_e_i(3,1) = d2e_i_dx1_dx3;
    H_e_i(1,4) = d2e_i_dx1_dx4; H_e_i(4,1) = d2e_i_dx1_dx4;

    H_e_i(2,3) = d2e_i_dx2_dx3; H_e_i(3,2) = d2e_i_dx2_dx3;
    H_e_i(2,4) = d2e_i_dx2_dx4; H_e_i(4,2) = d2e_i_dx2_dx4;

    H_e_i(3,4) = d2e_i_dx3_dx4; H_e_i(4,3) = d2e_i_dx3_dx4;

    % Add contribution to the total Hessian (2 * (J_e_i * J_e_i' + e_i * H_e_i))
    H = H + 2 * (J_e_i * J_e_i' + e_i * H_e_i);
end

end
