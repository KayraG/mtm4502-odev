function f = func(x)
% FUNC calculates the objective function value for the Kowalik Problem.
%
%   f = func(x)
%
%   Input:
%       x - A 4x1 column vector representing the current point (x1, x2, x3, x4).
%
%   Output:
%       f - The scalar objective function value.

% Define the constants ai and bi for the Kowalik Problem (from Table 7 in [1])
a_vals = [0.1957, 0.1947, 0.1735, 0.16, 0.0844, 0.0627, 0.0456, 0.0342, 0.0323, 0.0235, 0.0246];
b_vals = [0.25, 0.50, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0];

% Initialize the sum of squares
sum_of_squares = 0;

% Loop through i from 1 to 11
for i = 1:11
    ai = a_vals(i);
    bi = b_vals(i);

    % Calculate the denominator D
    D = 1 + x(3) * bi + x(4) * bi^2;

    % Calculate the numerator N
    N = x(1) * (1 + x(2) * bi);

    % Calculate the term inside the square
    term = ai - (N / D);

    % Add the square of the term to the sum
    sum_of_squares = sum_of_squares + term^2;
end

% The objective function value is the sum of squares
f = sum_of_squares;

end
