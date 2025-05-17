function y = func(x)
    % Kowalik function
    % x is a 4x1 vector [x1; x2; x3; x4]

    % Coefficients (given in literature for Kowalik function)
    a = [0.1957; 0.1947; 0.1735; 0.1600; 0.0844; 0.0627; 0.0456; 0.0342; 0.0323; 0.0235; 0.0246];

    b = [0.25; 0.5; 1.0; 2.0; 4.0; 6.0; 8.0; 10.0; 12.0; 14.0; 16.0];

    x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);

    y = 0;
    for i = 1:11
        numerator = x1 * (1 + x2 * b(i));
        denominator = 1 + x3 * b(i) + x4 * b(i)^2;
        y = y + (a(i) - numerator / denominator)^2;
    end
end
