% Problem 2 (Harmonic Series)
s20 = 20;
s200 = 200;
% 1)
A5 = sum(1./(1:s20));
% 2)
A6 = sum(1./(1:s200));
% 3)
n = 0;
% Not and int bur a double
f = 0.0;
while f <= 10
    variable = 1/n;
    n = n + 1;
    f = f + variable;
end
A7 = n + 1;
A8 = f + variable;

% 4)
n = 0;
% Not and int bur a double
f = 0.0;
while f <= 20
    variable = 1/n;
    n = n + 1;
    f = f + variable;
end
A9 = n + 1;
A10 = f + variable;