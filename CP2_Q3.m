% Problem 1 (Make the Matrix)
demoninator = [linspace(1,21,21);linspace(2,42,21);linspace(3,63,21);linspace(4,84,21);
       linspace(5,105,21);linspace(6,126,21);linspace(7,147,21);linspace(8,168,21);
       linspace(9,189,21);linspace(10,210,21);linspace(11,231,21);linspace(12,252,21);
       linspace(13,273,21);linspace(14,294,21);linspace(15,315,21);linspace(16,336,21);
       linspace(17,357,21);linspace(18,378,21);linspace(19,399,21);linspace(20,420,21);];
A = 1./demoninator;

% 1)
A1 = A;
% 2)
A2 = A;
A2(:,16)=linspace(0,0,20);
A2(15,:)=linspace(0,0,21);
% 3)
A3 = A2(18:20,17:21);
% 4)
A4 = A2(:,10);

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


% Problem 3 (Logistic Map)

% 1)
r = 2.75;
% x at 0 is 2 (it's an array cause we need to store multiple values in it)
A11 = [0.2];
% x`n+1 = rx`n(1-x`n)
% They asked us to iterate 100 times 
times_iterated = 99;
% n = 99
for n = 1:times_iterated
    % x`n+1 = rx`n(1-x`n) 
    % Ignore that it changes size each iteration of the for loop
    % (It has to change each time cause we are adding to it)
    A11 = [A11,[r*A11(n)*(1-A11(n))]];
end
% Make the x-axis values
x = 1:(times_iterated + 1);
% Graph it (Just uncomment if you wanna see the graph)
% plot(x,A11)
% The behavior looks like a sink
behaviour = 1;
x_std = std(A11);
% Save as row vector
A12 = [x_std,behaviour];


% 2)
r = 3.25;
% x at 0 is 2
A13 = [0.2];
% x`n+1 = rx`n(1-x`n)
% n = 99
for n = 1:times_iterated
    % x`n+1 = rx`n(1-x`n) where n = c
    A13 = [A13,[r*A13(n)*(1-A13(n))]];
end
% Graph it (Just uncomment if you wanna see the graph)
% plot(x,A13)
% The behavior is a periodic orbit
behaviour2 = 2;
x_std = std(A13);
% Save as row vector
A14 = [x_std,behaviour2];


% 3)
r = 3.75;
% x at 0 is 2
A15 = [0.2];
% x`n+1 = rx`n(1-x`n)
% n = 99
for n = 1:times_iterated
    % x`n+1 = rx`n(1-x`n) where n = c
    A15 = [A15,[r*A15(n)*(1-A15(n))]];
end
% Graph it (Just uncomment if you wanna see the graph)
% plot(x,A15)
% The behavior is chaotic (Looks a bit periodic at first but not really if
% you look closer)
behaviour3 = 3;
x_std = std(A15);
% Save as row vector
A16 = [x_std,behaviour3];