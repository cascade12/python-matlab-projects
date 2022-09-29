% Define the given fuction w tag time, t
c = @(t) 1.3*(exp(-t/11)-exp(-4*t/3));

% Define the derivative of the given fuction
c_dv = @(t) 1.3*(-exp(-t/11)/11+4*exp(-4*t/3)/3);

% Define a time vector that steps by 0.5
t = 1:0.5:20; 

% After plotting we can find the local maxima
min = 1;
max = 3;

% find the root of c_dv

error_tolerance = 10^-8;
root = bisec(c_dv, max, min, error_tolerance);

A1 = root;
A2 = c(root);
A3 = abs(c_dv(A1));




% Static vars
error_tolerance = 10^-8;
x0 = 2;
max_iteration = 2500;

% Declare all fuctions
f = @(x) x^2;
f_dx = @(x) 2 * x;
f_2 = @(x) x^500;
f_dx_2 = @(x) 500 * (x^499);
f_3 = @(x) x^1000;
f_dx_3 = @(x) 1000 * (x^999);

[iter_num, root] = newtons_method(f, f_dx, x0, max_iteration, error_tolerance);
A4 = iter_num;
A5 = root;
[iter_num1, root1] = newtons_method(f_2, f_dx_2, x0, max_iteration, error_tolerance);
A6 = iter_num1;
A7 = root1;
[iter_num2, root2] = newtons_method(f_3, f_dx_3, x0, max_iteration, error_tolerance);
A8 = iter_num2;
A9 = root2;

% Declare bisection method
function root = bisec(f ,a , b, tolerance)
% Inputs: derivative function = f, max = a, min = b, 
% tolerance = error tolerance

% Outputs: midpoint = root including the tolerance


    if f(a) * f(b) >= 0
        error("f(a) * f(b) has to be less than zero")
    end

    % Keep looping til we find it
    while true
        % Find mid point
        root = (a + b)/2; 
        % Get value at midpoint
        root_val = f(root); 
        
        % check for tolerance        
        if abs(root_val) < tolerance
            % Leave loop if found
            break; 
               
        else
            % check if the sign of fval is equal to 
            % the sign of the function value at min
            if sign(root_val) == sign(f(b)) 
                % Change lower limit
                b = root; 
            else 
                % Change upper limit
                a = root;
            end
        end
        
    end
    
end

%Matlab function for Newton Method
function [iter_num, root] = newtons_method(func, func_dx, x0, max_iteration, error_tolerance)
% Inputs: 
% func = fuction we take in
% func_dx = derivative of the fuction we take in
% x0 = inital guess
% max_iteration = maximum iterations
% error_tolerance = error tolerance

% Outputs:
% iter_num = iteration number
% root = final iteration

xx = x0;            % initial guess

% Loop for all intial guesses
    for iter_num = 1 : max_iteration
        % Newton Raphson Formula
        x2 = double(xx - (func(xx)./func_dx(xx))); 
        rt(iter_num + 1) = x2;
        cc = abs(rt(iter_num) - rt(iter_num + 1));
        xx = x2;
        if cc == 0
            break
        end
        %
    end
    root = xx;
end