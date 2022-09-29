
%%%  Problem 1
%%%  First go to the end of your m-file to create your Jacobi and
%%%  Gauss-Seidel functions.

%%% Once you have created your functions return here.
%%% Initialize your matrix A and RHS b
    A=[1.1,0.2,-0.2,0.5;0.2,0.9,0.5,0.3;0.1,0,1,0.4;0.1,0.1,0.1,1.2];
    b=[1;0;1;0];

%%% Use your Jacobi and Gauss-Seidel functions to find A1 through A4.
% Different sizes for A1 - A4
    e=10^-2;
    [Tj_2,Ej_2] = Jacobi(A, b, e);
    [Tgs_2,Egs_2]=Gauss_Seidel(A,b,e);
    
    e=10^-4;
    [Tj_4,Ej_4] = Jacobi(A, b, e);
    [Tgs_4,Egs_4]=Gauss_Seidel(A,b,e);
    
    e=10^-6;
    [Tj_6,Ej_6] = Jacobi(A, b, e);
    [Tgs_6,Egs_6]=Gauss_Seidel(A,b,e);
    
    e=10^-8;
    [Tj_8,Ej_8] = Jacobi(A, b, e);
    [Tgs_8,Egs_8]=Gauss_Seidel(A,b,e);
    
    % Jacobi values
    A1=[Tj_2 Tj_4,Tj_6,Tj_8];
    A2=[Ej_2 Ej_4,Ej_6,Ej_8];
    % Gauss-Seidel values
    A3=[Tgs_2 Tgs_4 Tgs_6 Tgs_8 ];
    A4=[Egs_2 Egs_4 Egs_6 Egs_8 ];


%%%  Problem 2
%%%  Initialize your Day 0 vector x
    S = 0.9;
    I = 0.09;
    R = 0.01;
    x = [S;I;R];

      
%%%  Part 1: without a vaccine
%%%  Make sure to have p = 0
%%%  Initialize the SIR matrix M, and save it as A5
    A5 = x.*[(R / 10000) + (S - (S / 200));S/200 + (I - I/1000);I / 1000 + (R - (R / 10000))];




%%%  Create a loop to find the day that the number of infected
%%%  individuals hits 50% and another loop for the steady state of the
%%%  infected population
%%%  There is a way to put everything under one loop if you make clever use
%%%  of conditionals

% When it reaches 50% infected
D0 = 0;
    iteration_number = 10000;
    for p = 0:iteration_number
        if I < .5
            % 1/1000 I to R
            % 1/200 S turn into I
            % 1/10 000 R to S
            new_s = (R / 10000);
            new_i = (S / 200);
            new_r = (I / 1000);
        
            true_s = new_s + (S - new_i);
            true_i = new_i + (I - new_r);
            true_r = new_r + (R - new_s);

            S = true_s;
            I = true_i;
            R = true_r;
        else
            D0 = p;
            break;
        end
    end
% Steady State
iteration_number = 10000;
    for p = 0:iteration_number
            new_s = (R / 10000);
            new_i = (S / 200);
            new_r = (I / 1000);
        
            true_s = new_s + (S - new_i);
            true_i = new_i + (I - new_r);
            true_r = new_r + (R - new_s);

            S = true_s;
            I = true_i;
            R = true_r;
    end
    F0 = I;

%%% Save the days and steady state in a row vector A6


A6 = [D0,F0];


%%%  Reinitialize your Day 0 vector x

    S = 0.9;
    I = 0.09;
    R = 0.01;
    x = [S;I;R];
       

%%%  Part 2: with a vaccine
%%%  Make sure to have p = 2/1000
%%%  Initialize the SIR matrix M, and save it as A7


p = 2/1000;



%%%  Create a loop to find the day that the number of infected
%%%  individuals hits 50% and another loop for the steady state of the
%%%  infected population
%%%  There is a way to put everything under one loop if you make clever use
%%%  of conditionals

% When it reaches 50% infected
D1 = 0;
    iteration_number = 1000;
    for iter = 0:iteration_number
        if I < .5
            % 1/1000 I to R
            % 1/200 S turn into I
            % 1/10 000 R to S
            new_s = (R / 10000) - p;
            new_i = (S / 200);
            new_r = (I / 1000) + p;
        
            true_s = new_s + (S - new_i);
            true_i = new_i + (I - new_r);
            true_r = new_r + (R - new_s);

            S = true_s;
            I = true_i;
            R = true_r;
        else
            D1 = iter;
            break;
        end
    end
% Steady State
iteration_number = 10000;
    for iter = 0:iteration_number
            new_s = (R / 10000) - p;
            new_i = (S / 200);
            new_r = (I / 1000) + p;
        
            true_s = new_s + (S - new_i);
            true_i = new_i + (I - new_r);
            true_r = new_r + (R - new_s);

            S = true_s;
            I = true_i;
            R = true_r;
    end

%%% Save the days and steady state in a row vector A8

A8 = [D1, F1];
 
 
%%%  Problem 3
  
%%%  Initialize your 114x114 tridiagonal matrix A

n = 114; % size of matrix
A = diag(2*ones(1,n))... %Main diagonal with entries 2
-diag(ones(1,n-1),1)... %Upper diagonal with entries -1
-diag(ones(1,n-1),-1); % Lower diagonal with entries -1

A9 = A;

%%%  Initialize your 114x1 RHS column vector rho

rho = 2*(1-cos(53*pi/115))*sin(53*pi*(1:114)/115); % Creating the vector Pj

A10 = rho;

%%%  Implement Jacobi's method for this system.
%%%  Don't use the function you created before because that was designed for
%%%  a specific task, and will not work here.

eps = 10^-5;
[A11, A12, A13] = Jacobi_Q3(A,rho,eps);


%%%  Create a column vector phi that contains the exact solution given in
%%%  the assignment file



%%%  Save the difference of the Jacobi solution and the exact solution as
%%%  A13.  Use the maximal entry in absolute value to calculate this error.




%%%  Implement Gauss-Seidel for this system.
%%%  Don't use the function you created before because that was designed for
%%%  a specific task, and will not work here.



eps = 10^-5;
[A14, A15, A16] = Gauss_Seidel_Q3(A,rho,eps);


%%%  Save the difference of the Gauss-Seidel solution and the exact solution as
%%%  A13.  Use the maximal entry in absolute value to calculate this error.



  

%%% Jacobi and Gauss Seidel Iteration functions
%%% Create your functions here
%%% Both functions will need two outputs and three inputs
%%% The code within the function will be very similar to
%%% Week 4 coding lecture 2

% Function to run Jacobi algorithm
function [T,E] = Jacobi(A, b, eps)
 % Input:
 % A is a matrix
 % b is the right hand side vector in the system Ax =b
 % eps is the tolerance
 %
 % Output:
 % T is the iteration at which the largest entry of abs(Ay-b)<eps
 % E is the largest entry in abs(Ay-b)
 
 % Start with error checking
 % Check if b is a row vector
 [row,~] = size(b);
 if row == 1
    b = b'; % convert to column vector
 end
 
 % Check to see if its a square
 [row,col] = size(A);
 if row ~= col
     error("A must be a square matrix")
 end
 
 
% set diagonals to 0
D = diag(A);
remove_diagonals = zeros(size(A));
for i = 1:length(b)
    remove_diagonals(i,i) = D(i);   
end
R = A - remove_diagonals;

% Now the the unworking matrixes are removed, start a method

% Step 1: y = zeros(n, 1)
y = zeros(row,1);

% Iteration Number
interation_number = 200; 

% Initialize storage to zero
E = 0;

% Step 2: Find for each T = ...
% Find the entry number
for t = 1:interation_number
    % Step 3: Update according to Jacobi
      y = (b - R * y).*(D.^-1);
      % Step 4: Break when max|Ay-b| < E
      if max(abs(A * y - b)) < eps
          % Store the value as E
          E = max(abs(A * y - b));
          % Store where its broken
          T = t;
          break; 
      end
   
end

end

% Function to run Gauss Seidel algorithm
function [T,E] = Gauss_Seidel(A,b,eps)
     % Inputs: 
     % A is a matrix
     % b is the right hand side vector in the system Ax =b
     % eps is the tolerance
     %
     % Outputs:
     % T is the iteration at which the largest entry of abs(Ay-b)<eps
     % E is the largest entry in abs(Ay-b)
     
     % get the size of the matrix A
    [col,row] = size(A);
    
    % Start with error checking
    % Check if b is a column vector
    [col1,~] = size(b);
    if col1 == 1
        % Convert to column vector
        b = b'; 
    end
    
   % Step 1: Set y = zeros(n,1)
   y = zeros(row,1);
    
   % Gauss seidel algorithm implementation
   
   % Extract the upper triangular matrix from A and initalize U
   U = zeros(size(A)); 
   % Fill matrix
   for i = 1:col
     for j = 1:row
       if i ~= j && j > i
         U(i,j) = A(i,j);
       end
     end
   end
   
         
   % Get the lower triangular matrix L by
   % Removing the upper triangular matrix from A
   L = A - U;

   iteration_number = 200; 
    
    % Gauss seidel algorithm implementation
    % Step 2: For each T = ...
    for row = 1:iteration_number
        next_y = L\(b - U * y);
        y = next_y;
        % Check and break
        % Step 4: Break when max|Ay - b| < E
        if max(abs(A * next_y - b)) < eps
            % Store iteration time
            T = row; 
            % Store breaking value
            E = max(abs(A * next_y - b));
            break; 
        end
        
    end
end

function [T,E,max] = Jacobi_Q3(A,rho,eps)
    n = 114;
    % Create variable to store unknowns
    phi = zeros(n,1); % column vector
    % Decomposing the A matrix in to D,L and U
    
    % D is the diagonal matrix of A
    D = diag(diag(A)); 
    % Lower diagonal matrix of A
    L = tril(A,-1); 
    % Upper diagonal matrix of A
    U = triu(A,1); 
    
    % Calculate M and c for later
    M = D \ (-L - U); 
    c = D \ rho'; % Calculating c in the given equation
    
    ext = 1;
    while ext % loop will continue until ext becomes 0
        phi_old = phi; 
        phi = (M * phi + c); 
        % Check accurate w tolerance of 10^-5
        if( norm(phi - phi_old) <= eps) 
            % Store value
            T = phi; 
            % Store number of loop
            E = ext; 
            % Break out
            ext = 0; 
        else 
            % Iterate through matrix
            ext = ext + 1; 
        end
    end
    max = 0;
end

function [T,E,max] = Gauss_Seidel_Q3(A,rho,eps)
    n = 114;
    % Create variable to store unknowns
    phi = zeros(n,1); % column vector
    % Decomposing the A matrix in to D,L and U
    
    % D is the diagonal matrix of A
    D = diag(diag(A)); 
    % Lower diagonal matrix of A
    L = tril(A,-1); 
    % Upper diagonal matrix of A
    U = triu(A,1); 
    
    % Calculate M and c for later
    M = D \ (-L - U); 
    c = D \ rho'; % Calculating c in the given equation
    ext = 1; 
    % Loop til you find it
    while ext 
        % save current values to calculate error later
        phi_old = phi; 
        % Gauss-Seidel iteration
        phi = (M * phi + c); 
        % Check accurate w tolerance of 10^-5
        if( norm(phi - phi_old) <= eps) 
            % Store value
            T = phi; 
            % Store number of loop
            E = ext; 
            % Break out
            ext = 0; 
        else 
            % Iterate through matrix
            ext = ext + 1; 
        end
    end
    max = 0;
end



