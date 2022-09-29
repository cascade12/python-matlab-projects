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