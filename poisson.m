% The exact solution of the governing PDE
clear all;

xl= 0; xr= 1; % x domain
yl= 0; yr= 1; % y domain

J= 100; % number of points/divisions in both x- and y-directions

h = (xr-xl)/J; % mesh size

% build up the coefficient matrix
nr = (J-1)^2; % order of the matrix
matA = zeros(nr,nr);

for i=1:nr
    matA(i,i)= -4;
    if i+1 <= nr & mod(i,J-1) ~= 0
        matA(i,i+1)= 1;
    end    
    if i+J-1 <= nr
        matA(i,i+J-1)= 1;
    end
    if i-1 >= 1 & mod(i-1,J-1) ~= 0
        matA(i,i-1)= 1;
    end
    if i-(J-1) >= 1
        matA(i,i-(J-1))= 1;
    end
end

% build up the right-hand side vector
for j=1:J-1
    y(j) = j*h;
    for i=1:J-1
        x(i) = i*h;
        [fij]=feval(@srcF,x(i),y(j)); % evaluate f(xi,yj)
        vecF((J-1)*(j-1)+i)=h^2*fij;
        if i == 1
            vecF((J-1)*(j-1)+i) = vecF((J-1)*(j-1)+i) - feval(@Bound, x(i), y(j));% boundary effect changes
        end
        if j == 1
            vecF((J-1)*(j-1)+i) = vecF((J-1)*(j-1)+i) - feval(@Bound, x(i), y(j));% boundary effect changes
        end
        if i == (J - 1)
            vecF((J-1)*(j-1)+i) = vecF((J-1)*(j-1)+i) - feval(@Bound, x(i), y(j));% boundary effect changes
        end
        if j == (J - 1)
            vecF((J-1)*(j-1)+i) = vecF((J-1)*(j-1)+i) - feval(@Bound, x(i), y(j));% boundary effect changes
        end
        
        [uij]=feval(@exU,x(i),y(j)); % evaluate exact solution
        ExU(i,j)=uij;
    end
end

% solve the system
vecU = matA\vecF'; % vecU is of order (J-1)^2
%vecU is a vector, it is converted into matrix
for j=1:J-1
    for i=1:J-1
        U2d(i,j)=vecU((J-1)*(j-1)+i); % change into 2-D array
    end
end

% Calulation of Root Mean Absolute Error {e_delta}
error = 0;
for i = 1:J-1
    for j = 1: J-1
        error = error + abs(ExU(i,j) - U2d(i,j))^2;
    end
end
e_delta = sqrt(h^2*error);


disp('Root Mean Square Absolute Error = '), e_delta,
figure(1);
surf(x,y,U2d);
title('Numerical solution');
figure(2);
surf(x,y,ExU);
title('Exact solution');

% The function to compute boundary effects
function uu=Bound(x,y)
    uu = x^3 + y^3;
end

% The exact solution of the governing PDE
function uu=exU(x,y)
    uu = x^3 + y^3;
end

% The RHS function f for the governing PDE
function uu=srcF(x,y)
    uu = 6*x + 6*y;
end
