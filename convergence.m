% Script file to understand the error and convergence of the approximate
% solution of the the Elliptical PDE
clear all;

% Calculating the log errors for grids 10, 20, 30,.. 100
gridDivisions = [10: 10: 100];
h = gridDivisions.^(-1); % h equals the mesh size i.e DelX or DelY
logdelX = log(h); %log(delX) where delX = 1/grids
for i= 1:10
    Error(i) = PoissonSolver(gridDivisions(i)); % stores log of the error for ith grid
    message = [num2str(i), ': Grid Size = ', num2str(gridDivisions(i)), ' Average Error = ',num2str(Error(i))];
    disp(message);
end

% Log-Log Plot of Error against Mesh size
loglog(h, Error);
grid on;
title('LogLog Plot of Error vs Mesh Size')
xlabel('\DeltaX');
ylabel('e_{\Delta}');

% compute order of error using the slope of log(error) against log(h)
linefit = polyfit(log10(h),log10(Error),1); % fitting a straight line
order = linefit(1);
disp(['Order is = ',num2str(order)]);
% order comes out be 0.9119 or approximately 1


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

% Computes the approximate solution for grid size h
% parameter h: grid size
% returns the error {e_delta}
function error = PoissonSolver(Number)

    xl= 0; xr= 1; % x domain
    yl= 0; yr= 1; % y domain

    J = Number; % number of points in both x- and y-directions

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
                vecF((J-1)*(j-1)+i) = vecF((J-1)*(j-1)+i) - feval(@Bound, x(i), y(j));
            end
            if j == 1
                vecF((J-1)*(j-1)+i) = vecF((J-1)*(j-1)+i) - feval(@Bound, x(i), y(j));
            end
            if i == (J - 1)
                vecF((J-1)*(j-1)+i) = vecF((J-1)*(j-1)+i) - feval(@Bound, x(i), y(j));
            end
            if j == (J - 1)
                vecF((J-1)*(j-1)+i) = vecF((J-1)*(j-1)+i) - feval(@Bound, x(i), y(j));
            end
        
            [uij]=feval(@exU,x(i),y(j)); % evaluate exact solution
            ExU(i,j)=uij;
        end
    end

    % solve the system
    vecU = matA\vecF'; % vecU is of order (J-1)^2
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
    error = sqrt(h^2*error);

end

