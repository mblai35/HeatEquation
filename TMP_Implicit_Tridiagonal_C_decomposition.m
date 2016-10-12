clear all
clc
%% Define variables:
% Temperature - T is a matrix where each row is reserved for spatial 
% (x - position) variable and column is reserved for a temporal (t - time) 
% variable 
t_f = 80;                % Time for last data collection
dt = 0.1;                % time descretization
dx = 0.1;                % space descretization
t = 0:dt:t_f;            % time variable
x = 0:dx:1.5;            % space intervals
N = round(t_f/dt);       % number of timesteps 
J = round(1.5/dx);       % number of space-steps
T(1:J+1,1:N+1) = 0;      % Initialiazing Temperature matrix

%% Boundary conditions:
T(1:J+1,1) =  %163.5;%174.03;        % at time = 0 min, temp is in F
%Geeta's data:
% T(1,1:N+1) = 81.09.*exp(-0.09036.*t)+92.93.*exp(-0.002168.*t);  % variation with time at x=0 in (Interior)
% T(J+1,1:N+1) = 80.35.*exp(-0.1156.*t)+93.69.*exp(-0.002442.*t); % variation with time at x=1.5 in (Exterior)

%Xiukun's data:
% T(1,1:N+1) = 85.6.*exp(-0.09233.*t)+90.35.*exp(-0.002332.*t);  % variation with time at x=0 in (Interior)
% T(J+1,1:N+1) = 63.76.*exp(-0.2109.*t)+96.21.*exp(-0.003575.*t); % variation with time at x=1.5 in (Exterior)

%Mallory's data:
T(1,1:N+1) = 90.49.*exp(-0.06361.*t)+80.13.*exp(-0.001023.*t);  % variation with time at x=0 in (Interior)
T(J+1,1:N+1) = 73.97.*exp(-0.08249.*t)+81.56.*exp(-0.001303.*t); % variation with time at x=1.5 in (Exterior)

%% Solving the tridiagonal system of equations:
% System of equations:
% -nu*U_j-1 + (1+2*nu)*U_j - nu*U_j+1 = U_j^n+1, j = 2,3,...,J
% [A]{U} = {U^0}

nu = dt/dx^2;
A = full(gallery('tridiag',J-1,-nu,1+2*nu,-nu)); % Defining tridiagonal matrix
[l,m] = size(A);
U = zeros(l);
% Cholesky decomposition:
for i = 1:l
       sum = 0;
       for j = 1:i-1
            sum = sum + U(j,i)^2;
       end     
       U(i,i) = sqrt(A(i,i) - sum);
       
       sum2 = 0;
       for n = i+1:l
           for k = 1:i     
            sum2 = sum2 + U(k,i)*U(k,n);
           end
       U(i,n) = (A(i,n) - sum2)/U(i,i);
       end
end
L = U' ;

for k = 2:(N+1)
    T_rhs = T(2:J,k-1);                          % Declaring RHS of the system of equations
    T_rhs(1) = T_rhs(1) + nu*T(1,k);
    T_rhs(J-1) = T_rhs(J-1) + nu*T(J+1,k);
    b = T_rhs;
    % Calculating  {s} in [U]{s} = {b}
    s = zeros(l,1);
    
    for h = 1:l
        sum3 = 0;
        for o = 1:h-1
            sum3 = sum3 + L(h,o)*s(o);
        end
        s(h) = (b(h) - sum3)/L(h,h);
    end

    % Calculating {x} in [L]{x}={s}
    x = zeros(l,1);
    for p = l:-1:1
        sum4 = 0;
        for w = p+1:l
            sum4 = sum4 + L(w,p)*x(w);
        end
        x(p) = (s(p) - sum4)/U(p,p);
    end
    
    T(2:J,k) = x;
end

surf(T,'EdgeColor','none')