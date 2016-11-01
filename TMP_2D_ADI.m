% Solving 2D Heat Transfer problem with ADI algorithm:
clear all
clc
%% Defining variables:
t_f = 80;                % Time for last data collection
dt = 0.001;                % time descretization
dx = 0.1;                % x-space descretization
dy = 0.1;                % y-space descretization
t = 0:dt:t_f;            % time variable
x = 0:dx:1.5;            % space intervals
y = 0:dy:0.5;            % space intervals
N = round(t_f/dt);       % number of timesteps 
J = round(1.5/dx);       % number of space-steps
K = round(0.5/dy);       % number of space-steps
T(1:J+1,1:K+1,1:N+1) = 0;      % Initialiazing Temperature matrix
T_dummy = T(:,:,1);      % Initializing dummy time variable

%% Initial BCs:
T(1:J+1,1:K+1,1) = 174.2;

%Geeta's data:
for k = 1:K+1
    T(1,k,1:N+1) = 81.09.*exp(-0.09036.*t)+92.93.*exp(-0.002168.*t);  % variation with time at x=0 in (Interior)
    T(J+1,k,1:N+1) = 80.35.*exp(-0.1156.*t)+93.69.*exp(-0.002442.*t); % variation with time at x=1.5 in (Exterior)
end

for j = 2:J+1
    T(j,1,1:N+1) = T(1,1,1:N+1);
    T(j,K+1,1:N+1) = T(1,K+1,1:N+1);
end

%% Evolution of T(temp) with time:

nu_x = dt/dx^2;
nu_y = dt/dy^2;

Ax = full(gallery('tridiag',J-1,-nu_x/2,1+nu_x,-nu_x/2));  % Defining tridiagonal matrix
Ay = full(gallery('tridiag',K-1,-nu_y/2,1+nu_y,-nu_y/2));  % Defining tridiagonal matrix

for n = 1:100           % Time loop
    %Populating grid for t+1/2
%     n=2;
    for k = 2:K
        for j1 = 2:J
            rhs_x(j1-1) = (1-nu_y)*T(j1,k,n)+(nu_y/2)*T(j1,k-1,n)+(nu_y/2)*T(j1,k+1,n);
        end
        T_dummy(2:J,k) = Ax\rhs_x';
        T_dummy(1,k) = T(1,k,n);
        T_dummy(J+1,k) = T(J,k,n);
    end
    for j = 2:J
        for k1 = 2:K
            rhs_y(k1-1) = (1+nu_x)*T_dummy(j,k1-1)-(nu_x/2)*(T_dummy(j-1,k1-1)+T_dummy(j+1,k1-1));
        end
        T(j,2:K,n+1) = Ay\rhs_y';
    end
end