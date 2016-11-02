% Solving 2D Heat Transfer problem with ADI algorithm:
clear all
clc
%% Defining variables:
t_f = 100;                % Time for last data collection
dt = 1;                % time descretization
dx = 0.05;                % x-space descretization
dy = 0.05;                % y-space descretization
t = 0:dt:t_f;            % time variable
x = 0:dx:1.5;            % space intervals in x
y = 0:dy:0.5;            % space intervals in y
N = round(t_f/dt);       % number of timesteps 
I = round(1.5/dx);       % number of space-steps in x
J = round(0.5/dy);       % number of space-steps in y
T(1:I+1,1:J+1,1:N+1) = 0;      % Initialiazing Temperature matrix
T_dummy = T(:,:,1);      % Initializing dummy time variable

%% Initial BCs:
% Parameters for temperature variation, x=0, t>0
a1 = 89.09;
b1 = -0.09036;
c1 = 92.93;
d1 = -0.002168;

% Parameters for temperature variation, x=1.5 in, t>0
a2 = 80.35;
b2 = -0.1156;
c2 = 93.69;
d2 = -0.002442;

% Parabolic fit (2D) - constant along y-direction (for now):
B = a1+c1;
A = ((a2+c2) - B)/(1.5^2);

for j = 1:J+1
    T(1,j,1:N+1) = a1.*exp(b1.*t) + c1.*exp(d1.*t);
    T(I+1,j,1:N+1) = a2.*exp(b2.*t) + c2.*exp(d2.*t);
    T(1:I+1,j,1) = A.*x.^2 + B;
end

B = T(1,1,1:N+1);
A = (T(I+1,1,1:N+1)-B)/(1.5^2);
for n = 1:N+1
    T(:,1,n) = A(n)*x.^2 + B(n);
    T(:,J+1,n) = A(n)*x.^2 + B(n);
end

%% Evolution of T(temp) with time:
a = 1; % 10*1e-6 * (39.37)^2 * 60;
nux = a*dt/dx^2;
nuy = a*dt/dy^2;

Ax = full(gallery('tridiag',I-1,-nux/2,1+nux,-nux/2)); % Defining tridiagonal matrix
Ay = full(gallery('tridiag',J-1,-nuy/2,1+nuy,-nuy/2)); % Defining tridiagonal matrix

for n = 1:N
    T_dummy = (T(:,:,n)+T(:,:,n+1))/2;  % NOTE: Boundary values for T_dummy 
    ... should technically be interpolated to half time step:(T_n+T_n+1)/2
    for j = 2:J    
        rhs1 = (1-nuy).*T(2:I,j,n) + (nuy/2).*(T(2:I,j-1,n)+T(2:I,j+1,n));
        rhs1(1) = rhs1(1)+(nux/2)*T_dummy(1,j);
        rhs1(I-1) = rhs1(I-1)+(nux/2)*T_dummy(I+1,j);
        T_dummy(2:I,j) = Ax\rhs1;
    end 
    
    for i = 2:I
        rhs2 = (1-nux).*T_dummy(i,2:J) + (nux/2).*(T_dummy(i+1,2:J)+T_dummy(i-1,2:J));
        rhs2(1) = rhs2(1)+(nuy/2)*T(i,1,n+1);
        rhs2(J-1) = rhs2(J-1)+(nuy/2)*T(i,J+1,n+1);
        T(i,2:J,n+1) = Ay\rhs2';
    end
    surf(y,x,T(:,:,n),'EdgeColor','none')
    colorbar;
    caxis([75 175]);
    axis([0, 0.5, 0, 1.5, 75, 175]);
    drawnow;
    pause(.01);
end