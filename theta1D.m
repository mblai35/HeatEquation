clear; close all; clc;

% Initial settings
len   = 3 ; % inches
T     = 70; % minutes
alpha = 1 %10*1e-6 * (39.37)^2 * 60; % inch/min (Thermal Diffusivity)


% Boundary Model Function: f(t) = a1 * exp(a2*x) + a3 * exp(a4*x)
par_cnt = [85.6  -.09233 90.35 -.002332];   % Boundary model function parameters
par_bnd = [63.76 -.2109  96.21 -.003575];   % Center   model function parameters

BndMod  = @(t) par_bnd(1) * exp( par_bnd(2) * t ) + ...
               par_bnd(3) * exp( par_bnd(4) * t );
           
CntMod  = @(t) par_cnt(1) * exp( par_cnt(2) * t ) + ...
               par_cnt(3) * exp( par_cnt(4) * t );
        
           
% Initial Data are linear interpolation of CntMod(0) and BndMod(0)
IniMod  = @(x) (BndMod(0)-CntMod(0))*4/(len^2)*(x-len/2).^2+CntMod(0);

           
% Initializing

dx = .2; dt = .01;
theta = 1/2-(dx)^2/(12*dt); 

% Check stability
mu = dt/(dx^2)*alpha;
if ( mu * (1-theta) > 1/2 )
    fprintf('Not stable!');
    return;
end

x = 0:dx:len;
t = 0:dt:T;
U = IniMod(x);


% Parameters for Thomas Algorithm
a = mu * theta;
b = 2*mu*theta + 1;
c = a;
D = @(U) [nan, U(1:end-2)*mu*(1-theta)+U(2:end-1)*(1-2*mu*(1-theta))+U(3:end)*mu*(1-theta), nan];
e = zeros(size(x));
f = e;
Unew = zeros(size(U));
plot(x,U);
axis([0,3,77,180]);
ax = axes;


% Calculation
for it = 2 : numel(t)
    d = D(U);
    Unew(1)   = BndMod(t(it));
    Unew(end) = Unew(1);
    f(1)      = Unew(1);
    for ix = 2 : numel(x)-1
        e(ix) = c / (b - a*e(ix-1));
        f(ix) = (d(ix)+a*f(ix-1)) / (b - a*e(ix-1));
    end
    
    for ix = numel(x)-1 : -1 : 2
        Unew(ix) = f(ix) + e(ix) * Unew(ix+1);
    end
    plot(x,Unew);
    axis([0,3,77,180]);
    drawnow;
    pause(.01);
    U=Unew;
end


