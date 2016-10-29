clear; close all; clc;
cd ~/HeatEquation/HeatEquation/

% Initial settings
length = 3 ; % inches
height = .5; % inches
T     = 70; % minutes
alpha = 10*1e-6 * (39.37)^2 * 60; % inch^2/min (Thermal Diffusivity)


% Boundary Model Function: f(t) = a1 * exp(a2*x) + a3 * exp(a4*x)
par_cnt = [85.6  -.09233 90.35 -.002332];   % Boundary model function parameters
par_bnd = [63.76 -.2109  96.21 -.003575];   % Center   model function parameters

BndMod  = @(t) par_bnd(1) * exp( par_bnd(2) * t ) + ...
               par_bnd(3) * exp( par_bnd(4) * t );
           
CntMod  = @(t) par_cnt(1) * exp( par_cnt(2) * t ) + ...
               par_cnt(3) * exp( par_cnt(4) * t );
        
           
% Initializing
t0 = -1.97;
dx = .01; 
dy = .01;
dt = .001;
dth = dt/2;
thetax = 1; 
thetay = 1;


% Initial Data are linear interpolation of CntMod(0) and BndMod(0)
IniMod  = @(x,y) BndMod(t0)*ones(size(x));


% Check stability
mux = dt/(dx^2)*alpha;
muy = dt/(dy^2)*alpha;
conx = mux * (1-thetax);
cony = muy * (1-thetay);
if ( conx > 1/2 || cony > 1/2 )
    fprintf('Not stable!');
    return;
end

x = 0:dx:length;
y = 0:dy:height;
[X,Y]=meshgrid(x,y);
t = t0:dt/2:T;
U = IniMod(X,Y);


% Parameters for Thomas Algorithm
ax = mux * thetax;
bx = 2*mux*thetax + 1;
cx = ax;
ay = muy * thetay;
by = 2*muy*thetay + 1;
cy = ay;
Dx = @(U) [U(1,:)*nan; U(1:end-2,:)*cony+U(2:end-1,:)*(1-2*cony)+U(3:end,:)*cony; U(1,:)*nan];
Dy = @(U) [U(:,1)*nan, U(:,1:end-2)*conx+U(:,2:end-1)*(1-2*conx)+U(:,3:end)*conx, U(:,1)*nan];
ex = zeros(size(X));
fx = ex;
ey = ex;
fy = ey;
Unew = zeros(size(U));
mesh(X,Y,U);
axis([0, length, 0, height, 180, 195]);
colorbar;
drawnow;
% axis([0,3,77,180]);


% Calculation
it = 2;
while (it <= numel(t))
    % First half time step
    d = Dx(U);
    Unew(:,1)   = BndMod(t(it));
    Unew(:,end) = Unew(:,1);
    fx(:,1)     = Unew(:,1);
    
    for ix = 2 : numel(x)-1
        ex(2:end-1,ix) = cx ./ (bx - ax*ex(2:end-1,ix-1));
        fx(2:end-1,ix) = (d(2:end-1,ix)+ax*fx(2:end-1,ix-1)) ./ (bx - ax*ex(2:end-1,ix-1));
    end
    
    for ix = numel(x)-1 : -1 : 2
        Unew(2:end-1,ix) = fx(2:end-1,ix) + ex(2:end-1,ix) .* Unew(2:end-1,ix+1);
    end
    U = Unew;
    it = it + 1;
    
    % Next half time step
    d = Dy(U);
    Unew(1,:) = BndMod(t(it));
    Unew(end,:) = Unew(1,:);
    fy(1,:) = Unew(1,:);
    
    for iy = 2 : numel(y)-1
        ey(iy,2:end-1) = cy ./ (by - ay * ey(iy-1,2:end-1));
        fy(iy,2:end-1) = (d(iy,2:end-1) + ax*fy(iy-1,2:end-1)) ./ (by - ay * ey(iy-1,2:end-1));
    end
    
    for iy = numel(y)-1 : -1 : 2
        Unew(iy,2:end-1) = fy(iy,2:end-1) + ey(iy,2:end-1) .* Unew(iy+1,2:end-1);
    end
    U = Unew;
    it = it + 1;
    mesh(X,Y,U);
    colorbar;
    axis([0, length, 0, height, 180, 195]);
    drawnow;
    pause(.01);
end


