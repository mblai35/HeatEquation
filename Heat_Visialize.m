function Heat_Visialize( nD, timestep )

switch(nargin)
    case 0
        nD = 1;
        timestep = 0.1;
    case 1
        timestep = 0.1;
    case 2
    otherwise
        fprintf( 'Wrong number of input.\n' );
end

filename = ['theta' num2str(nD) 'D'];
file_info = fopen( [filename '.txt'], 'r' );
file_data = fopen( [filename '.bin'], 'r' );

Info=fscanf( file_info, '%d %e', 4 );

ny = 1;
dy = 1;
t0 = 0;

i = 1;
nx = Info(i); i = i+1;
dx = Info(i); i = i+1;

if nD == 2
    Info = [Info; fscanf( file_info, '%d %e %e' )];
    ny = Info(i); i = i+1;
    dy = Info(i); i = i+1;
    t0 = Info(7);
end

nt = Info(i); i = i+1;
dt = Info(i);

if timestep < dt
    timestep = dt;
end

x = (0:nx-1) * dx;
y = (0:ny-1) * dy;
t = (0:nt-1) * dt - t0;
[X,Y]=meshgrid(x,y);

[Data, nData] = fread( file_data, 'double' );

if nData ~= nx*ny*nt
    fprintf( 'Wrong number of data.\n' );
    return;
end

Data = reshape( Data, nx*ny, nt )';

for it = 1:floor(timestep/dt):nt
    figure(1)
    if ny == 1
        figure(1);
        plot( x, Data(it,:) );
        axis([0,3,77,180]);
        drawnow;
        pause(timestep);
    else
        figure(1);
        iData = reshape(Data(it,:),nx,ny)';
        surf(X,Y,iData);
        drawnow;
        pause(timestep);
    end
end





