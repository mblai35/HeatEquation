function Heat_Visialize( nD, timestep )

CntData=[176
122
100
90
84
81
79
77];
tData=0:10:70;

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

if file_info == -1 || file_data == -1
    fprintf( 'Error opening files\n' );
    return;
end

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
t = (0:nt-1) * dt + t0;
[X,Y]=meshgrid(x,y);

[Data, nData] = fread( file_data, 'double' );

if nData ~= nx*ny*nt
    fprintf( 'Wrong number of data.\n' );
    return;
end

Data = reshape( Data, nx*ny, nt )';
indxyData=nx * round((ny+1)/2) - round((nx-1)/2);
Cnt = Data( : , indxyData );
figure(1)
plot(t,Cnt,tData,CntData,'x');
    
for it = 1:floor(timestep/dt):nt
    if ny == 1
        figure(2);
        plot( x, Data(it,:) );
        axis([0,3,77,180]);
        drawnow;
        pause(timestep);
    else
        figure(2);
        iData = reshape(Data(it,:), nx, ny)';
        surf(X,Y,iData);
        axis([0,3,0,.5,77,195]);
        caxis([77,195]);
        drawnow;
        pause(timestep);
    end
end








