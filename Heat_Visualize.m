function Heat_Visualize( nD, timestep )

mkmovie=1;

    
CntData_Mallory = [171 126 104 93 84 79 77 75 75 73 73];
CntData_Geeta = [174.2 122.9 103.1 93.2 87.8 82.4 78.8];
CntData_Xiukun = [176 122 100 90 84 81 79 77];
tData_Mallory = 0:10:100;
tData_Geeta = [0 10 20 30 40 50 80];
tData_Xiukun = 0:10:70;

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

switch t(end)
    case 70
        CntData=CntData_Xiukun;
        tData = tData_Xiukun;
    case 80
        CntData=CntData_Geeta;
        tData = tData_Geeta;
    case 100
        CntData=CntData_Mallory;
        tData = tData_Mallory;
    otherwise
        CntData=[];
        tData = [];
end

[Data, nData] = fread( file_data, 'double' );

minData = min(Data(:));
maxData = max(Data(:));

if nData ~= nx*ny*nt
    fprintf( 'Wrong number of data.\n' );
    return;
end

Data = reshape( Data, nx*ny, nt )';
indxyData=nx * round((ny+1)/2) - round((nx-1)/2);
Cnt = Data( : , indxyData );
figure(1)
plot(t,Cnt,tData,CntData,'x');
xlabel('Time (min)');
ylabel('Temperature (°F)');
legend('Model','Real');
   

if mkmovie
    vidobj = VideoWriter(['theta',num2str(nD),'D'],'MPEG-4');
    vidobj.FrameRate = t(end)/timestep/10;
    vidobj.Quality = 100;
    open(vidobj);
end


for it = 1:ceil(timestep/dt):nt/8*3
    if ny == 1
        figure(2);
        plot( x, Data(it,:) );
        xlabel('x (inch)');
        ylabel('Temperature (°F)');
        title(['t=',num2str((it-1)*dt),' min']);
        axis([0,3,72,180]);
        drawnow;
        pause(timestep);
    else
        figure(2);
        iData = reshape(Data(it,:), nx, ny)';
        surf(X,Y,iData);
        xlabel('x (inch)'); ylabel('y (inch)'); zlabel('Temperature (°F)');
        title(['t=',num2str((it-1)*dt),' min']);
        axis([x(1),x(end),y(1),y(end),minData,maxData]);
        caxis([72,180]);
        drawnow;
    end
    
    if mkmovie
        writeVideo(vidobj,getframe(gcf));
    end
end

if mkmovie
    close(vidobj);
end






