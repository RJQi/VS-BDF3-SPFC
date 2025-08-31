function Para = spacemesh(L,M)
% Generate space mesh for Fourier spectral method
h = L/M;
x = 0:h:L-h; 
y = x;
z = x;
[X,Y,Z] = meshgrid(x,y,z);
dx = 1i*2*pi/L*[0:M/2-1,0,-M/2+1:-1];
dy = dx;
dz = dy;
[Dx,Dy,Dz] = meshgrid(dx,dy,dz);
DDx = Dx.^2; DDy = Dy.^2; DDz = Dz.^2;
DLap = DDx + DDy + DDz;
DDLap = DLap.^2;
DDDLap = DLap.^3;

Para.M = M;
Para.h = h;
Para.X = X;
Para.Y = Y;
Para.Z = Z;
Para.Dx = Dx;
Para.Dy = Dy;
Para.Dz = Dz;
Para.DLap = DLap;
Para.DDLap = DDLap;
Para.DDDLap = DDDLap;
end