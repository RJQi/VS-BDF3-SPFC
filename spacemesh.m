function Para = spacemesh(L,M)
% Generate space mesh for Fourier spectral method
h = L/M;
x = 0:h:L-h; y = x;
[X,Y] = meshgrid(x,y);
dx = 1i*2*pi/L*[0:M/2-1,0,-M/2+1:-1];
dy = dx;
[Dx,Dy] = meshgrid(dx,dy);
DDx = Dx.^2; DDy = Dy.^2;
DLap = DDx + DDy;
DDLap = DLap.^2;
DDDLap = DLap.^3;

Para.M = M;
Para.h = h;
Para.X = X;
Para.Y = Y;
Para.Dx = Dx;
Para.Dy = Dy;
Para.DLap = DLap;
Para.DDLap = DDLap;
Para.DDDLap = DDDLap;

% h = L/M;
% x = 0:h:L-h; y = x;
% [X,Y] = meshgrid(x,y);
% dx1 = 1i*2*pi/L*[0:M/2-1,0,-M/2+1:-1];
% dy1 = dx1;
% [Dx,Dy] = meshgrid(dx1,dy1);
% dx2 = 1i*2*pi/L*[0:M/2,-M/2+1:-1];
% dy2 = dx2;
% [Dx2,Dy2] = meshgrid(dx2,dy2);
% DDx = Dx2.^2; DDy = Dy2.^2;
% DLap = DDx + DDy;
% DDLap = DLap.^2;
% DDDLap = DLap.^3;
% 
% Para.M = M;
% Para.h = h;
% Para.X = X;
% Para.Y = Y;
% Para.Dx = Dx;
% Para.Dy = Dy;
% Para.DLap = DLap;
% Para.DDLap = DDLap;
% Para.DDDLap = DDDLap;
end