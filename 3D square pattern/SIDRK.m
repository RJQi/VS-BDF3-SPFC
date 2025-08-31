function u1 = SIDRK(Spara,tau,a,u0)
% single diagonal implicit RK method

%% ------ mesh setting ------

M = Spara.M;
Dx = Spara.Dx;
Dy = Spara.Dy;
Dz = Spara.Dz;
DLap = Spara.DLap;
DDLap = Spara.DDLap;
DDDLap = Spara.DDDLap;

gam = (3-sqrt(3))/6;
a11 = gam; 
a21 = 1-2*gam; a22 = gam;
b1 = 1/2; b2 = 1/2;
% c1 = gam; c2 = 1-gam;

%% ------ solve u1 ------
Tol = 1e-13;
% ------ solve u1 ------
tol = 1;
A = 1/tau/a11 - a*DLap - 2*DDLap - DDDLap;
F1 = fftn(u0)/tau/a11;
u_temp1 = u0;
while tol > Tol
    fftu_temp = fftn(u_temp1(1:M,1:M,1:M));
    gradient_x = ifftn(Dx.*fftu_temp);
    gradient_y = ifftn(Dy.*fftu_temp);
    gradient_z = ifftn(Dz.*fftu_temp);
    gradient_norm = gradient_x.^2 + gradient_y.^2 + gradient_z.^2;
    norm_times_gradient_x = gradient_norm.*gradient_x;
    norm_times_gradient_y = gradient_norm.*gradient_y;
    norm_times_gradient_z = gradient_norm.*gradient_z;
    gradient = ifftn(Dx.*fftn(norm_times_gradient_x(1:M,1:M,1:M))) + ...
        ifftn(Dy.*fftn(norm_times_gradient_y(1:M,1:M,1:M))) + ...
        ifftn(Dz.*fftn(norm_times_gradient_z(1:M,1:M,1:M)));

    Fnl = -DLap.*fftn(gradient);
    F = F1 + Fnl;
    u_temp2 = ifftn(F./A);
%     tol = pagenorm(u_temp1-u_temp2,inf);
    tol = max(max(max(abs(u_temp1-u_temp2))));
    u_temp1 = u_temp2;
end
k1h = u_temp1;
k1 = (k1h - u0)/tau/a11;

tol = 1;
A = 1/tau/a22 - a*DLap - 2*DDLap - DDDLap;

F1 = fftn(u0+tau*a21*k1)/tau/a22;
u_temp1 = u0;
while tol > Tol
    fftu_temp = fftn(u_temp1(1:M,1:M,1:M));
    gradient_x = ifftn(Dx.*fftu_temp);
    gradient_y = ifftn(Dy.*fftu_temp);
    gradient_z = ifftn(Dz.*fftu_temp);
    gradient_norm = gradient_x.^2 + gradient_y.^2 + gradient_z.^2;
    norm_times_gradient_x = gradient_norm.*gradient_x;
    norm_times_gradient_y = gradient_norm.*gradient_y;
    norm_times_gradient_z = gradient_norm.*gradient_z;
    gradient = ifftn(Dx.*fftn(norm_times_gradient_x(1:M,1:M,1:M))) + ...
        ifftn(Dy.*fftn(norm_times_gradient_y(1:M,1:M,1:M))) + ...
        ifftn(Dz.*fftn(norm_times_gradient_z(1:M,1:M,1:M)));
    Fnl = -DLap.*fftn(gradient);
    F = F1 + Fnl;
    u_temp2 = ifftn(F./A);
%     tol = pagenorm(u_temp1-u_temp2,inf);
    tol = max(max(max(abs(u_temp1-u_temp2))));
    u_temp1 = u_temp2;
end
k2h = u_temp1;
k2 = (k2h - (u0 + tau*a21*k1))/tau/a22;
u1 = u0 + tau*(b1*k1 + b2*k2);
u1(M+1,:,:) = u1(1,:,:);
u1(:,M+1,:) = u1(:,1,:);
u1(:,:,M+1) = u1(:,:,1);


end