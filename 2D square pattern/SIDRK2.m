function u1 = SIDRK2(Spara,tau,a,u0)
% single diagonal implicit RK method

%% ------ mesh setting ------

M = Spara.M;
Dx = Spara.Dx;
Dy = Spara.Dy;
DLap = Spara.DLap;
DDLap = Spara.DDLap;
DDDLap = Spara.DDDLap;

gam = (3+sqrt(3))/6;
a11 = gam; 
a21 = 1-2*gam; a22 = gam;
b1 = 1/2; b2 = 1/2;
% c1 = gam; c2 = 1-gam;

%% ------ solve u1 ------
Tol = 1e-12;
% ------ solve u1 ------
tol = 1;
A = 1/tau/a11 - a*DLap - 2*DDLap - DDDLap;
F1 = fft2(u0)/tau/a11;
u_temp1 = u0;
while tol > Tol
    fftu_temp = fft2(u_temp1(1:M,1:M));
    gradient_x = ifft2(Dx.*fftu_temp);
    gradient_y = ifft2(Dy.*fftu_temp);
    gradient_norm = gradient_x.^2 + gradient_y.^2;
    norm_times_gradient_x = gradient_norm.*gradient_x;
    norm_times_gradient_y = gradient_norm.*gradient_y;
    gradient = ifft2(Dx.*fft2(norm_times_gradient_x(1:M,1:M))) + ...
        ifft2(Dy.*fft2(norm_times_gradient_y(1:M,1:M)));

    Fnl = -DLap.*fft2(gradient);
    F = F1 + Fnl;
    u_temp2 = ifft2(F./A);
    tol = norm(u_temp1-u_temp2,inf);
    u_temp1 = u_temp2;
end
k1h = u_temp1;
k1 = (k1h - u0)/tau/a11;

tol = 1;
A = 1/tau/a22 - a*DLap - 2*DDLap - DDDLap;

F1 = fft2(u0+tau*a21*k1)/tau/a22;
u_temp1 = u0;
while tol > Tol
    fftu_temp = fft2(u_temp1(1:M,1:M));
    gradient_x = ifft2(Dx.*fftu_temp);
    gradient_y = ifft2(Dy.*fftu_temp);
    gradient_norm = gradient_x.^2 + gradient_y.^2;
    norm_times_gradient_x = gradient_norm.*gradient_x;
    norm_times_gradient_y = gradient_norm.*gradient_y;
    gradient = ifft2(Dx.*fft2(norm_times_gradient_x(1:M,1:M))) + ...
        ifft2(Dy.*fft2(norm_times_gradient_y(1:M,1:M)));
    Fnl = -DLap.*fft2(gradient);
    F = F1 + Fnl;
    u_temp2 = ifft2(F./A);
    tol = norm(u_temp1-u_temp2,inf);
    u_temp1 = u_temp2;
end
k2h = u_temp1;
k2 = (k2h - (u0 + tau*a21*k1))/tau/a22;
u1 = u0 + tau*(b1*k1 + b2*k2);
u1(M+1,:) = u1(1,:);
u1(:,M+1) = u1(:,1);

end