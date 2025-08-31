function [t,tau,Energy,Masserror,usnap] = BDF3coarsening(Spara,T,taumin,taumax,rhomax,rhomin,beta,a,u0,snaptime)
%
%

%% ------ mesh setting ------
 
M = Spara.M;
h = Spara.h;
Dx = Spara.Dx;
Dy = Spara.Dy;
Dz = Spara.Dz;
DLap = Spara.DLap;
DDLap = Spara.DDLap;
DDDLap = Spara.DDDLap;
usnap = zeros(M+1,M+1,M+1,length(snaptime));
s = 1;
%% ------ solve scheme ------
a0f = @(x,y) (2*x+1)./(x+1)+x.*y./(x.*y+y+1);
a1f = @(x,y) -x./(x+1)-x.*y./(x.*y+y+1)-x.*y.^2.*(x+1)./(x.*y+y+1)./(y+1);
a2f = @(x,y) x.*y.^2.*(x+1)./(x.*y+y+1)./(y+1);
t(1) = 0; 
tau(1) = taumin;
t(2) = taumin; 
tau(2) = taumin; 
t(3) = 2*taumin;

u(:,:,:,1) = u0;
Mass0 = h^3*sum(sum(sum(u(1:M,1:M,1:M,1)))); 
[Energy(1), Masserror(1)] = EnergyMass(Spara,a,Mass0,u0);

u(:,:,:,2) = SIDRK(Spara,taumin/100,a,u(1:M,1:M,1:M,1));
[Energy(2), Masserror(2)] = EnergyMass(Spara,a,Mass0,u(:,:,:,2));

u(:,:,:,3) = SIDRK(Spara,taumin/100,a,u(1:M,1:M,1:M,2));
[Energy(3), Masserror(3)] = EnergyMass(Spara,a,Mass0,u(:,:,:,3));

t_currunt = t(3);
k = 3;
Tol = 1e-10;
tic;
while t_currunt < T

    Eold = Energy(k-1);
    Enew = Energy(k);
    cphi = (Enew-Eold)/tau(k-1);
    nphi = cphi^2;
    taub = sqrt(1+beta*nphi);
    tauada = max(taumin,taumax/taub);
    tau(k) = max(min(tauada,rhomax*tau(k-1)),rhomin*tau(k-1));

    rho1 = tau(k)/tau(k-1);
    rho2 = tau(k-1)/tau(k-2);
    a0 = a0f(rho1,rho2);
    a1 = a1f(rho1,rho2);
    a2 = a2f(rho1,rho2);
    A = a0/tau(k) - a*DLap - 2*DDLap - DDDLap;
    
    F1 = a0/tau(k)*fftn(u(1:M,1:M,1:M,3)) - a1/tau(k-1)*fftn(u(1:M,1:M,1:M,3)-u(1:M,1:M,1:M,2))...
        - a2/tau(k-2)*fftn(u(1:M,1:M,1:M,2)-u(1:M,1:M,1:M,1));

    tol = 1;
    u_temp1 = u(1:M,1:M,1:M,3);
    while tol > Tol
        fftu_temp = fftn(u_temp1(1:M,1:M,1:M));
        gradient_x = ifftn(Dx.*fftu_temp);
        gradient_y = ifftn(Dy.*fftu_temp);
        gradient_z = ifftn(Dz.*fftu_temp);
        gradient_norm = gradient_x.^2 + gradient_y.^2 + gradient_z.^2;
        norm_times_gradient_x = gradient_norm.*gradient_x;
        norm_times_gradient_y = gradient_norm.*gradient_y;
        norm_times_gradient_z = gradient_norm.*gradient_z;
        gradient = Dx.*fftn(norm_times_gradient_x(1:M,1:M,1:M)) + ...
            Dy.*fftn(norm_times_gradient_y(1:M,1:M,1:M)) + ...
            Dz.*fftn(norm_times_gradient_z(1:M,1:M,1:M));
 
        Fnl = -DLap.*gradient;
        F = F1 + Fnl;
        u_temp2 = ifftn(F./A);
%         tol = pagenorm(u_temp1-u_temp2,inf);
        tol = max(max(max(abs(u_temp1-u_temp2))));
        u_temp1 = u_temp2;
    end
    u(:,:,:,1) = u(:,:,:,2);
    u(:,:,:,2) = u(:,:,:,3);
    u(1:M,1:M,1:M,3) = u_temp1;
    u(M+1,:,:,3) = u(1,:,:,3);
    u(:,M+1,:,3) = u(:,1,:,3);
    u(:,:,M+1,3) = u(:,:,1,3);
    t_currunt = t_currunt + tau(k);
    t(k+1) = t_currunt;

    [Energy(k+1), Masserror(k+1)] = EnergyMass(Spara,a,Mass0,u(:,:,:,3));
    if s <=length(snaptime)
        if t_currunt >= snaptime(s) && t_currunt <= snaptime(s) + taumax
            usnap(:,:,:,s)=u(:,:,:,3);
            s = s+1;
        end
    end
    X = [tau(k) t_currunt Energy(k+1) tol];
    k = k+1;
    disp(X);
end
toc;
end