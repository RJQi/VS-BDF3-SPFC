function [t,tau,Energy,Masserror,usnap] = CSBDF3coarsening(Spara,T,taumin,taumax,rhomax,rhomin,beta,a,u0,snaptime,Sp,index)
%

%% ------ mesh setting ------
 
M = Spara.M;
h = Spara.h;
Dx = Spara.Dx;
Dy = Spara.Dy;
DLap = Spara.DLap;
DDLap = Spara.DDLap;
DDDLap = Spara.DDDLap;
usnap = zeros(M+1,M+1,length(snaptime));
s = 1;
epsi = 1-a;
taum = taumax;

%% ------ solve scheme ------
a0f = @(x,y) (2*x+1)./(x+1)+x.*y./(x.*y+y+1);
a1f = @(x,y) -x./(x+1)-x.*y./(x.*y+y+1)-x.*y.^2.*(x+1)./(x.*y+y+1)./(y+1);
a2f = @(x,y) x.*y.^2.*(x+1)./(x.*y+y+1)./(y+1);
t(1) = 0; 
tau(1) = taumin;
t(2) = taumin; 
tau(2) = taumin; 
t(3) = 2*taumin;
u(:,:,1) = u0;
Mass0 = h^2*sum(sum(u(1:M,1:M,1)));
[Energy(1), Masserror(1)] = EnergyMass(Spara,a,Mass0,u0);
u(:,:,2) = SIDRK2(Spara,taumin,a,u(1:M,1:M,1));
[Energy(2), Masserror(2)] = EnergyMass(Spara,a,Mass0,u(:,:,2));
u(:,:,3) = SIDRK2(Spara,taumin,a,u(1:M,1:M,2));
[Energy(3), Masserror(3)] = EnergyMass(Spara,a,Mass0,u(:,:,3));


t_currunt = t(3);
k = 3;
Tol = 1e-12;
tic;
while t_currunt < T

    if index == 0
        Eold = Energy(k-1);
        Enew = Energy(k);
        cphi = (Enew-Eold)/tau(k-1);
        nphi = cphi^2;
        taub = sqrt(1+beta*nphi);
        tauada = max(taumin,taumax/taub);
        tau(k) = max(min(tauada,rhomax*tau(k-1)),rhomin*tau(k-1));
    else
        tau(k) = 0.005;
    end
%     taum = tau(k);
  

    rho1 = tau(k)/tau(k-1);
    rho2 = tau(k-1)/tau(k-2);
    a0 = a0f(rho1,rho2);
    a1 = a1f(rho1,rho2);
    a2 = a2f(rho1,rho2);
%     A = a0/tau(k) - a*DLap - 2*DDLap - DDDLap;
    A = a0/tau(k) - DLap - (2-Sp*taum^3*a0/tau(k))*DDLap - DDDLap;
    
%     F1 = a0/tau(k)*fft2(u(1:M,1:M,3)) - a1/tau(k-1)*fft2(u(1:M,1:M,3)-u(1:M,1:M,2))...
%         - a2/tau(k-2)*fft2(u(1:M,1:M,2)-u(1:M,1:M,1));
    if k == 3
        F1 = a0/tau(k)*fft2(u(1:M,1:M,3)) - a1/tau(k-1)*fft2(u(1:M,1:M,3)-u(1:M,1:M,2))...
            - a2/tau(k-2)*fft2(u(1:M,1:M,2)-u(1:M,1:M,1))...
            + Sp*taum^3*a0/tau(k)*DDLap.*fft2(u(1:M,1:M,3));
    elseif k == 4
        F1 = a0/tau(k)*fft2(u(1:M,1:M,3)) - a1/tau(k-1)*fft2(u(1:M,1:M,3)-u(1:M,1:M,2))...
            - a2/tau(k-2)*fft2(u(1:M,1:M,2)-u(1:M,1:M,1))...
            + Sp*taum^3*a0/tau(k)*DDLap.*fft2(u(1:M,1:M,3))...
            - Sp*taum^3*a1/tau(k-1)*DDLap.*fft2(u(1:M,1:M,3)-u(1:M,1:M,2));
    else
        F1 = a0/tau(k)*fft2(u(1:M,1:M,3)) - a1/tau(k-1)*fft2(u(1:M,1:M,3)-u(1:M,1:M,2))...
            - a2/tau(k-2)*fft2(u(1:M,1:M,2)-u(1:M,1:M,1))...
            + Sp*taum^3*a0/tau(k)*DDLap.*fft2(u(1:M,1:M,3))...
            - Sp*taum^3*a1/tau(k-1)*DDLap.*fft2(u(1:M,1:M,3)-u(1:M,1:M,2))...
            - Sp*taum^3*a2/tau(k-2)*DDLap.*fft2(u(1:M,1:M,2)-u(1:M,1:M,1));
    end 
    k1 = 1 + 1/rho1;
    k2 = 1 + 1/rho1 + 1/rho1/rho2;
    c1 = k1*k2/(k1-1)/(k2-1);
    c2 = k2/(k1-1)/(k1-k2);
    c3 = -k1/(k1-k2)/(k2-1);
    u_cs = c1*u(1:M,1:M,3) + c2*u(1:M,1:M,2) + c3*u(1:M,1:M,1);
    F1 = F1 - epsi*DLap.*fft2(u_cs);

    tol = 1;
    u_temp1 = u(1:M,1:M,3);
    while tol > Tol
        fftu_temp = fft2(u_temp1);
        gradient_x = ifft2(Dx.*fftu_temp);
        gradient_y = ifft2(Dy.*fftu_temp);
        gradient_norm = gradient_x.^2 + gradient_y.^2;
        norm_times_gradient_x = gradient_norm.*gradient_x;
        norm_times_gradient_y = gradient_norm.*gradient_y;
        gradient = Dx.*fft2(norm_times_gradient_x(1:M,1:M)) + ...
            Dy.*fft2(norm_times_gradient_y(1:M,1:M));
 
        Fnl = -DLap.*gradient;
        F = F1 + Fnl;
        u_temp2 = ifft2(F./A);
        tol = norm(u_temp1-u_temp2,inf);
        u_temp1 = u_temp2;
    end
    u(:,:,1) = u(:,:,2);
    u(:,:,2) = u(:,:,3);
    u(1:M,1:M,3) = u_temp1;
    u(M+1,:,3) = u(1,:,3);
    u(:,M+1,3) = u(:,1,3);
    t_currunt = t_currunt + tau(k);
    t(k+1) = t_currunt;
    [Energy(k+1), Masserror(k+1)] = EnergyMass(Spara,a,Mass0,u(:,:,3));
    if s <=length(snaptime)
        if t_currunt >= snaptime(s) && t_currunt <= snaptime(s) + taumax
            usnap(:,:,s)=u(:,:,3);
            s = s+1;
        end
    end
    % X = [tau(k) t_currunt Energy(k+1) tol];
    X = [t_currunt];
    k = k+1;
    disp(X);
end
toc;
end