function u = BDF3accuracy(Spara,t,tau,rho,a,f,u0)
% Variable-step BDF3 scheme for Square phase field crystal model;
% Third order scheme for first two-level solution;
% Fourier spectral method in space

%% ------ mesh setting ------
 
N = length(tau);
M = Spara.M;
X = Spara.X;
Y = Spara.Y;
Dx = Spara.Dx;
Dy = Spara.Dy;
DLap = Spara.DLap;
DDLap = Spara.DDLap;
DDDLap = Spara.DDDLap;

%% ------ initial ------
a0f = @(x,y) (2*x+1)./(x+1)+x.*y./(x.*y+y+1);
a1f = @(x,y) -x./(x+1)-x.*y./(x.*y+y+1)-x.*y.^2.*(x+1)./(x.*y+y+1)./(y+1);
a2f = @(x,y) x.*y.^2.*(x+1)./(x.*y+y+1)./(y+1);
u = zeros(M+1,M+1,N+1); 
u(1:M,1:M,1) = u0(X,Y); u(M+1,:,1) = u(1,:,1); u(:,M+1,1) = u(:,1,1);
% u(1:M+1,1:M+1,1) = 0.5;
%     u(M/2,M/2,1) = 10;
%     u(M/4,M/2,1) = 10;
u(:,:,2) = SIDRK(Spara,tau(1),a,t(1),u(1:M,1:M,1),f);
u(:,:,3) = SIDRK(Spara,tau(2),a,t(2),u(1:M,1:M,2),f);

%% ------ solve scheme ------
Tol = 1e-12;
for k=3:N
    tol = 1;
    a0 = a0f(rho(k),rho(k-1));
    a1 = a1f(rho(k),rho(k-1));
    a2 = a2f(rho(k),rho(k-1));
    A = a0/tau(k) - a*DLap - 2*DDLap - DDDLap;
    force = f(X,Y,t(k+1));
    F1 = a0/tau(k)*fft2(u(1:M,1:M,k)) - a1/tau(k-1)*fft2(u(1:M,1:M,k)-u(1:M,1:M,k-1))...
        - a2/tau(k-2)*fft2(u(1:M,1:M,k-1)-u(1:M,1:M,k-2)) + fft2(force(1:M,1:M));

    u_temp1 = u(1:M,1:M,k);
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
        tol = norm(u_temp1-u_temp2,"inf");
        u_temp1 = u_temp2;
    end
    u(1:M,1:M,k+1) = u_temp1;
    u(M+1,:,k+1) = u(1,:,k+1);
    u(:,M+1,k+1) = u(:,1,k+1); 
end

end