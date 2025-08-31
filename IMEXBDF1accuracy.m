function u = IMEXBDF1accuracy(Spara,t,tau,rho,a,f,u0,S)
% IMEX BDF1 scheme for Square phase field crystal model;
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
epsi = 1-a;
taum = max(tau);

%% ------ initial ------
u = zeros(M+1,M+1,N+1); 
u(1:M,1:M,1) = u0(X,Y); 
u(M+1,:,1) = u(1,:,1); u(:,M+1,1) = u(:,1,1);

%% ------ solve scheme ------
Tol = 1e-12;
for k = 1:N
    tol = 1;
    A = 1/tau(k) - a*DLap - 2*DDLap - DDDLap;
    force = f(X,Y,t(k+1));
    F1 = 1/tau(k)*fft2(u(1:M,1:M,k)) + fft2(force(1:M,1:M));
    
    gradient_x = ifft2(Dx.*fft2(u(1:M,1:M,k)));
    gradient_y = ifft2(Dy.*fft2(u(1:M,1:M,k)));
    gradient_norm = gradient_x.^2 + gradient_y.^2;
    norm_times_gradient_x = gradient_norm.*gradient_x;
    norm_times_gradient_y = gradient_norm.*gradient_y;
    gradient = Dx.*fft2(norm_times_gradient_x(1:M,1:M)) + ...
        Dy.*fft2(norm_times_gradient_y(1:M,1:M));

    Fnl = -DLap.*gradient;
    F = F1 + Fnl;

    % u_temp1 = u(1:M,1:M,k);
    % while tol > Tol    
    %     u_temp2 = ifft2(F./A);
    %     tol = norm(u_temp1-u_temp2,"inf")
    %     u_temp1 = u_temp2;
    % end
    u(1:M,1:M,k+1) = ifft2(F./A);
    u(M+1,:,k+1) = u(1,:,k+1);
    u(:,M+1,k+1) = u(:,1,k+1); 
end

end