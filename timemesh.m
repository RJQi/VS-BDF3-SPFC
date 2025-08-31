function [t,tau,rho,maxtau] = timemesh(T,N,index,mu)
% time step size generate
% tau = zeros(1,N+1);
if index == 1
    rn = rand(1,N);
    tau = T*rn/sum(rn);
elseif index == 2
    % mu = 1.7;
    % tau(1) = 2/(1+mu)/N;
    % for i = 2:N
    %      if mod(i,2) == 0
    %          tau(i) = mu*tau(1);
    %      else
    %          tau(i) = tau(1);
    %      end
    % end
    % mu = 1.7;
    tau(1) = 4*T/(1+2*mu+mu^2)/N;
    tau(2) = mu*tau(1);
    tau(3) = mu^2*tau(1);
    tau(4) = mu*tau(1);
    for i = 0:(N-8)/4
        tau(4*i+5:4*i+8) = tau(1:4);
    end
    % for i = 2:N
    %     if mod(i,2) == 0
    %         tau(i) = mu*tau(1);
    %     elseif mod(i,3) == 0
    %         tau(i) = mu^2*tau(1);
    %     else
    %         tau(i) = tau(1);
    %     end
    % end
elseif index == 3
    
    if N == 40
        load tauref.mat;
        tau = tauref;
    elseif N == 80
        load tau80.mat;
        tau = tau80;
%         tau(1:N/2) = tauref(1:N/2)/2;
%         tau(N/2+1:N) = tau(1:N/2);
    elseif N == 160
        load tau160.mat;
        tau = tau160;
    elseif N == 320
        load tau320.mat;
        tau = tau320;
    else
        load tau640.mat;
        tau = tau640;
    end

%     mu = 1.2;
%     tau(1) = T*(1-mu)/(1-mu^(N/2))/(1+1/mu);
%     for i = 2:N/2
%         tau(i) = mu*tau(i-1);
%     end
%     for i = N/2+1:N
%         tau(i) = tau(i-1)/mu;
%     end
else
    tau = ones(1,N)*T/N;
end
rho(1) = 0;
rho(2:N) = tau(2:N)./tau(1:N-1);
t(1) = 0;
t(2:N+1) = cumsum(tau);
maxtau = max(tau);
end