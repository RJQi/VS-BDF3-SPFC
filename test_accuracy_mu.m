% spacemesh.m: generate uniform spatial meshes
% f_uexactf.m: define forcing term, exact solution, initial value
% timemesh.m: generate nonuniform/uniform temproal step sizes

% clear;
clearvars -except Uexact1 Uexact2;
L = 2*pi;
M = 128;
T = 1;
N0 = 16;

K = 5;
S = 3;

Spara = spacemesh(L,M);
% index1 = input('时间步长: 随机步长-1, 周期步长-2, 增强步长-3, 均匀步长-4:');
index2 = input('精确解: 1, 无精确解: 2');
[f,uexactf,u0] = f_uexactf(index2);


% reference solution with tau = 1/10000

[t,tau,rho,~] = timemesh(T,10000,4);
a1 = 0.55;
a2 = 0.85;
% Uexact1 = BDF3accuracy(Spara,t,tau,rho,a1,f,u0);
% Uexact2 = BDF3accuracy(Spara,t,tau,rho,a2,f,u0);

for i = 1:K
    N = N0*2^(i-1);

    % mu = 2;
    mu = 2;
    [t_per,tau_per,rho_per,maxtau_per(i)] = timemesh(T,N,2,mu);
    a = a1;
    uCSBDF3 = CSBDF3accuracy(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    uCSBDF32 = CSBDF3accuracy2(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    error(1,i) = Error(Spara.h,uCSBDF3,Uexact1,index2);
    error(2,i) = Error(Spara.h,uCSBDF32,Uexact1,index2);

    a = a2;
    uCSBDF3 = CSBDF3accuracy(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    uCSBDF32 = CSBDF3accuracy2(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    error(3,i) = Error(Spara.h,uCSBDF3,Uexact2,index2);
    error(4,i) = Error(Spara.h,uCSBDF32,Uexact2,index2);

    % mu = 3
    mu = 3;
    [t_per,tau_per,rho_per,maxtau_per(i)] = timemesh(T,N,2,mu);
    a = a1;
    uCSBDF3 = CSBDF3accuracy(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    uCSBDF32 = CSBDF3accuracy2(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    error(5,i) = Error(Spara.h,uCSBDF3,Uexact1,index2);
    error(6,i) = Error(Spara.h,uCSBDF32,Uexact1,index2);

    a = a2;
    uCSBDF3 = CSBDF3accuracy(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    uCSBDF32 = CSBDF3accuracy2(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    error(7,i) = Error(Spara.h,uCSBDF3,Uexact2,index2);
    error(8,i) = Error(Spara.h,uCSBDF32,Uexact2,index2);


end

N = N0*2.^((1:K)-1);

bdf3rate(1,:) = log(error(1,1:end-1)./error(1,2:end))./log(maxtau_per(1:end-1)./maxtau_per(2:end));
bdf3rate(2,:) = log(error(2,1:end-1)./error(2,2:end))./log(maxtau_per(1:end-1)./maxtau_per(2:end));
bdf3rate(3,:) = log(error(3,1:end-1)./error(3,2:end))./log(maxtau_per(1:end-1)./maxtau_per(2:end));
bdf3rate(4,:) = log(error(4,1:end-1)./error(4,2:end))./log(maxtau_per(1:end-1)./maxtau_per(2:end));