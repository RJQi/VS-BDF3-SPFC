% spacemesh.m: generate uniform spatial meshes
% f_uexactf.m: define forcing term, exact solution, initial value
% timemesh.m: generate nonuniform/uniform temproal step sizes

clear;
L = 2*pi;
M = 128;
T = 5;
N0 = 2;

K = 4;
% S = 1;

Spara = spacemesh(L,M);
index1 = input('时间步长: 随机步长-1, 周期步长-2, 增强步长-3, 均匀步长-4:');
index2 = input('精确解: 1, 无精确解: 2');
[f,uexactf,u0] = f_uexactf(index2);

a = 0.65;
% reference solution with tau = 1/10000
[t,tau,rho,~] = timemesh(T,10000,4);
Uexact = BDF3accuracy(Spara,t,tau,rho,a,f,u0);
A = [0.95 0.85 0.75 0.65];

mu = 1.7;
for i = 1:K
    N = N0*2^(i-1);
    % N = N0;
    % a = A(i);

    S = 35;
    [t_per,tau_per,rho_per,maxtau_per(i)] = timemesh(T,N,index1,mu);
    uCSBDF3 = CSBDF3accuracy(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    % uCSBDF32 = CSBDF3accuracy2(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    error(1,i) = Error(Spara.h,uCSBDF3,Uexact,index2);
    % error(2,i) = Error(Spara.h,uCSBDF32,Uexact,index2);

    S = 10;
    [t_per,tau_per,rho_per,maxtau_per(i)] = timemesh(T,N,index1,mu);
    uCSBDF3 = CSBDF3accuracy(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    % uCSBDF32 = CSBDF3accuracy2(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    error(2,i) = Error(Spara.h,uCSBDF3,Uexact,index2);
    % error(4,i) = Error(Spara.h,uCSBDF32,Uexact,index2);

    S = 5;
    [t_per,tau_per,rho_per,maxtau_per(i)] = timemesh(T,N,index1,mu);
    uCSBDF3 = CSBDF3accuracy(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    % uCSBDF32 = CSBDF3accuracy2(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    error(3,i) = Error(Spara.h,uCSBDF3,Uexact,index2);
    % error(6,i) = Error(Spara.h,uCSBDF32,Uexact,index2);

    S = 1;
    [t_per,tau_per,rho_per,maxtau_per(i)] = timemesh(T,N,index1,mu);
    uCSBDF3 = CSBDF3accuracy(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    % uCSBDF32 = CSBDF3accuracy2(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    error(4,i) = Error(Spara.h,uCSBDF3,Uexact,index2);
    % error(8,i) = Error(Spara.h,uCSBDF32,Uexact,index2);

    S = 0;
    [t_per,tau_per,rho_per,maxtau_per(i)] = timemesh(T,N,index1,mu);
    uCSBDF3 = CSBDF3accuracy(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    % uCSBDF32 = CSBDF3accuracy2(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    error(5,i) = Error(Spara.h,uCSBDF3,Uexact,index2);
    % error(10,i) = Error(Spara.h,uCSBDF32,Uexact,index2);

end

N = N0*2.^((1:K)-1);

bdf3rate(1,:) = log(error(1,1:end-1)./error(1,2:end))./log(N(1:end-1)./N(2:end));
bdf3rate(2,:) = log(error(2,1:end-1)./error(2,2:end))./log(N(1:end-1)./N(2:end));
bdf3rate(3,:) = log(error(3,1:end-1)./error(3,2:end))./log(N(1:end-1)./N(2:end));
bdf3rate(4,:) = log(error(4,1:end-1)./error(4,2:end))./log(N(1:end-1)./N(2:end));