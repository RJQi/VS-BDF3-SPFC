% spacemesh.m: generate uniform spatial meshes
% f_uexactf.m: define forcing term, exact solution, initial value
% timemesh.m: generate nonuniform/uniform temproal step sizes

clear;
%clearvars -except Uexacta1 Uexacta2;
L = 2*pi;
M = 128;
T = 1.5;
N0 = 12;

K = 5;
S = 1;
mu = 1.7;

Spara = spacemesh(L,M);
%index1 = input('时间步长: 随机步长-1, 周期步长-2, 增强步长-3, 均匀步长-4:');
index2 = input('精确解: 1, 无精确解: 2');
[f,uexactf,u0] = f_uexactf(index2);

a1 = 0.55;
a2 = 0.85;
% reference solution with tau = 1/5120
[t,tau,rho,~] = timemesh(T,10000,4);
a = a1;
% Uexacta1 = BDF3accuracy(Spara,t,tau,rho,a,f,u0);
a = a2;
% Uexacta2 = BDF3accuracy(Spara,t,tau,rho,a,f,u0);


for i = 1:K
    N = N0*2^(i-1);
    % random stepsizes
    [t_random,tau_random,rho_random,maxtau_random(i)] = timemesh(T,N,1,mu);
    max_rho(i) = max(rho_random);
    min_rho(i) = min(rho_random(2:end));
    num_great(i) = sum(rho_random>1.7319);
    num_less(i) = sum(rho_random<0.5);
    a = a1;
    uCSBDF3 = CSBDF3accuracy(Spara,t_random,tau_random,rho_random,a,f,u0,S);
    uCSBDF32 = CSBDF3accuracy2(Spara,t_random,tau_random,rho_random,a,f,u0,S);
    error(1,i) = Error(Spara.h,uCSBDF3,Uexacta1,index2);
    error(2,i) = Error(Spara.h,uCSBDF32,Uexacta1,index2);

    a = a2;
    uCSBDF3 = CSBDF3accuracy(Spara,t_random,tau_random,rho_random,a,f,u0,S);
    uCSBDF32 = CSBDF3accuracy2(Spara,t_random,tau_random,rho_random,a,f,u0,S);
    error(3,i) = Error(Spara.h,uCSBDF3,Uexacta2,index2);
    error(4,i) = Error(Spara.h,uCSBDF32,Uexacta2,index2);

    % perodic stepsizes
    [t_per,tau_per,rho_per,maxtau_per(i)] = timemesh(T,N,2,mu);
    a = a1;
    uCSBDF3 = CSBDF3accuracy(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    uCSBDF32 = CSBDF3accuracy2(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    error(5,i) = Error(Spara.h,uCSBDF3,Uexacta1,index2);
    error(6,i) = Error(Spara.h,uCSBDF32,Uexacta1,index2);

    a = a2;
    uCSBDF3 = CSBDF3accuracy(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    uCSBDF32 = CSBDF3accuracy2(Spara,t_per,tau_per,rho_per,a,f,u0,S);
    error(7,i) = Error(Spara.h,uCSBDF3,Uexacta2,index2);
    error(8,i) = Error(Spara.h,uCSBDF32,Uexacta2,index2);

end

N = N0*2.^((1:K)-1);

bdf3rate(1,:) = log(error(1,1:end-1)./error(1,2:end))./log(N(1:end-1)./N(2:end));
bdf3rate(2,:) = log(error(2,1:end-1)./error(2,2:end))./log(N(1:end-1)./N(2:end));
bdf3rate(3,:) = log(error(3,1:end-1)./error(3,2:end))./log(N(1:end-1)./N(2:end));
bdf3rate(4,:) = log(error(4,1:end-1)./error(4,2:end))./log(N(1:end-1)./N(2:end));