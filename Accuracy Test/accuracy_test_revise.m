% spacemesh.m: generate uniform spatial meshes
% f_uexactf.m: define forcing term, exact solution, initial value
% timemesh.m: generate nonuniform/uniform temproal step sizes

clear;
% clearvars -except Uexact;
L = 2*pi;
M = 128;
T = 1;
N0 = 20;
N1 = [150 1000 8000];
N2 = [20 40 80];

K = 3;
S = 1;
mu = 1.7;

Spara = spacemesh(L,M);
index1 = input('时间步长: 随机步长-1, 周期步长-2, 增强步长-3, 均匀步长-4:');
index2 = input('光滑解: 2, 不光滑解: 3');
[f,uexactf,u0] = f_uexactf(index2);

a = 0.1;
% reference solution with tau = 1/5120
[t,tau,rho,~] = timemesh(T,10000,4);
Uexact = BDF3accuracy(Spara,t,tau,rho,a,f,u0);


for i = 1:K
    tic;
    % N01 = 150;
    N = N1(i);
    [t_random,tau_random,rho_random,maxtau_random(i)] = timemesh(T,N,index1,mu);
    uIMEXBDF1 = IMEXBDF1accuracy(Spara,t_random,tau_random,rho_random,a,f,u0,S);
    error(1,i) = Error(Spara.h,uIMEXBDF1,Uexact,index2);
    toc
end

for i = 1:K
    % N = N0*2^(i-1);
    tic;
    N = N2(i);
    [t_random,tau_random,rho_random,maxtau_random(i)] = timemesh(T,N,index1,mu);
    [uCSBDF3,num_iter1] = CSBDF3accuracy(Spara,t_random,tau_random,rho_random,a,f,u0,S);
    error(2,i) = Error(Spara.h,uCSBDF3,Uexact,index2);
    sum(num_iter1);
    toc
end

% tic;
% for i = 1:K
%     % N = N0*2^(i-1);
%     N = N2(i);
%     [t_random,tau_random,rho_random,maxtau_random(i)] = timemesh(T,N,index1,mu);
%     [uCSBDF32,num_iter2] = CSBDF3accuracy2(Spara,t_random,tau_random,rho_random,a,f,u0,S);   
%     error(3,i) = Error(Spara.h,uCSBDF32,Uexact,index2);
%     sum(num_iter2)
% end
% t3 = toc;


N = N0*2.^((1:K)-1);

rate(1,:) = -log(error(1,1:end-1)./error(1,2:end))./log(N1(1:end-1)./N1(2:end));
rate(2,:) = -log(error(2,1:end-1)./error(2,2:end))./log(N2(1:end-1)./N2(2:end));
% rate(3,:) = -log(error(3,1:end-1)./error(3,2:end))./log(N2(1:end-1)./N2(2:end));
