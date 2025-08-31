clear;
format long;
L = 128;
M = 96;
T = 1000;
a = 0.75;

taumin = 5e-3;
taumax = 1;
rhomax = 1.73;
rhomin = 0.5;

Spara = spacemesh(L,M);
index = input('初值条件: 初值1, 初值2, 初值3, 初值4:');
u0 = initial(M,index);

snaptime = [1 5 10 50 100 200 300 400 500 600 700 800 900 1000];
% snaptime = [0.1 0.2 0.5 0.7 1.0];
beta = 10;
tic;
[t,tau,Energy,Masserror,usnap] = BDF3coarsening(Spara,T,taumin,taumax,rhomax,rhomin,beta,a,u0,snaptime);
totaltime = toc;

save('data3D-1.mat');

