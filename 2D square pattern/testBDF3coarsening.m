clear;
% format long;
L = 400;
M = 256;
T = 5000;
a = 0.95;


rhomax = 1.73;
rhomin = 0.5;
Sp = 1;

Spara = spacemesh(L,M);
index = input('初值条件: 初值1, 初值2, 初值3, 初值4:');
% savename = input('The name of save:');

u0 = initial(M,index);
% load('initialvalue0226.mat')
% snaptime = [1 20 40 60 80 100];
% snaptime = [1 5 10 50 100 200 300 400 500 600 700 800 900 1000];
% snaptime = [1 2 5 10 30 50 100 200 400 800 1000 1200 1600 2000];
snaptime = [1 2 5 10 20 30 40 50 100 200 300 400 1000 1500 2000 3000 4000 5000];
% snaptime = [1 2 5 10 30 50 100 200 400 700 1000 1500 2000 2500 3000 3500 ...
%     4000 4500 5000 6000 7000 8000 9000 10000 12000 15000];
beta = 100;

% tic;
% taumin = 1e-2;
% taumax = 5e-1;
% [t1,tau1,Energy1,Masserror1,usnap1] = CSBDF3coarsening(Spara,T,taumin,taumax,rhomax,rhomin,beta,a,u0,snaptime,Sp,0);
% totaltime1 = toc;

tic;
taumin = 1e-2/2;
taumax = 1e-1/2;
[t4,tau4,Energy4,Masserror4,usnap4] = CSBDF3coarsening(Spara,T,taumin,taumax,rhomax,rhomin,beta,a,u0,snaptime,Sp,0);
totaltime4 = toc;

% tic;
% taumin = 0.01;
% taumax = 0.01;
% [t3,tau3,Energy3,Masserror3,usnap3] = CSBDF3coarsening(Spara,T,taumin,taumax,rhomax,rhomin,beta,a,u0,snaptime,Sp,1);
% totaltime3 = toc;
% [E1,Masserror1] = EnergyMass(Spara,a,u1);


% save('new0225.mat');
% % figure;
% % plot(t(1:end-1),Energy(1:end-1));
% % figure;
% % plot(t(1:end-1),Masserror(1:end-1));
% % figure;
% % plot(t(1:end-1),tau);
% % 

% snaptime = [1 8 9 11 13 19];
x = 0:Spara.h:L; y = x;
[X,Y] = meshgrid(x,y);
for k = 1:length(snaptime)
    figure;
%     un_Snap(:,:,k)=u(:,:,tn(k));
    colormap hsv;
    pcolor(X,Y,usnap4(:,:,k));
    shading interp
    set(gca,'position',[0 0 1 1]); % 这个命令是去掉多余的留白
    axis off
%     zf = num2str(k);
%     print(zf,'-dpng');
end

% mycolor=[0 0 0;1 1 1];
x = 0:Spara.h:L; y = x;
[X,Y] = meshgrid(x,y);
figure
colormap(hsv);
pcolor(X,Y,usnap1(:,:,end));
shading interp
set(gca,'position',[0 0 1 1]); % 这个命令是去掉多余的留白
    axis off
figure
colormap(hsv);
pcolor(X,Y,usnap2(:,:,end));
shading interp
set(gca,'position',[0 0 1 1]); % 这个命令是去掉多余的留白
    axis off
figure
colormap(hsv);
pcolor(X,Y,usnap3(:,:,end));
shading interp
set(gca,'position',[0 0 1 1]); % 这个命令是去掉多余的留白
    axis off