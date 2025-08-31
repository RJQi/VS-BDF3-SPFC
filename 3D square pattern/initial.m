function u0 = initial(M,index)
%define initial value


if index == 1
    u0(1:M+1,1:M+1) = 0.5;
    u0(M/2,M/2) = 10;
    u0(M/4,M/2) = 10;
%     u0(1:M,1:M) = 0.1*rand(M,M)-0.05;
%     u0(M/2,M/2) = 10;
%     u0(M+1,:) = u0(1,:); u0(:,M+1) = u0(:,1);
elseif index == 2
    u0 = 0.285+0.4*rand(M+1,M+1,M+1)-0.2;  %generate nellsite
elseif index == 3
    load('u0ref.mat');
    % last 20:52
    u0 = 0.285*ones(M+1,M+1,M+1);
    u0(18:50,18:50,64:96) = u0ref;
    u0(64:96,64:96,18:50) = u0ref;
else
    u0 = 2*rand(M+1,M+1)-1;  
end

end