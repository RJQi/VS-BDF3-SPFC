function  Snapshot(X,Y,t,u)
% function Snapshot(X,Y,t,un_Snap)
%% ====显示不同时刻数值解的俯视图===
n=length(t);
time = [1 2 5 10];
time_length = length(time);
for k=1:n
    s = t(k);
    for j = 1:time_length
        if s >= time(j) && s <= time(j) + 0.1
            tn(j)=k;
        end
    end
end
% tn(time_length+1) = n;
%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:time_length
    figure;
    un_Snap(:,:,k)=u(:,:,tn(k));
    colormap jet;
    pcolor(X,Y,un_Snap(:,:,k));
    shading interp
    set(gca,'position',[0 0 1 1]); % 这个命令是去掉多余的留白
    axis off
%     zf = num2str(k);
%     print(zf,'-dpng');
end
%%%%%%%%%%%%%%%%%%%%%%%%%