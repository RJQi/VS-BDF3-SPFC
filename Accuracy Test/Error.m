function maxl2error = Error(h,u,Uexact,index)
%calculate maximum L2 norm error

N = size(u,3)-1;
% L2error = zeros(1,N+1);
if index == 1
    for k = 1:N+1
        L2error(k) = sqrt(h^2*sum(sum((Uexact(:,:,k)-u(:,:,k)).^2)));
    end
else
    L2error = sqrt(h^2*sum(sum((Uexact(:,:,end)-u(:,:,end)).^2)));
%     r = 5120/N;
%     for k = 1:N+1
%         L2error(k) = sqrt(h^2*sum(sum((Uexact(:,:,k+8)-u(:,:,k)).^2)));
%     end
end
maxl2error = max(L2error);
end