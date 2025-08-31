function [f,uexact,u0] = f_uexactf(index)
%define forcing term, exact solution and initial value

if index == 1
    uexact = @(x,y,t) cos(pi*x).*cos(pi*y).*exp(-t);
    u0 = @(x,y) uexact(x,y,0);
    f = @(x,y,t) exp(-3*t).*cos(pi*x).*cos(pi*y).*(exp(2*t)*(-1-2*kap*pi^2+4*kap*epsi^2*pi^4)...
        +6*kap*pi^2*cos(pi*x).^2.*cos(2*pi*y)-6*kap*pi^2*cos(pi*y).^2.*sin(pi*x).^2);

%     uexact = @(x,y,t) cos(x).*cos(y).*exp(-t);
%     u0 = @(x,y) uexact(x,y,0);
%     f = @(x,y,t) exp(-3*t).*cos(x).*cos(y).*(exp(2*t)*(-1-2*kap+4*kap*epsi^2)...
%         +6*kap*cos(x).^2.*cos(2*y)-6*kap*cos(y).^2.*sin(x).^2);
elseif index == 2
    uexact = [];
    u0 = @(x,y) cos(x).*cos(y); 
    % u0 = @(x,y) sin(x*pi*2).*cos(y*pi*2)/2/pi; 
    f = @(x,y,t) 0*x+0*y+0*t;
else 
    uexact = [];
    u0 = @(x,y) cos(x).*cos(y); 
    f = @(x,y,t) 0*x+0*y+0*t;
end

end