function [E,Masserror] = EnergyMass(Spara,a,Mass0,u)

M = Spara.M;
h = Spara.h;
DLap = Spara.DLap;
Dx = Spara.Dx;
Dy = Spara.Dy;
Dz = Spara.Dz;

ui = u(1:M,1:M,1:M);
fftui = fftn(ui);
gradient_norm4 = ((ifftn(Dx.*fftui)).^2 + (ifftn(Dy.*fftui)).^2 +...
    (ifftn(Dz.*fftui)).^2).^2;
E = 1/4*h^3*sum(sum(sum(gradient_norm4.^2))) +...
    a/2*h^3*sum(sum(sum(ui.^2))) +...
    1/2*h^3*sum(sum(sum(  2*ui.*ifftn(DLap.*fftui) + ...
    ifftn(DLap.*fftui).*ifftn(DLap.*fftui) )));
Masserror = h^3*sum(sum(sum(ui))) - Mass0;



end