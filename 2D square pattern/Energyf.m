function [E] = Energyf(Spara,a,u)

M = Spara.M;
h = Spara.h;
DLap = Spara.DLap;
Dx = Spara.Dx;
Dy = Spara.Dy;
% epsi = 1-a;

% N = size(u,3);
% E = zeros(1,N);
% Masserror = zeros(1,N);
% Mass0 = h^2*sum(sum(u(1:M,1:M,1)));

ui = u(1:M,1:M);
gradient_norm4 = ((ifft2(Dx.*fft2(ui))).^2 + (ifft2(Dy.*fft2(ui))).^2).^2;
E = 1/4*h^2*sum(sum(gradient_norm4)) +...
    a/2*h^2*sum(sum(ui.^2)) +...
    1/2*h^2*sum(sum(  2*ui.*ifft2(DLap.*fft2(ui)) + ...
    ifft2(DLap.*fft2(ui)).*ifft2(DLap.*fft2(ui)) ));
% Masserror = h^2*sum(sum(ui)) - Mass0;
end

