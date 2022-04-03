function [U1,U2] = clle(U1,U2,linop1,linop2,nlop1,nlop2,MAXTAU) %Coupled LLE
h = 1e-3; %Step size
ctau = 0; %Current tau. Runs up to MAXTAU


% linop1 = @(h,uf1) h*(0.5i*sign(beta2)*kaxis.^2 + 1i*(1i-alpha)).*exp(uf1) + 1i*h*sqrt(P1);
% nlop1 = @(h,u1,u2) h*1i*(gs*abs(u1).^2 + gc*abs(u2).^2).*exp(u1);
%
% linop2 = @(h,uf2) h*(0.5i*sign(beta2)*kaxis.^2 + 1i*(1i-alpha)).*exp(uf2) + 1i*h*sqrt(P2);
% nlop2 = @(h,u1,u2) h*1i*(gc*abs(u1).^2 + gs*abs(u2).^2).*exp(u2);

while ctau<MAXTAU
    U1ori = U1;
    U2ori = U2;
    U1F = fft(nlop1(h/2,U1,U2));
    U1 = ifft(linop1(h,U1F));
    U1 = nlop1(h/2,U1,U2);
    
    U2F = fft(nlop2(h/2,U1,U2));
    U2 = ifft(linop2(h,U2F));
    U2 = nlop2(h/2,U1,U2);

    ctau = ctau + h;
end

end