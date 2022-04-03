%Code to solve coupled LLEs based on PRA
%Author: Pradyoth Shandilya
%Contact: shandilya@umbc.edu

close all
clearvars
clc
MAXTAU = 150; %Total propagation "distance" of the equations


NT = 2^10; %Discretization
beta2 = -0.01; %beta2(bar) 
alpha = 2;
P1 = 0.5;
P2 = 1.5; %alpha - P1
gc = 1;
gs = 1;

K = 8*sqrt(abs(beta2));

SIGMA = pi/sqrt(abs(beta2)); %Max value of sigma (see fig 1). Equivalent to T in NLSE

dsig = 2*SIGMA/NT;
sigma = (-NT/2:1:NT/2-1).'*dsig; %The sigma axis (scaled azimuthal coordinate)

dk = 1*pi/SIGMA;
kaxis = dk * [(0:NT/2-1) (-NT/2:-1)].'; %Wavenumber (scaled) grid

%Initial conditions are defined based on equations 8, 10 and 11 of the [1]

U0_1 = 1i*sqrt(P1)*ones(NT,1); %Eq 8 - CW(1)
U0_2 = 1i*sqrt(P2)*ones(NT,1); %Eq 8 - CW(2)

%The following is from eq 11. 
omega_mp = -1i*(1 - sqrt(-K^4/4 + alpha*K^2/2*(-sign(beta2)+1))); %Eq 9
Delta = 2*(1 - 1i*omega_mp);

C = 0.1;

x1 = C;
y1 = C*sign(beta2)*K^2 / Delta; 
x2 = -(C*sign(beta2)*K^2*sqrt(P1*P2))/(Delta^2/4 + K^4/4 + sign(beta2)*K^2*P2);
y2 = -(C*K^4*sqrt(P1*P2)/Delta)/(Delta^2/4 + K^4/4 + sign(beta2)*K^2*P2);

U1 = U0_1 + real(x1*exp(1i*K*sigma)) + 1i*real(y1*exp(1i*K*sigma)); 
U2 = U0_2 + real(x2*exp(1i*K*sigma)) + 1i*real(y2*exp(1i*K*sigma)); 


%---------------------------------------------------------------------------------------------------

%Solving the coupled LLEs (Eq 4)

%Define operators:
linop1 = @(h,uf1) exp(h*(0.5i*sign(beta2)*kaxis.^2 + 1i*(1i-alpha))).*(uf1 + fft(U0_1)./(0.5i*sign(beta2)*kaxis.^2 + 1i*(1i-alpha))) - fft(U0_1)./(0.5i*sign(beta2)*kaxis.^2 + 1i*(1i-alpha));
nlop1 = @(h,u1,u2) exp(h*1i*(gs*abs(u1).^2 + gc*abs(u2).^2)).*(u1);

linop2 = @(h,uf2) exp(h*(0.5i*sign(beta2)*kaxis.^2 + 1i*(1i-alpha))).*(uf2 + fft(U0_2)./(0.5i*sign(beta2)*kaxis.^2 + 1i*(1i-alpha))) - fft(U0_2)./(0.5i*sign(beta2)*kaxis.^2 + 1i*(1i-alpha));
nlop2 = @(h,u1,u2) exp(h*1i*(gc*abs(u1).^2 + gs*abs(u2).^2)).*(u2);

[U1out,U2out] = clle(U1,U2,linop1,linop2,nlop1,nlop2,MAXTAU);

figure(1); plot(sigma,abs(U1out).^2,'-','color','red','LineWidth',2); 
hold on; 
plot(sigma,abs(U1).^2,'--','color','red','LineWidth',2); 
plot(sigma,abs(U2).^2,'--','color','blue','LineWidth',2); 
plot(sigma,abs(U2out).^2,'-','color','blue','LineWidth',2)
xlim([sigma(1),sigma(end)])
xlabel('\sigma','FontSize',16)
legend('|U^{(1)}(\tau=150)|^2','|U^{(1)}(\tau=0)|^2','|U^{(2)}(\tau=0)|^2','|U^{(2)}(\tau=150)|^2')
Spectrum1 = (fft(abs(U1out).^2));
Spectrum2 = real(fft(abs(U2out).^2));


Spectrum_nz1 = Spectrum1/max(Spectrum2);
modenum1 = linspace(-NT/2,NT/2,NT);
figure(2); plot(modenum1,fftshift(real(10*log10(Spectrum_nz1))),'color','red','LineWidth',2)
xlabel('Mode Number $m - \bar{m}$','Interpreter','Latex')
ylim([-80 0])
ylabel('F_t[|U^(1)|^2]')

Spectrum_nz2 = Spectrum2/max(Spectrum2);
Spectrum_dB2 = 10*log10(Spectrum_nz2/10);
modenum2 = linspace(-NT/2,NT/2,NT);
figure(3); plot(modenum2,fftshift(real(10*log10(Spectrum_nz2))),'color','blue','LineWidth',2)
xlabel('Mode Number $m - \bar{m}$','Interpreter','Latex')
ylim([-80 0])
ylabel('F_t[|U^(1)|^2]')