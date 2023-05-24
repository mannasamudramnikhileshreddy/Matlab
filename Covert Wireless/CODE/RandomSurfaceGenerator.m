
%% References: 

%%          Dolcetti, G., and Krynkin, A. (2020), Matlab Codes for Two-Dimensional
%%      Scattering Surface Reconstruction using Broadband Data, Zenodo

%%          Dolcetti, G., Alkmim, M., Cuenca, J., De Ryck, L., Krynkin, A. (2020), 
%%      Robust Surface Shape Inversion with a Linear Microphones Array, submitted to 
%%      Journal of Sound and Vibrations

%% Author: Giulio Dolcetti g.dolcetti@sheffield.ac.uk
%% Date: 17-07-2020
%% Distribution: Creative Commons Attribution 4.0 International

%%-------------------------------------------------------------------------


%% Random Surface Generator:
%% creates random realisations of a surface with power function spectrum 
%% S(k) = S_0 (k/\kappa_0)^{-\alpha_0} for \kappa_0 < k < \kappa_1, and 
%% S(k) = S_0                          for k < \kappa_0

%% INPUT:
%% sigma_0:         surface standard deviation (m)
%% kappa_0:         saturation wavenumber (rad/m)
%% kappa_1:         wavenumber cut-off (rad/m)
%% alpha_0:         power spectrum slope (-)
%% Dx:              spatial grid size (m)
%% N:               number of grid points (-)

%% OUTPUT
%% x:               surface x-co-ordinate (m)   (N x 1)
%% z:               surface z-co-ordinate (m)   (N x 1)


function [x,z,Cs_AnalysE] = RandomSurfaceGenerator(sigma_0,kappa_0,kappa_1,alpha_0,Dx,N)
phfi=pi/18;
R=10*10^2;
Gtx=2/(1-cos(phfi/2));
Grx=2/(1-cos(phfi/2));
ks = 2*pi/Dx;
Dk = ks/N;
k = (0:N-1)'*Dk;
H=1;
n=1000;
alpha=4;
d1=3*(n^(1/2*alpha));
d=5*10^2;n=1000;
K=0.01;
f=500*10^8;
% d1=3*(n^(1/2*alpha));
rb=0.1*10^2;

% amplitude spectrum
Z = (k/kappa_0).^(-alpha_0/2);
Z(k<kappa_0) = 1;
Z(k>kappa_1) = 0;
Z(k==0) = 0;
% spectrum normalisation
Z = Z/sqrt(trapz(k,abs(Z).^2)/(2*pi));
 Cs_AnalysE=Grx/Gtx;
% random coefficients
RA = randn(size(Z));
RP = rand(size(Z))*2*pi;
A=H*Gtx*Grx;
Prx=A*d^-2*exp(-K*d);
Ptx=12;
c=3*10^8;
h=6.63 * 10^-34;
kb=1.38064852*10^-23;
T=213;
H=Ptx*(c^2/16*pi^2*f^2);
S=(h*f)/(exp(h*f/kb*T)-1);
GW=3.66;
GB=0.16;
r=10;
r1=12;
for lambda=0:0.01:1
Ib=A*lambda*phfi*exp(lambda*(rb)^2)*[expint(R*(K+lambda*rb))-expint(-rb*(K+lambda*rb))];
end
for lambda=0:0.01:0.1
Iw=A*lambda*phfi*exp(lambda*(rb)^2)*[expint(R*(K+lambda*rb))-expint(-rb*(K+lambda*rb))];
end
SINRB=(A*(d^-2)*exp(-K*d)*GB)/(S+Ib);
SINRW=(A*(d1^-2)*exp(-K*d1)*GW)/(S+Iw);

for i=1:11
  Cs_Analys(i)=log(1+SINRB)-log(1+SINRW)/log(1+SINRB); 
 end

x = (0:N-1)'*Dx; x = x-mean(x);
z = sigma_0*real((ifft(Z.*RA.*exp(1i*RP)*N/2,[],1,'symmetric')));







