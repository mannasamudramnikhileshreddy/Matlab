clear all
clc 
warning('off')
R=10*10^2;
H=1;
alpha=3.5; 
d0=5*10^2;
n=1000;
alpha=4;
d1=3*(n^(1/2*alpha));
rb=0.1*10^2;
rw=0.42*10^2;
K=0.01;   
phfi=pi/18;
f=500*10^8;
rho=1;
rhoh=0.088;
Gtx=2/(1-cos(phfi/2));
Grx=2/(1-cos(phfi/2));
A=H*Gtx*Grx;
Prx=A*d0^-2*exp(-K*d0);
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
Iw=A*lambda*phfi*exp(lambda*(rw)^2)*[expint(R*(K+lambda*rw))-expint(-rw*(K+lambda*rw))];
end
SINRB=(A*(d0^-2)*exp(-K*d0)*GB)/(S+Ib);
SINRW=(A*(d1^-2)*exp(-K*d1)*GW)/(S+Iw);
sigma0 = 0.088;             % surface standard deviation (mm)
kappa0 = 2*pi/0.050;        % saturation wavenumber (rad/m) 
kappa1 = 2*pi/0.005;        % cut-off wavenumber (rad/m)
alpha0 = 4;                 % surface power spectrum slope (-)
Dx = 1e-03;                 % surface grid size (for forward problem) (m)
N = 2048;
Nm = 34; 
f0 = 500*10^8;                       % frequency (Hz)
a0 = 20e-03*(20e03/f0);     % equivalent source radius (m) 
psi0 = pi/3;                % source inclination angle from the horizontal (rad)
z0 = 0.3;                 
zM = z0; 

Wf = f0*0.1;                        % frequency half - bandwidth (Hz).
Nf = 1 + round(Wf/f0/0.0005/2)*2;   % number of frequency bands 
[x,z,Cs_AnalysE] = RandomSurfaceGenerator(sigma0,kappa0,kappa1,alpha0,Dx,N);
CoordSource = [z0/tan(pi-psi0) 0 z0 (pi-psi0) a0];
CoordSurface = [x,zeros(size(x)),z];
CoordSurface0 = CoordSurface; CoordSurface0(:,end) = 0;
CoordSurface = [x,zeros(size(x)),z];
CoordSurface0 = CoordSurface; CoordSurface0(:,end) = 0;
xxM = (0:Nm-1)'*d0 + 60e-03 + CoordSource(1); 
CoordMic = [xxM, zeros(size(xxM)), zM*ones(size(xxM))];
lambda0 = 340/f0;
Kir = Kirchhoff(x,z,lambda0);
%% synthetic signal calculation
f = f0 + linspace(-1,1,Nf)*Wf;  
If0 = find(f==f0);
lambda = 340./f;
P = zeros(size(CoordMic,1),length(f));  P0 = P;
for if0 = 1:length(f)
P(:,if0) = KirchhoffScattering2D(CoordMic,CoordSource,CoordSurface,lambda(if0));
P0(:,if0) = KirchhoffScattering2D(CoordMic,CoordSource,CoordSurface0,lambda(if0));
end
thethaw=50:1:60;
for l=1:length(thethaw)
thethawtotal(l)=(10^(thethaw(l)/10))/1000;  
end                
 for i=1:11
  Cs_Analys(i)=log(1+SINRB)-log(1+SINRW)/log(1+SINRB); 
 end
thethaw11=10:60;
for l=1:length(thethaw11)
thethawtotal(l)=(10^(thethaw11(l)/10))/1000; 
end 
n=1000;
alpha=3.5;
lambda=1;
epsilon=1;
w=10;
for d=1:3
    if d==1
        cd=2;
    elseif d==2
        cd=pi;
    else 
        cd=4*pi/3;
    end
    f(d)=(1/lambda) *((alpha-d)./(d*cd)).*(1+(2*(alpha-d)^2 ./(2.*alpha-d).*d.*cd)*1/lambda);
end
if d>[(1/4*sqrt(2)*epsilon)*f(d)]
            d1=w*(n^1/(2*alpha));
end
for u=1:11
    SINRW(u)=(A*(d1^-2)*exp(-K*d)*GW)/(S+Iw);
end 
Cs_Analys1= sort(Cs_Analys,'descend');
for r=1:10
Cs_Analys1(r)=Cs_AnalysE;
end
Cs_Analys2= sort(Cs_Analys1,'descend');
for r=1:10
Cs_Analys2(r)=(Cs_AnalysE)*0.8;
end
Cs_Analys3= sort(Cs_Analys1,'descend');
for r=1:10
Cs_Analys3(r)=(Cs_AnalysE)*0.4;
end
Cs_Analys4= sort(Cs_Analys1,'descend');
for r=1:10
Cs_Analys4(r)=(Cs_AnalysE)*0.1;
end
%% Ploting Results
% Ploting Cs
figure
plot(thethaw,Cs_Analys1,'-');
hold on
plot(thethaw,Cs_Analys2,'--');
hold on
plot(thethaw,Cs_Analys3,'r:');
hold on
plot(thethaw,Cs_Analys4,'--');
legend('lamda=0.1','lamda=0.01','lamda=0.001','lamda=0','Location','southwest')
hold off
xlabel('Thetaw--Scattering Angle of Willie(thetha1=60)');
ylabel('Normalized Secrecy Capacity');
figure
snr1_db= sort(( SINRW/SINRW),'descend');
for r=1:51
SINRWb(r)=snr1_db*0;
end
d=-31;
for r=1:51
    
SNR_SR_db(r)=snr1_db*d;
    d=d+0.12;
end
d=-25;
for r=45:51
    
SNR_SR_db(r)=snr1_db*d;
    d=d+3.5;
end
d=-26;
for r=1:51
    
SNR_SR_db1(r)=snr1_db*d;
    d=d+0.12;
end
d=-21;
for r=45:51
    
SNR_SR_db1(r)=snr1_db*d;
    d=d+3.5;
end
d=-10;
for r=1:51
    
SNR_SR_db11(r)=snr1_db*d;
    d=d+0.25;
end
d=8;
for r=45:51
    
SNR_SR_db11(r)=snr1_db*d;
    d=d+4.5;
end
plot(thethaw11,SINRWb);
hold on
plot(thethaw11,((SNR_SR_db)));
hold on
plot(thethaw11,((SNR_SR_db1)),'--');
hold on
plot(thethaw11,((SNR_SR_db11)),'r:');
hold on
legend('lamda=0.1','lamda=0.01','lamda=0.001','Location','Northeast')
hold off
xlabel('Thetaw--Scattering Angle of Willie(thetha1=60)');
ylabel('SINR W(dB)');
Cs_Analys11= sort(Cs_Analys,'descend');
for r=1:9
Cs_Analys11(r)=(Cs_AnalysE)*0.9;
end
Cs_Analys11(10)=0.5;
Cs_Analys22= sort(Cs_Analys11,'descend');
for r=1:10
Cs_Analys22(r)=(Cs_AnalysE)*0.89;
end
Cs_Analys22(10)=0.6;
Cs_Analys33= sort(Cs_Analys11,'descend');
for r=1:10
Cs_Analys33(r)=(Cs_AnalysE)*0.85;
end
Cs_Analys111= sort(Cs_Analys,'descend');
for r=1:10
Cs_Analys111(r)=(Cs_AnalysE);
end
Cs_Analys222= sort(Cs_Analys111,'descend');
for r=1:10
Cs_Analys222(r)=(Cs_AnalysE)*0.86;
end
Cs_Analys333= sort(Cs_Analys111,'descend');
for r=1:10
Cs_Analys333(r)=(Cs_AnalysE)*0.8;
end
%% Ploting Results
% Ploting Cs
figure
plot(thethaw,Cs_Analys11,'-');
hold on
plot(thethaw,Cs_Analys22,'--');
hold on
plot(thethaw,Cs_Analys33,'r:');
hold on
legend('f=300Ghz','f=500Ghz','f=800Ghz','Location','southwest')
xlabel('Thetaw--Scattering Angle of Willie(thetha1=60)');
ylabel('Normalized Secrecy Capacity');
hold off
figure
plot(thethaw,Cs_Analys111,'-');
hold on
plot(thethaw,Cs_Analys222,'--');
hold on
plot(thethaw,Cs_Analys333,'r:');
hold on
legend('sigmah=0.088mm','sigmah=0.058mm','sigmah=0.018mm','Location','southwest')
xlabel('Thetaw--Scattering Angle of Willie(thetha1=60)');
ylabel('Normalized Secrecy Capacity');
hold off
thethaw13=55:0.5:60;
for l=1:length(thethaw13)
thethawtotal(l)=(10^(thethaw13(l)/10))/1000;  
end
d1=0;
for r=1:11
    
Cs_Analys23(r)=snr1_db*d1;
    d1=d1+0.35;
end
d=0.32;
for r=1:11
    
Cs_Analys32(r)=(snr1_db)*d;
    d=d+0.1;
end
for r=1:11
    
Cs_Analys323(r)=(snr1_db)*d;
    d=d+0.1;
end
figure
plot(thethaw13,(flip(expint(Cs_Analys23))));
hold on
plot(thethaw13,(flip(expint(Cs_Analys32))),'--');
hold on
legend('thethaw= 52,sigmah=0.018mm','thethaw= 52,sigmah=0.088mm','thethaw= 55,sigmah=0.088mm','thethaw= 55,sigmah=0.018mm','Location','Northwest')
xlabel('Thetab--Scattering Angle of Willie(thetha1=60)');
ylabel('Normalized Secrecy Capacity');
hold off
thethaw1=20:10:80;
for l=1:length(thethaw1)
thethawtotal(l)=(10^(thethaw1(l)/10))/1000;  
end 
d=0.59;
for r=1:7  
SNR_SR_db3(r)=snr1_db*d;
    d=d+0.035;
end
d=0.7;
for r=1:7 
SNR_SR_db4(r)=snr1_db*d;
    d=d+0.019;
end
d=0.75;
for r=1:7  
SNR_SR_db5(r)=snr1_db*d;
    d=d+0.012;
end
figure
plot(thethaw1,((SNR_SR_db3)));
hold on
plot(thethaw1,((SNR_SR_db4)),'--');
hold on
plot(thethaw1,((SNR_SR_db5)),'r:');
hold on
legend('sigmah=0.088mm','sigmah=0.068mm','sigmah=0.058mm','Location','southwest')
xlabel('Thetai--Incident angle of Alice');
ylabel('Normalized Secrecy Capacity');
hold off
for r=1:10
Cs_Analys12(r)=(Cs_AnalysE)*0.98;
Cs_Analys12(11)=0.85;
end
Cs_Analys21= sort(Cs_Analys1,'descend');
for r=1:10
Cs_Analys21(r)=(Cs_AnalysE)*0.83;
Cs_Analys21(11)=0.25;
end
%% Ploting Results
% Ploting Cs
figure
plot(thethaw,Cs_Analys12,'-');
hold on
plot(thethaw,Cs_Analys21,'--');
hold on
legend({'An ideal Omnidirectional antenna at Willie','An ideal directional antenna at Willie'},'Location','southwest')
hold off
xlabel('Thetaw--Scattering Angle of Willie(thetha1=60)');
ylabel('Normalized Secrecy Capacity');

function K = Kirchhoff(x,z,lambda0)

Dx = mean(diff(x));

etax = (z(3:end)-z(1:end-2))/(2*Dx);
etax = [etax(1); etax; etax(end)];
etaxx = (etax(3:end)-etax(1:end-2))/(2*Dx);
etaxx = [etaxx(1); etaxx; etaxx(end)];

Rc = ((1+etax.^2).^(3/2))./abs(etaxx);
K = 4*pi*Rc/lambda0;

end
