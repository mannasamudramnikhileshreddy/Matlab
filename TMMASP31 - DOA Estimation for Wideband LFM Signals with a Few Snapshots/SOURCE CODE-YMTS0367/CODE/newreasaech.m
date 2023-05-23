tic
clc; 
clear; 
close all;
warning off;
T=1;
t=-T:0.01:T;
a1=5;%amplitude of a signal
fc=2.4*10^6;%centre frequency
mu=8*10^13;%chirp rate
X=a1*exp((1i*2*pi*fc*t)+(1i*pi*mu*(t).^2));%LFM signal emitted by tx
U=16;%range of frequency
%FrFT starts%%%%%
a=1;%transform order of FrFT
p1=0;q=0; 
 for t=-T:0.01:T
  p1=p1+1;
    for u=-U:0.1:U
        q=q+1;
        x1(p1,q)=(sqrt(1-1i*cot(a*pi/2)))*exp(1i*pi*((u^2)*cot(a*pi/2)+(t^2)*cot(a*pi/2)-2*u*t*(1/sin(a*pi/2))));%rotation angle with alpha=p*pi/2;
    end
    q=0;
 end 
u=-U:0.1:U;
x2=0.01*(X*x1);
t=-T:0.01:T;
figure;plot(t,X);xlabel('Time');ylabel('Amplitude');title('Input Function');
figure;plot(u,abs(x2));xlabel('frequency');ylabel('Amplitude');title('Magnitude Response of FRFT');
p = 301; % Number of time snapshots
fs = 600*10^6; % Sampling frequency
fc = 12*10^6; % Center frequency of narrowband sources
M = 180; % Number of array elements, i.e., sensors or antennas
N = 3; % Number of sources
sVar = 1; % Variance of the amplitude of the sources
mu=8*10^13;
T=1;
t=-T:0.01:T;
u=16;
%eqn 8 in the paper:
s1=sqrt(p)/fs;
alphad=atan(mu*s1.^2);
Y1=(1i*2*pi*fc*cos(alphad)*u);
C=1*cos(alphad)*sqrt(1+1i*tan(alphad))*exp(-1i*pi*fc*fc*sin(alphad)*cos(alphad));
s=C*exp(Y1);
doa = [70,72,110]; %DOAs
cSpeed = 3*10^8 ; % Speed of light
dist = 7; % Sensors (i.e., antennas) spacing in meters
% Constructing the Steering matrix
A = zeros(M, N);
for k = 1:N
    A(:, k) = exp(-1i*2*pi*fc*dist*cosd(doa(k))*(1/cSpeed)*[0:M-1]');   
end
noiseCoeff = 1; % Variance of added noise
x = A*s + sqrt(noiseCoeff)*randn(M, 1); % Sensor signals
% Estimating the covariance matrix of the sensor array %%%%%%%%%%%
R = (x*x')/p; % Empirical covariance of the antenna data
%  Finding the noise subspace and estimating the DOAs %%%%%%%%%%%%%
[V, D] = eig(R);
noiseSub = V(:, 1:M-N); % Noise subspace of R
theta = 0:1:180; %Peak search
a = zeros(M, length(theta));
res = zeros(length(theta), 1);
for i = 1:length(theta)
    a(:, i) = exp(-1i*2*pi*fc*dist*cosd(i)*(1/cSpeed)*[0:M-1]');
    res(i, 1) = 1/(norm(a(:, i)'*noiseSub).^2);
end
[resSorted, orgInd] = sort(res, 'descend');
DOAs = orgInd(1:N, 1);
DOAs1 = [];
DOAs2 = [];
DOAs3 = [];
p = 301; % Number of time snapshots
fs = 600*10^6; % Sampling frequency
fc = 12*10^6; 
M = 180; % Number of array elements, i.e., sensors or antennas
N = 3; % Number of sources
sVar = 1; % Variance of the amplitude of the sources
doa = [70,72,110]; %DOAs
cSpeed = 3*10^8 ; % Speed of light
mu=8*10^13;
dist = 7;
sVar = 1;
u=16;
s1=sqrt(p)/fs;
alphad=atan(mu*s1.^2);
Y1=(1i*2*pi*fc*cos(alphad)*u);
C=1*cos(alphad)*sqrt(1+1i*tan(alphad))*exp(-1i*pi*fc*fc*sin(alphad)*cos(alphad)*repmat([1:p]/fs, N, 1));
s=C*exp(Y1);
A = zeros(M, N);
for k = 1:N
    A(:, k) = exp(-1i*2*pi*fc*dist*cosd(doa(k))*(1/cSpeed)*[0:M-1]'); 
end    
noiseV = [0.2,0.1,100];
for i = 1:size(noiseV, 1)
    noiseCoeff = noiseV(i); % Variance of added noise
    
    x = A*s + sqrt(noiseCoeff)*randn(M, 1); % Sensor signals
    R = (x*x')/p; % Empirical covariance of the antenna data
    [V, D] = eig(R);
    noiseSub = V(:,1:M-N); % Noise subspace of R
    theta = 0:1:180; %Peak search
    a = zeros(M, length(theta));
    res = zeros(length(theta), 1);
    for j = 1:length(theta)
        a(:, j) = exp(-1i*2*pi*fc*dist*cosd(j)*(1/cSpeed)*[0:M-1]');
        res(j, 1) = 1/(norm(a(:, j)'*noiseSub).^2);
    end
    [resSorted, orgInd] = sort(res, 'descend');
    DOAs1 = [DOAs1,DOAs2,DOAs3 orgInd(1:N, 1)]; 
    figure
    plot(res);
    xlabel('Angle (deg)');
    ylabel('power');
    grid;

end
wavelength = 1; % normalized wavelength
d_0 = wavelength / 2;
designs = {...
    design_array_1d('nested', [2 12], d_0, 'Nested (2, 12)') ...
    design_array_1d('nested', [3  9], d_0, 'Nested (3, 9)') ...
};
n_designs = length(designs);

power_source = 3;
n_snaphots = 301;

n_grid = 180;
SNRs = linspace(0, 30, n_grid);

doas1 = linspace(-6, 6, 8);

MSEs_SNR_ana1 = zeros(n_designs, n_grid);
for dd = 1:n_designs
    design = designs{dd};
    A = steering_matrix(design, wavelength, doas1);
    for ii = 1:n_grid
        power_noise = power_source*10^(-SNRs(ii)/10);
        MSEs_SNR_ana1(dd, ii) = mean(ecov_coarray_mspice_1d(design, wavelength, ...
                doas1, power_source, power_noise, n_snaphots, 'DiagonalsOnly'));
    end
end
figure;
semilogy(SNRs, rad2deg(sqrt(MSEs_SNR_ana1)));
xlabel('SNR (dB)'); ylabel('RMSE (deg)'); grid on;

