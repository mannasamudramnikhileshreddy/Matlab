
clc
close all
clear all
L = 50;
snr_dB = -10; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf = 0.01:0.01:1; % Pf = Probability of False Alarm
%% Simulation to plot Probability of Detection (Pd) vs. Probability of False Alarm (Pf) 
    for m = 1:length(Pf)
    m;
    i = 0;
for kk=1:10000 % Number of Monte Caarlo Simulations
 n = randn(1,L); %AWGN noise with mean 0 and variance 1
  sig(kk) =var(n);
 s = sqrt(snr).*randn(1,L); % Real valued Gaussina Primary User Signal 
 y = s + n; % Received signal at SU
   P =var(y);
 energy = abs(y).^2; % Energy of received signal over N samples
 fin =(1/L).*sum(energy); % Test Statistic 
 thresh(m) = (qfuncinv(Pf(m))./sqrt(L))+ 1; % Theoretical value of Threshold
 if(fin >= thresh(m))  
     i = i+1;
 end
end
Pd(m) = i/kk; 
 end
thou=100;
l=thou:1:10000;
c1=P/2.*((P+sig.^2).*sig.^2);
c2=0.5.*log(sig.^2/P+sig.^2);
sum5=0;
for r=thou:l+1
%     zeta=abs(((l+2-r).*c2)/c1);
sum5=sum5+[gammainc(50,50)]./(gamma(l+2-r)/2);
end
sum6=0;
for j=thou:l
for r=thou:j
sum6=sum6+[gammainc(0.5,0.5)]./(gamma(j-r+1)/2);
% sum6=sum6+[gammainc((j-r+1)/2),zeta/2*P(l)+sig(l)^2]/(gamma(j-r+1)/2)
sum7=prod(sum6);
end
end
Pd1=(1-sum5).*(sum7);
sum8=0;
for l=thou:10000
sum8=abs(sum8)+abs(Pd1);
end
%analytical
subplot(1,3,1)
plot(Pf, Pd)

hold on
%% Theroretical ecpression of Probability of Detection; refer above reference.

thresh = (qfuncinv(Pf)./sqrt(L))+ 1;
Pd_the = qfunc(((thresh - (snr + 1)).*sqrt(L))./(sqrt(2).*(snr + 1)));
subplot(1,3,1)
plot(Pf, Pd_the, 'r')
% numerical
title('thou+20');
xlabel('Pf');
ylabel('Pd');
legend('analytical','numerical');
hold off
L = 200;
snr_dB = -10; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf1 = 0.01:0.01:1; % Pf = Probability of False Alarm
%% Simulation to plot Probability of Detection (Pd) vs. Probability of False Alarm (Pf) 
for m = 1:length(Pf1)
    m
    i = 0;
for kk=1:10000 % Number of Monte Caarlo Simulations
 n = randn(1,L); %AWGN noise with mean 0 and variance 1
 s = sqrt(snr).*randn(1,L); % Real valued Gaussina Primary User Signal 
 y = s + n; % Received signal at SU
 energy = abs(y).^2; % Energy of received signal over N samples
 fin =(1/L).*sum(energy); % Test Statistic 
 thresh(m) = (qfuncinv(Pf1(m))./sqrt(L))+ 1; % Theoretical value of Threshold
 if(fin >= thresh(m)) 
     i = i+1;
 end
end
Pd(m) = i/kk; 
end
%analytical
subplot(1,3,2)
plot(Pf1, Pd)
hold on
%% Theroretical ecpression of Probability of Detection; refer above reference.

thresh = (qfuncinv(Pf1)./sqrt(L))+ 1;
Pd_the = qfunc(((thresh - (snr + 1)).*sqrt(L))./(sqrt(2).*(snr + 1)));
subplot(1,3,2)
plot(Pf1, Pd_the, 'r')
% numerical
title('thou+40');
xlabel('Pf');
ylabel('Pd');
legend('analytical','numerical');
hold off
L = 1000;
snr_dB = -10; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf2 = 0.01:0.01:1; % Pf = Probability of False Alarm
%% Simulation to plot Probability of Detection (Pd) vs. Probability of False Alarm (Pf) 
for m = 1:length(Pf2)
    m
    i = 0;
for kk=1:10000 % Number of Monte Caarlo Simulations
 n = randn(1,L); %AWGN noise with mean 0 and variance 1
 s = sqrt(snr).*randn(1,L); % Real valued Gaussina Primary User Signal 
 y = s + n; % Received signal at SU
 energy = abs(y).^2; % Energy of received signal over N samples
 fin =(1/L).*sum(energy); % Test Statistic 
 thresh(m) = (qfuncinv(Pf2(m))./sqrt(L))+ 1; % Theoretical value of Threshold
 if(fin >= thresh(m)) 
     i = i+1;
 end
end
Pd(m) = i/kk; 
end
%analytical
subplot(1,3,3)
plot(Pf2, Pd)+
hold on
%% Theroretical ecpression of Probability of Detection; refer above reference.

thresh = (qfuncinv(Pf1)./sqrt(L))+ 1;
Pd_the = qfunc(((thresh - (snr + 1)).*sqrt(L))./(sqrt(2).*(snr + 1)));
subplot(1,3,3)
plot(Pf2, Pd_the, 'r')
% numerical
title('thou+60');
xlabel('Pf');
ylabel('Pd');
legend('analytical','numerical');
hold 

L = 5;
snr_dB = 3; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf3 = 0.01:0.01:1; % Pf = Probability of False Alarm
%% Simulation to plot Probability of Detection (Pd) vs. Probability of False Alarm (Pf) 
for m = 1:length(Pf3)
    m
    i = 0;
for kk=1:10000 % Number of Monte Caarlo Simulations
 n = randn(1,L); %AWGN noise with mean 0 and variance 1
 s = sqrt(snr).*randn(1,L); % Real valued Gaussina Primary User Signal 
 y = s + n; % Received signal at SU
 energy = abs(y).^2; % Energy of received signal over N samples
 fin =(1/L).*sum(energy); % Test Statistic 
 thresh(m) = (qfuncinv(Pf3(m))./sqrt(L))+ 1; % Theoretical value of Threshold
 if(fin >= thresh(m))  
     i = i+1;
 end
end
Pd(m) = i/kk; 
end
%analytical
figure,
subplot(1,3,1)
plot(Pf3, Pd)

hold on
%% Theroretical ecpression of Probability of Detection; refer above reference.

thresh = (qfuncinv(Pf3)./sqrt(L))+ 1;
Pd_the = qfunc(((thresh - (snr + 1)).*sqrt(L))./(sqrt(2).*(snr + 1)));
subplot(1,3,1)
plot(Pf3, Pd_the, 'r')
% numerical
title('thou+20');
xlabel('Pf');
ylabel('Pd');
legend('analytical','numerical');
hold off
L = 10;
snr_dB = 3; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf3 = 0.01:0.01:1; % Pf = Probability of False Alarm
%% Simulation to plot Probability of Detection (Pd) vs. Probability of False Alarm (Pf) 
for m = 1:length(Pf3)
    m
    i = 0;
for kk=1:10000 % Number of Monte Caarlo Simulations
 n = randn(1,L); %AWGN noise with mean 0 and variance 1
 s = sqrt(snr).*randn(1,L); % Real valued Gaussina Primary User Signal 
 y = s + n; % Received signal at SU
 energy = abs(y).^2; % Energy of received signal over N samples
 fin =(1/L).*sum(energy); % Test Statistic 
 thresh(m) = (qfuncinv(Pf3(m))./sqrt(L))+ 1; % Theoretical value of Threshold
 if(fin >= thresh(m))  
     i = i+1;
 end
end
Pd(m) = i/kk; 
end
%analytical

subplot(1,3,2)
plot(Pf3, Pd)

hold on
%% Theroretical ecpression of Probability of Detection; refer above reference.

thresh = (qfuncinv(Pf3)./sqrt(L))+ 1;
Pd_the = qfunc(((thresh - (snr + 1)).*sqrt(L))./(sqrt(2).*(snr + 1)));
subplot(1,3,2)
plot(Pf3, Pd_the, 'r')
% numerical
title('thou+40');
xlabel('Pf');
ylabel('Pd');
legend('analytical','numerical');
hold off
L = 100;
snr_dB = 3; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf3 = 0.01:0.01:1; % Pf = Probability of False Alarm
%% Simulation to plot Probability of Detection (Pd) vs. Probability of False Alarm (Pf) 
for m = 1:length(Pf3)
    m
    i = 0;
for kk=1:10000 % Number of Monte Caarlo Simulations
 n = randn(1,L); %AWGN noise with mean 0 and variance 1
 s = sqrt(snr).*randn(1,L); % Real valued Gaussina Primary User Signal 
 y = s + n; % Received signal at SU
 energy = abs(y).^2; % Energy of received signal over N samples
 fin =(1/L).*sum(energy); % Test Statistic
 thresh(m) = (qfuncinv(Pf3(m))./sqrt(L))+ 1; % Theoretical value of Threshold
 if(fin >= thresh(m)) 
     i = i+1;
 end
end
Pd(m) = i/kk; 
end
%analytical

subplot(1,3,3)
plot(Pf3, Pd)

hold on
%% Theroretical ecpression of Probability of Detection; refer above reference.

thresh = (qfuncinv(Pf3)./sqrt(L))+ 1;
Pd_the = qfunc(((thresh - (snr + 1)).*sqrt(L))./(sqrt(2).*(snr + 1)));
subplot(1,3,3)
plot(Pf3, Pd_the, 'r')
% numerical
title('thou+60');
xlabel('Pf');
ylabel('Pd');
legend('analytical','numerical');
hold on





