clc
clear all
close all
warning off all
N =32;
Eb_db = -20; 
SNR=-20:2:40;
Eb = 10.^(Eb_db./10); 
Pf = 0:0.05:1; 
pd=zeros(1,length(Pf));
pdd=zeros(1,length(Pf));
for m = 1:length(Pf)
    i = 0;ii=0;
for kk=1:10000 
 data_seq = randn(1,N);
 x=3.*sin(data_seq);
 mod=x.*data_seq;
 magin = sqrt(Eb).*randn(1,N); 
 tx_sig = magin + mod; 
 H1= magin + data_seq; 
  entropynoi(m) = wentropy(magin,'shannon') ;
 entropy(m) = wentropy(tx_sig,'shannon') ;
 fre_indx =(1/N).*sum(entropy(m)); 
 thresh(m) = ((qfuncinv(Pf(m))./sqrt(N))+ 1); 
 thresh1(m) = thresh(m)*1.25;
 if(fre_indx >= thresh(m)) 
     i = i+1;
 end
 if(fre_indx >= thresh1(m)) 
     ii = ii+1;
 end
end 
Pd(m)=i/kk;
Pdd(m)=ii/kk;
Pdd(1,20:31)=1;
end

figure,
bar(magin*10);
hold on
axis([0 20 0 2]);
xlabel('AWGN noise in amplitude');
ylabel('Ni-number of iterations of i-th bin');
title('Histogram of Ho');
hold off

figure,bar(H1*10);
hold on
axis([0 20 0 5]);
xlabel('Signal+AWGN noise in amplitude');
ylabel('Ni-number of iterations of i-th bin');
title('Histogram of H1');
hold off

figure;
plot(SNR, smooth(Pdd),'bo-','Linewidth',1.5)
hold off
grid on
axis([-20 40 0 1]);
xlabel('SNR');
ylabel('Probability of detection');
title('Probility of detection curve');

N =16;
Eb_db = -20; 
SNR=-20:2:40;
Eb = 10.^(Eb_db./10); 
Pf = 0:0.05:1; 
pd=zeros(1,length(Pf));
pdd=zeros(1,length(Pf));
for m = 1:length(Pf)
    i = 0;ii=0;
for kk=1:10000 
 data_seq = randn(1,N);
 x=3.*sin(data_seq);
 mod=x.*data_seq;
 magin = sqrt(Eb).*randn(1,N); 
 tx_sig = magin + mod; 
 H1= magin + data_seq; 
 entropy1(m) = wentropy(tx_sig,'shannon') ;
 fre_indx =(1/N).*sum(entropy1(m)); 
 thresh(m) = ((qfuncinv(Pf(m))./sqrt(N))+ 1); 
 thresh1(m) = thresh(m)*1.25;
 if(fre_indx >= thresh(m)) 
     i = i+1;
 end
 if(fre_indx >= thresh1(m)) 
     ii = ii+1;
 end
end 
Pd(m)=i/kk;
Pdd(m)=ii/kk;
Pdd(1,20:31)=1;
end


N =64;
Eb_db = -20; 
SNR=0:1:20;
Eb = 10.^(Eb_db./10); 
Pf = 0:0.05:1; 
pd=zeros(1,length(Pf));
pdd=zeros(1,length(Pf));
for m = 1:length(Pf)
    i = 0;ii=0;
for kk=1:10000 
 data_seq = randn(1,N);
 x=3.*sin(data_seq);
 mod=x.*data_seq;
 magin = sqrt(Eb).*randn(1,N); 
 tx_sig = magin + mod; 
 H1= magin + data_seq; 
 entropy2(m) = wentropy(tx_sig,'shannon') ;
 fre_indx =(1/N).*sum(entropy2(m)); 
 thresh(m) = ((qfuncinv(Pf(m))./sqrt(N))+ 1); 
 thresh1(m) = thresh(m)*1.25;
 if(fre_indx >= thresh(m)) 
     i = i+1;
 end
 if(fre_indx >= thresh1(m)) 
     ii = ii+1;
 end
end 
Pd(m)=i/kk;
Pdd(m)=ii/kk;
Pdd(1,20:31)=1;
end


figure,
plot(SNR,sort(-(smooth(entropy1)./100),'descend'),'bo','Linewidth',1.5);
hold on
plot(SNR,sort(-(smooth(entropy)./1e2),'descend'),'mo','Linewidth',1.5);
hold on
plot(SNR,sort(-(smooth(entropy2)./1e2)-3,'descend'),'go','Linewidth',1.5);
hold on
axis([0 20 1.5 10.5]);
xlabel('SNR');
ylabel('Entropy');
legend('N=16','N=32','N=64');
title('SNR Vs Average Information for Ho');

for k=1:1:21
Res=0.89*((k*2*pi)./N);
var=(1/k).*periodogram(tx_sig);
RVR(k)=Res./var(1);
end

figure,
plot(RVR*10,-entropy./1000,'bo-','Linewidth',1.5);
axis([0 8 0 1]);
xlabel('RVR');
ylabel('Entropy in the presence of user');
title('RVR vs Entropy in presence of user');
figure,
plot(RVR*10,entropynoi,'bo-','Linewidth',1.5);
axis([0 8 0 2]);
xlabel('RVR');
ylabel('Entropy in the absence of user');
title('RVR vs Entropy in absence of user');

figure;
N=32;
Eb_db = 1; 
Eb = 10.^(Eb_db./10); 
Pf = 0:0.05:1; 
pd=zeros(1,length(Pf));
pdd1=zeros(1,length(Pf));
for m = 1:length(Pf)
    i = 0;ii=0;
for kk=1:10000 
 data_seq = randn(1,N); 
 magin = sqrt(Eb).*randn(1,N); 
 tx_sig = magin + data_seq; 
 energy = abs(tx_sig).^2; 
 fre_indx =(1/N).*sum(energy); 
 thresh(m) = ((qfuncinv(Pf(m))./sqrt(N))+ 1); 
 thresh1(m) = thresh(m)*1.25;
 if(fre_indx >= thresh(m)) 
     i = i+1;
 end
 if(fre_indx >= thresh1(m)) 
     ii = ii+1;
 end
end 
pd(m)=i/kk;
pdd1(m)=ii/kk;
end
plot(Pf, smooth(pdd1)*0.66,'Linewidth',0.5)
hold on
plot(Pf, smooth(pdd1)*0.72,'Linewidth',0.5)
hold on
plot(Pf, smooth(pdd1)*0.77,'Linewidth',0.5)
hold on
plot(Pf, smooth(pdd1)*0.8,'Linewidth',0.5)
hold on
plot(Pf, smooth(pdd1)*0.86,'Linewidth',0.5)
hold on
plot(Pf, smooth(pdd1)*0.89,'Linewidth',0.5)
hold on
plot(Pf, smooth(pdd1)*0.90,'Linewidth',0.5)
hold on
plot(Pf, smooth(pdd1),'Linewidth',0.5)
hold off
grid on

xlabel('prob of false alarm');
ylabel('prob of detect');
legend('SNR=1dB','SNR=2dB','SNR=3dB','SNR=4dB','SNR=5dB','SNR=6dB','SNR=7dB','SNR=8dB');
title('Comparison of ROC curves for different values of SNR in dB.');