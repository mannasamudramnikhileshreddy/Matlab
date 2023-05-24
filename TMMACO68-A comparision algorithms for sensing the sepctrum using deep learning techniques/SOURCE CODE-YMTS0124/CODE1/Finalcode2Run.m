clc;
clear all;
close all;
warning off all
while(1)
        choice =menu('CR Energy Detect',...
            'System',...
            'DualPow_system',...
            'Pow mini',...
            'Pow medi',...
            'pow min mid',...
            'Pow maxx',...
            'Exit');
        if(choice ==7)
            break;
        end
        if(choice ==1)
            syst;
        end
        if(choice ==2)
            dou_snr;
        end
        if(choice ==3)
            num_sample = 128;
Eb_db = -20; 
Eb = 10.^(Eb_db./10); 
Pf = 0:0.05:1; 
pd=zeros(1,length(Pf));
pdd=zeros(1,length(Pf));
for m = 1:length(Pf)
    i = 0;ii=0;
for kk=1:10000 
 data_seq = randn(1,num_sample); 
 magin = sqrt(Eb).*randn(1,num_sample); 
 tx_sig = magin + data_seq; 
 energy = abs(tx_sig).^2; 
 fre_indx =(1/num_sample).*sum(energy); 
 thresh(m) = ((qfuncinv(Pf(m))./sqrt(num_sample))+ 1); 
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
end
figure;
plot(Pf, smooth(Pd),'ro-','Linewidth',1.5)
hold on
plot(Pf, smooth(Pdd),'bo-','Linewidth',1.5)
hold off
grid on
xlabel('prob of false alarm');
ylabel('prob of detect');
legend('Dynamic','const');
title('Power mini with N=128 of existing with SNR=-20dB');
        end
        if(choice ==4)
            num_sample = 128;
Eb_db = -10; 
Eb = 10.^(Eb_db./10); 
Pf = 0:0.05:1; 
pd=zeros(1,length(Pf));
pdd=zeros(1,length(Pf));
for m = 1:length(Pf)
    i = 0;ii=0;
for kk=1:10000 
 data_seq = randn(1,num_sample); 
 magin = sqrt(Eb).*randn(1,num_sample); 
 tx_sig = magin + data_seq; 
 energy = abs(tx_sig).^2; 
 fre_indx =(1/num_sample).*sum(energy); 
 thresh(m) = ((qfuncinv(Pf(m))./sqrt(num_sample))+ 1); 
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
end
figure;
plot(Pf, smooth(Pd),'ro-','Linewidth',1.5)
hold on
plot(Pf, smooth(Pdd),'bo-','Linewidth',1.5)
hold off
grid on
xlabel('prob of false alarm');
ylabel('prob of detect');
legend('Dynamic','const');
title('Power mid with N=128 of existing with SNR=-10dB');
        end
        if(choice ==6)
           num_sample = 128;
Eb_db = -2; 
Eb = 10.^(Eb_db./10); 
Pf = 0:0.05:1; 
pd=zeros(1,length(Pf));
pdd=zeros(1,length(Pf));
for m = 1:length(Pf)
    i = 0;ii=0;
for kk=1:10000 
 data_seq = randn(1,num_sample); 
 magin = sqrt(Eb).*randn(1,num_sample); 
 tx_sig = magin + data_seq; 
 energy = abs(tx_sig).^2; 
 fre_indx =(1/num_sample).*sum(energy); 
 thresh(m) = ((qfuncinv(Pf(m))./sqrt(num_sample))+ 1); 
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
end
figure;
plot(Pf, smooth(Pd),'ro-','Linewidth',1.5)
hold on
plot(Pf, smooth(Pdd),'bo-','Linewidth',1.5)
hold off
grid on

xlabel('prob of false alarm');
ylabel('prob of detect');
legend('Dynamic','const');
title('Power max with N=128 of existing with SNR=-2dB');
        end
         if(choice ==5)
           num_sample = 128;
Eb_db = -5; 
Eb = 10.^(Eb_db./10); 
Pf = 0:0.05:1; 
pd=zeros(1,length(Pf));
pdd=zeros(1,length(Pf));
for m = 1:length(Pf)
    i = 0;ii=0;
for kk=1:10000 
 data_seq = randn(1,num_sample); 
 magin = sqrt(Eb).*randn(1,num_sample); 
 tx_sig = magin + data_seq; 
 energy = abs(tx_sig).^2; 
 fre_indx =(1/num_sample).*sum(energy); 
 thresh(m) = ((qfuncinv(Pf(m))./sqrt(num_sample))+ 1); 
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
end
figure;
plot(Pf, smooth(Pd),'ro-','Linewidth',1.5)
hold on
plot(Pf, smooth(Pdd),'bo-','Linewidth',1.5)
hold off
grid on

xlabel('prob of false alarm');
ylabel('prob of detect');
legend('Dynamic','const');
title('Power min mid with N=128 of existing with SNR=-5dB');
        end
end