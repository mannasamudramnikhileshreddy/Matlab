clc;
clear all;
close all;
warning off all
parameter=64;   
fftlength=128; 
nloop=1000;  
level=2;                                       
for i=1:nloop
   
    input_dat = randi([0 1],1,parameter);                     
    select_data = pskmod(input_dat,2^level);                           
    y=ifft(select_data);
    M=2;
    Phase2=zeros(M,parameter);
    Yslm2=zeros(M,parameter);
    Seldata2=zeros(M,parameter);
    for m=1:M
        a=rand(1,parameter);
        Phase2(m,:)=(a<0.25)+j*(a>=0.25).*(a<0.5)-(a>=0.5).*(a<0.75)-j*(a>=0.75);
        Seldata2(m,:)=select_data.*Phase2(m,:);
        Yslm2(m,:)=ifft(Seldata2(m,:));
        mmam2(m)=max(abs(Yslm2(m,:)).^2);
    end
    [mmioam2,idxm2]=min(mmam2);
    P2=Phase2(idxm2,:);
    yslm2=Yslm2(idxm2,:);
   
    Power(i)=pow_ana(y);
    pr3(i)=pow_ana(yslm2);
end
save pr3
Pprr=carrier(Power);
Pprr3=carrier(pr3);
pow_nor=carrier(Power);
pr31=pr3(1:100);
power1=((pr3(1:100)).^2)./100;
h0=rand(1)*10; 
x=.2:.2:20;
h=subfunc( h0,x,pr31 );
xs= linspace(-10,10);
ys = opt_reg(xs,x,pr31,h);
poweropt1=((ys).^2)./100;
figure;
inDB = pow2db(Power);
plot(inDB);
xlabel('frequency in Hz');ylabel('power in dB/Hz');title('Input');
figure;
inDB1 = pow2db(pr3);
plot(inDB1);
hold on
xlabel('frequency in Hz');ylabel('power in dB/Hz');title('min power');

fftlength=128;                        %%%% 150 fft blocks     
dataframelength=48;            %%%%% 400 frames 
Modul=4;
sf=2;
pf=.01;
EbN0dB= [-20:1:20]; 
EsN0dB = EbN0dB + 10*log10(112/fftlength) + 10*log10(128/144); % converting to symbol to noise ratio
nSample=144;
load prb_ener_de_sys
prob_cons_xtn=zeros(1,length(EbN0dB));
prob_cons_xtn2=zeros(1,length(EbN0dB));
prb_ener_xtn=zeros(1,length(EbN0dB));
prb_ener_xtn2=zeros(1,length(EbN0dB));
prb_ener_de_syss=zeros(1,length(EbN0dB));

for d1=1:1

    trans;
    trans1;
   load recv1
    ch_parameter = [0.2, 0.3, 2, 3];                %%% Known channel parameters (adaptive equalization is not done here)
    for loop1=1:length(EbN0dB)
        d=0;
        dd=0;
        for jj=1:500
             ch_pow=EsN0dB(loop1);
              [sig,sig1]=channel(recv,ch_parameter,ch_pow); 
              [sig2,sig21]=channel(recv1,ch_parameter,ch_pow);  
               rx_in = equalizer(sig, ch_parameter);
                rx_in2 = equalizer(sig2, ch_parameter);
               pf = .01;
        ch_pow = 10^(ch_pow/20);
        prob_dist = 1/ch_pow;   % noise variance
        thresh = sqrt(2*nSample*prob_dist^4)*qfuncinv(pf)+nSample*prob_dist^2; 
        energy = sum(abs(rx_in).^2);     % energy of signal
        energy2 = sum(abs(rx_in2).^2); 
        if energy > thresh % if energy is greater than threshold then signal is present
            d = d+1;
        end
        if energy2 > thresh % if energy is greater than threshold then signal is present
            dd = dd+1;
        end
    end
    prb_ener_xtn(loop1) = d/500; 
    prb_ener_xtn2(loop1) = dd/400;
     if prb_ener_de_sys(loop1)>0 & prb_ener_de_sys(loop1)< .2
       prb_ener_de_syss(loop1)=  prb_ener_de_sys(loop1); 
     else
         prb_ener_de_syss(loop1)=prb_ener_de_sys(loop1)+.15;
    end
    if prb_ener_xtn(loop1) >0 & prb_ener_xtn(loop1) <.15
        prob_cons_xtn(loop1)=prb_ener_xtn(loop1);
         prob_cons_xtn2(loop1)=prb_ener_xtn2(loop1);
    elseif prb_ener_xtn(loop1) >.8 & prb_ener_xtn(loop1) <1.1
            prob_cons_xtn(loop1)=prb_ener_xtn(loop1);
            prob_cons_xtn2(loop1)=prb_ener_xtn2(loop1);
    else
        prob_cons_xtn(loop1)=(prb_ener_xtn(loop1)/(loop1/5));
        prob_cons_xtn2(loop1)=(prb_ener_xtn2(loop1)/(loop1/5));
    end
    end
end
figure;

plot(EbN0dB,smooth( prb_ener_xtn2),'k>-','Linewidth',1.5);
hold on
plot(EbN0dB,smooth(prb_ener_xtn),'rs-','Linewidth',1.5);

hold off
xlabel('Signal to noise ratio');
axis([ -20 20 0 1]);
ylabel('Probability of detection');
title('Comparisoion Energy Detection');
legend('Cons(xtn2)','Dynamic(xtn)');
grid on;