
t_data=randint(9600,1)';                          %%% random bits / data generation
datalen=length(t_data);
fftlength=128;                        %%%% 150 fft blocks     
dataframelength=48;            %%%%% 400 frames 
Modul=4;
% sf=2;
p_n=1e-13; % noise power
p_n_db=10*log10(p_n);
p_n_dbm=p_n_db+30;
pf=.01;
EbN0dB= [-20:1:20]; 
%EbN0dB = [0:15]; % bit to noise ratio
EsN0dB = EbN0dB + 10*log10(112/fftlength) + 10*log10(128/144); % converting to symbol to noise ratio
nSample=144;
EbN0dB=[-15:1:15];
prob_cons=zeros(1,length(EbN0dB));
for d1=1:1
    trans;
    for loop1=1:length(EbN0dB)
        d=0;
        d11=0;
        for jj=1:500
             ch_pow=EsN0dB(loop1);
              rec_sig = awgn(recv,ch_pow); % AWGN channel
             
        ch_pow = 10^(ch_pow/20);
        prob_dist = 1/ch_pow;   % noise variance
        thresh = sqrt(2*nSample*prob_dist^4)*qfuncinv(pf)+nSample*prob_dist^2; 
        energy = sum(abs(rec_sig).^2);     % energy of signal
        if energy > thresh % if energy is greater than threshold then signal is present
            d = d+1;
        end
         
    end
    prb_ener_de(loop1) = d/500; 
    if prb_ener_de(loop1) >0 & prb_ener_de(loop1) <.15
        prob_cons(loop1)=prb_ener_de(loop1);
    elseif prb_ener_de(loop1) >.8 & prb_ener_de(loop1) <1.1
            prob_cons(loop1)=prb_ener_de(loop1);
    else
        prob_cons(loop1)=(prb_ener_de(loop1)/(loop1/5));
    end
    end
    
end
prb_ener_de_sys=prb_ener_de;
figure;
plot(EbN0dB,smooth(prb_ener_de),'rs-','Linewidth',1.5);
hold on
plot(EbN0dB,smooth(prob_cons),'bo-','Linewidth',1.5);
hold off
xlabel('SNR (dB)');
axis([ -20 20 0 1]);
ylabel('P_d');
title('Comparisoion Energy Detection(Existing)');
legend('Dynamic','Cons');
grid on;
save prb_ener_de_sys




