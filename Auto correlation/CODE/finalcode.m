clc;clear all;close all;warning off;
B=7;
W=1;
N=64;
K=3;
n=K*N:B*N-1;
% prototypefilter
h1=1;
h2=0.911438;
h3= 0.411438;
factech=1+2*(h1+h2+h3);
hef(1:3*N)=0;
for i=1:3*N-1
   hef(1+i)=1-2*h1*cos(pi*i/(2*N))+2*h2*cos(pi*i/N)-2*h3*cos(pi*i*3/(2*N));
end

hef=hef/factech;
h=hef; 
%-------
Frame=1;
 y=zeros(1,3*N+(Frame-1)*N/2);
%  Generate mapped symbol data
 s=sign(randn(B,N))+1j*sign(randn(B,N)); % Data symbols, QPSK 
%     % Signal initialization
 SR=zeros(B,N);
 SI=zeros(B,N);
%=====generation of Sb blocks  (b=1,2,...7)   
    for nblock=1:B/W % Repeat for the number of subframes in a frame  
        for w=1:W % Repeat for the number of FBMC symbols in single subframe  
            SR((nblock-1)*W+w,:)=real(s((nblock-1)*W+w,:)); % Re
            SI((nblock-1)*W+w,:)=imag(s((nblock-1)*W+w,:)); % Im
         
        end
    end

s1=zeros(1,N);
s2=zeros(1,N);
s3=zeros(1,N);
s4=zeros(1,N); 

for ntrame=1:2
    for b=1:7
% OQAM Modulator
    if rem(ntrame,2)==1
      s1(b,1:N)=SR(b,N);
      s2(b,2:N)=j*SI(b,N);
    else
     s3(b,1:N)=j*SI(b,N);
     s4(b,2:N)=SR(b,N);
    end
    end
g=s1+s3;
x=ifft(g);
% Duplication of the signal
x1=[x x x];


end

ys=x1.*h;
y1=ys(3*N:7*N-1);
% y1=zeros(1,B*N);
% for k=1:K  %1:3
%   for L=1:N  %1:64
%     for b=1:B  %1:7
%         for n=K*N:B*N   %192:448
%             yf(1,n)=(real(x1((b-k)*N + L)).* h(L + k*N));
%             y1(1,n)=y1(1,n)+yf(1,n);
%         end 
%     end 
% end
% end   
g1=s2+s4;

%ppn2 network
xd=ifft(g1);

% Duplication of the signal
xd2=[xd xd xd];
% hadamard product in ppn2 network
y2=xd2.*h;
% %delay of N/2
% yx = delayseq(y2,N/2);
%steady state frame of ppn2 network
yd=y2((3*N)+N/2:(7*N-1)+N/2);
%sum of both ppn1 and ppn2 networks
y=y1+yd;
% yx=y1';
yr = awgn(y,20,'measured');
yfft = [zeros(1,length(y)/2) fft(y) zeros(1,length(y)/2)]
yTx = ifft(yfft);
yTxFreq = fft(yTx',500);
yTxFreqAbs = abs(yTxFreq);
yTxFreqAbs = yTxFreqAbs/max(yTxFreqAbs);
yTxFreqAbsPwr = 20*log(yTxFreqAbs);
figure,
subplot(3,2,1); plot(h);xlim([0 length(h)]);xlabel('filter prototype with K=3');
subplot(3,2,2); stem(real(s1));xlabel('real of s1');
subplot(3,2,3); stem(imag(s2));xlabel('imag of s3');
subplot(3,2,4); plot(real(x));xlabel('ifft of real of of s1+s3');
subplot(3,2,5),plot(imag(x));xlabel('ifft of real of of s1+s3');
subplot(3,2,6); plot(real(y1));xlabel('in ppn1 multiplying s1+s3 signal with h');
figure,
subplot(3,2,1); stem(imag(s3));xlabel('imag of s2');set(gca,'yticklabel',[]);
subplot(3,2,2); stem(real(s4));xlabel('real of s4');set(gca,'yticklabel',[]);
subplot(3,2,3); plot(real(xd));xlabel('ifft of real of of s2+s4');set(gca,'yticklabel',[]);
subplot(3,2,4),plot(imag(xd));xlabel('ifft of imag of of s2+s4');set(gca,'yticklabel',[]);
subplot(3,2,5); plot(real(y2));xlabel('in ppn2 multiplying s2+s4 signal with h');set(gca,'yticklabel',[]);
subplot(3,2,6); plot(real(yd));xlabel('delay signal in ppn2');set(gca,'yticklabel',[]);

%------final y
figure,
subplot(2,2,1); plot(real(y));xlabel('final signal');set(gca,'yticklabel',[]);
subplot(2,2,2); plot(yTxFreqAbs);xlabel('freqabs of tx signal');xlim([0 length(yTxFreqAbs)]);ylim([0 1]);set(gca,'yticklabel',[]);
subplot(2,2,3);plot(yTxFreqAbsPwr);xlabel('freqabspwr of tx signal' );xlim([0 length(yTxFreqAbsPwr)]);ylim([-100 0]);set(gca,'yticklabel',[]);
subplot(2,2,4); plot(real(yr));xlabel('signal send through awgn');set(gca,'yticklabel',[]);
%demodulator of fbmc signal
% Demodulator
% sdemod=zeros(N,Frame);
% for ntrame=1:2
%     r=yr.*h; % We apply the filter
%     u=zeros(1,N);
%     for k=1:Frame
%         u=u+r(1,1+(k-1)*N:k*N);                  % Summation
%     end
%     u=u.';
%     sest=fft(u);
%     sest=sest/0.6863;                           % Normalization 
%     % OQAM demodulation
%     if rem(ntrame,2)==1
%         sdemod(1:2:N,1)=real(sest(1:2:N));
%         sdemod(2:2:N,1)=j*imag(sest(2:2:N));
%     else
%        sdemod(1:2:N,1)=j*imag(sest(1:2:N));
%        sdemod(2:2:N,1)=real(sest(2:2:N));
%     end
% end
% % ysfft = [zeros(1,length( sdemod)/2) sdemod zeros(1,length( sdemod)/2)];
% ysRx = ifft( sdemod);
% ysRxFreq = fft(ysRx,500);
% ysRxFreqAbs = abs(ysRxFreq);
% ysRxFreqAbs = ysRxFreqAbs/max(ysRxFreqAbs);
% ysRxFreqAbsPwr = 20*log(ysRxFreqAbs);
% figure,
% subplot(3,2,1); plot(ysRxFreqAbs);xlim([0 length(ysRxFreqAbs)]);xlabel('frequency absoluteof received signal');ylim([0 1]);set(gca,'yticklabel',[]);
% subplot(3,2,2);plot(ysRxFreqAbsPwr); xlim([0 length(ysRxFreqAbsPwr)]);xlabel('frequency abs pwr received signal');ylim([-100 0]);set(gca,'yticklabel',[]);
% % subplot(3,2,3); plot(real( sdemod));xlabel('real of demod');set(gca,'yticklabel',[]);
% % subplot(3,2,4); plot(imag( sdemod));xlabel('imag of demod');set(gca,'yticklabel',[]);
% subplot(3,2,3); plot(( real(sdemod)));xlabel(' sdemod');set(gca,'yticklabel',[]);
% subplot(3,2,4); plot(( sdemod));xlabel(' sdemod');set(gca,'yticklabel',[]);
% % y=y';

% fe= acf1(y,63);
% subplot(3,2,5); plot((fe));xlabel(' sdemod');set(gca,'yticklabel',[])
% a=x1b(2).*h3+x1b(3).*h2+x1b(4).*h1;
% b=x1b(3).*h3+x1b(4).*h2+x1b(5).*h1;
% c=x1b(4).*h3+x1b(5).*h2+x1b(6).*h1;
% d=x1b(5).*h3+x1b(6).*h2+x1b(7).*h1;
% 
% a1=x2b(2).*h3+x2b(3).*h2+x2b(4).*h1;
% b2=x2b(3).*h3+x2b(4).*h2+x2b(5).*h1;
% c3=x2b(4).*h3+x2b(5).*h2+x2b(6).*h1;
% d4=x2b(5).*h3+x2b(6).*h2+x2b(7).*h1;
% -------modified FBMC signal----------%
% ym=zeros(1,K*N);
% ym(1,N)=[y(1,3*N),y(1,3*N+1)...
%          y(7*N+3*N-1)];

% yp=y1(3*N:4*N-1)
% L=0:N-1;
% t=[0:(b-K)*N-1];
%      
%      for k=2:3
%          for L=0:63 
%          if h_k(L)> delta
%              I_k(L,delta)=1;
%          else
%              I_k(L,delta)=0;
%          end
%          a_k(L)=(1/h_k(L))*I_k(L,delta);
%          end
%      end 
%     z(t)=yp*a(t);    
%     fy= acf1(z(t),63);           
%           Rz1(n,N)=-44*(1+h1(n)*a3(n));
%          figure,
%         subplot(4,2,1); plot((Rz1(n,N)));xlabel('Rz1(n,N)');set(gca,'yticklabel',[]) ;

% rstream = RandStream.create('mt19937ar','seed',2009);

Ntrial = 1e4;             % number of Monte-Carlo trials
snrdb = 10;                % SNR in dB
snr = db2pow(snrdb);      % SNR in linear scale
spower = 1;               % signal power is 1
npower = spower/snr;           % noise power
namp = sqrt(npower/2);    % noise amplitude in each channel
s = ones(1,Ntrial);       % signal  

x = y + yr;
mf = 1;
y = mf'*x; 
z = real(y);
Pfa = 0.1;
snrthreshold = db2pow(npwgnthresh(Pfa, 1,'coherent'));
mfgain = mf'*mf;
% To match the equation in the text above
% npower - N
% mfgain - M
% snrthreshold - SNR
threshold = sqrt(npower*mfgain*snrthreshold);
Pd = sum(z>threshold)/Ntrial
% x = n;
y = mf'*x;
z = real(y);
Pfa = sum(z>threshold)/Ntrial
figure;
hold on
plot(rocsnr(snrdb,'SignalType','NonfluctuatingCoherent','MinPfa',1e-4), 'r-.');


Ntrial = 1e4;             % number of Monte-Carlo trials
snrdb = 10.1;                % SNR in dB
snr = db2pow(snrdb);      % SNR in linear scale
spower = 1;               % signal power is 1
npower = spower/snr;           % noise power
namp = sqrt(npower/2);    % noise amplitude in each channel
s = ones(1,Ntrial);       % signal  

x = y + yr;
mf = 1;
y = mf'*x; 
z = real(y);
Pfa = 0.1;
snrthreshold = db2pow(npwgnthresh(Pfa, 1,'coherent'));
mfgain = mf'*mf;
% To match the equation in the text above
% npower - N
% mfgain - M
% snrthreshold - SNR
threshold = sqrt(npower*mfgain*snrthreshold);
Pd = sum(z>threshold)/Ntrial
% x = n;
y = mf'*x;
z = real(y);
Pfa = sum(z>threshold)/Ntrial
plot(rocsnr(snrdb,'SignalType','NonfluctuatingCoherent','MinPfa',1e-4), 'b-.');
hold off
% % plot(FPs, TPs, 'linewidth', 3, 'marker', 'o', 'markersize',10,'markeredgecolor', 'k', 'markerfacecolor', 'g');
% % hold on; 
% % plot(ThresholdRange, ThresholdRange, 'r-.', 'linewidth', 3);
% % axis([0,1,0,1]);
figure;
MaxSNR=-25;
MinSNR=0;
 NumPoints=-5
delta = [0.001 0.01];
rocpfa(delta,'SignalType','NonfluctuatingCoherent','MinSNR',-5);

NRdB = [-25 -20 -15 -10 -5 0];
[Pd,Pfa] = rocsnr(snrdb,'SignalType','NonfluctuatingCoherent');
figure;
semilogx(Pfa,Pd)
grid on
xlabel('P_{fa}')
ylabel('P_d')
legend('SNR -25 dB','SNR -20 dB','SNR -15 dB','SNR -10 dB',  'location','northwest')
title('Receiver Operating Characteristic (ROC) Curves')
% SNRs = (6:2:12);
% rocsnr(SNRs,'SignalType','NonfluctuatingNoncoherent')



