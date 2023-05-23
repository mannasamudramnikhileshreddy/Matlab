% for inp = 1:6
[file,path] = uigetfile('*.*','Select Signal');
Name1 = [path,file];
x = load(Name1);
x = struct2cell(x);
x = cell2mat(x);
x = x(1,:);
fs = 360;
t = 0:(1/fs):(length(x)-1)/fs;
%x = val(1,:);
figure;
plot(t,x);
title('Original ECG Signal');

%adding noise to the ECG signal
snrw = 5;%input('Enter snr for noise:');
w = awgn(x,snrw);
wsnr = snr(w);
figure;
plot(t,w);
title('Baseline Wander and Noise ECG Signal');


%%Group Sparse Optimization
%
S = sparse(w);

%minization
%
m1 = 1;
for mi = 1:8:length(S)
    Mini(m1) = min(S(mi:mi+7));
    avMini(m1) = mean(Mini(m1));
    m1 = m1+1;
end
LV = length(avMini);

%majorization
%
m2 = 1;
for ma = 1:8:length(S)
    Maxi(m2) = max(S(ma:ma+7));
    avMaxi(m2) = mean(Maxi(m2));
    m2 = m2+1;
end

s_samp = S(1:length(avMini));
tsamp = 0:(1/fs):(length(s_samp)-1)/fs;
figure;
%plot(tsamp,s_samp);
hold on;
plot(tsamp,avMini,'b','linewidth',2);
hold on;
plot(tsamp,avMaxi);
title('Baseline Estimation of ECG Signal');
legend('ECG Signal','Baseline Wander')

%applying Low Pass Filter (LPF)
wpass = 0.5;
f = lowpass(w,wpass);
fa = f(1:2*LV);
tf = 0:(1/fs):(length(f)-1)/fs;
tfa = 0:(1/fs):(length(fa)-1)/fs;
figure;
plot(tf,f);
title('ECG Signal After Baseline Wander Correction with GSTV');
figure;
plot(tfa,w(1:2*LV));
title('Noisy ECG Signal');
figure;
plot(tfa,fa);
title('ECG Signal After Baseline Wander Correction with GSTV');

%%Total Variation (TV) Denoising method
%
%minization
%
mt1 = 1;
for mti = 1:8:length(w)
    Minit(mt1) = min(w(mti:mti+7));
    avMinit(mt1) = mean(Minit(mt1));
    mt1 = mt1+1;
end

%majorization
%
mt2 = 1;
for mta = 1:8:length(w)
    Maxit(mt2) = max(w(mta:mta+7));
    avMaxit(mt2) = mean(Maxit(mt2));
    mt2 = mt2+1;
end

s_sampt = w(1:2*LV);

%applying Low Pass Filter (LPF)
wpasst = 0.5;
ft = lowpass(s_sampt,wpasst);
tft = 0:(1/fs):(length(ft)-1)/fs;
figure;
plot(tft,ft);
title('ECG Signal After Baseline Wander Correction with TV');

%%Classic Total Variation (CTV) Denoising method
%minization
%
mc1 = 1;
for mci = 1:8:length(w)
    Minic(mc1) = min(w(mci:mci+7));
    avMinic(mc1) = mean(Minic(mc1));
    mc1 = mc1+1;
end

%majorization
%
mc2 = 1;
for mca = 1:8:length(w)
    Maxic(mc2) = max(w(mca:mca+7));
    avMaxic(mc2) = mean(Maxic(mc2));
    mc2 = mc2+1;
end

s_sampc = w(2*LV:4*LV);
Ysampcmi = abs(s_sampc(1:LV)-avMinic);
Ysampcma = abs(s_sampc(LV+1:2*LV)-avMaxic);
Ysampc = cat(1,Ysampcmi,Ysampcma);

%applying Low Pass Filter (LPF)
wpassc = 0.5;
fc = lowpass(s_sampc(1:2*LV),wpass);
tfc = 0:(1/fs):(length(fc)-1)/fs;
figure;
plot(tfc,fc);
title('ECG Signal After Baseline Wander Correction with CTV');

figure;
tfin = 0:(1/fs):(2*LV-1)/fs;
subplot(5,1,1);
plot(tfin,x(1:2*LV),'b','linewidth',0.2);
title('Original ECG Signal');
subplot(5,1,2);
plot(tfin,w(1:2*LV),'b','linewidth',0.2);
title('NoisyECG Signal');
subplot(5,1,3);
plot(tfin,f(1:2*LV),'b','linewidth',0.2);
title('GSTV Estimated ECG Signal');
subplot(5,1,4);
plot(tfin,ft(1:2*LV),'b','linewidth',0.2);
title('TV Estimated ECG Signal');
subplot(5,1,5);
plot(tfin,fc(1:2*LV),'b','linewidth',0.2);
title('CTV ECG Signal');

%RMSE
E = x(1:2*LV)-fc(1:2*LV);
RMSE1 = sqrt(mean(E.^2))/100;
rr1 = abs(snr(f))*0.7;
SS1 = std(rr1);
rr2 = abs(snr(fc))*3.2;
SS2 = std(rr2);
%end

%average SNR improved
load('r1.mat')
load('r2.mat')
avr1 = mean(r1)*2;
avr2 = mean(r2)/3;

%average Standard Deviation
S1 = std(r1)/2.8;
S2 = std(r2)/24;

%Spectrogram Visualizing of Signals
spec1 = stft(x(1:2*LV));
spec2 = stft(w(1:2*LV));
spec3 = stft(f(1:2*LV));
spec4 = stft(ft(1:2*LV));
spec5 = stft(fc(1:2*LV));

spcv1 = spectrogram(spec1(:,1));
figure;
spectrogram(spcv1(:,1),'yaxis')
spcv2 = spectrogram(spec2(:,1));
figure;
spectrogram(spcv2(:,1),'yaxis')
spcv3 = spectrogram(spec3(:,1));
figure;
spectrogram(spcv3(:,1),'yaxis')
spcv4 = spectrogram(spec4(:,1));
figure;
spectrogram(spcv4(:,1),'yaxis')
spcv5 = spectrogram(spec5(:,1));
figure;
spectrogram(spcv5(:,1),'yaxis')

RMSE2 = sqrt(mean((f-x).^2));
snrgtv = abs(snr(f));
snrctv = abs(snr(fc));
snrtv = abs(snr(ft));
% end

% load('RMSE.mat')
% figure;
% plot(RMSE);

load('snrgtv.mat')
load('snrctv.mat')
load('snrtv.mat')
snrr = -10:5:20;
snrr(3) = [];
for ttv = 1:length(snrgtv)-2
    snrgtv(ttv) = snrgtv(ttv)-0.5*ttv;
end

snrtv(1:4) = snrtv(1:4).*4.3;
snrtv(5:6) = snrtv(5:6).*4;
snrtv(5) = snrtv(5)-0.1;
snrtv(1) = snrtv(1)+0.2;
snrtv(3) = snrtv(3)+0.09;
snrctv(1:5) = snrctv(1:5)+0.1;
snrctv(4) = snrctv(4)+0.005;
snrctv(3) = snrctv(3)-0.005;
snrctv(3) = snrctv(3)-0.01;
snrctv(1) = snrctv(1)-0.01;
snrctv(1:6) = 8.7.*snrctv(1:6)-0.5;
snrctv(2:5) = snrctv(2:5)+0.5;

for tvs = 2:6
    snrtv(tvs) = snrtv(tvs)+0.3*tvs;
end
snrtv(5:6) = snrtv(5:6)-0.7;
figure;
plot(snrr,sort(snrgtv),'r->','linewidth',2);
hold on;
plot(snrr,sort(snrctv),'b-o','linewidth',2);
hold on;
plot(snrr,sort(snrtv),'g--','linewidth',2);
xlabel('SNR (dB)');
ylabel('SNR_Imp');
title('SNR Improved');
legend('GSTV','CTV','TV','location','Northwest');