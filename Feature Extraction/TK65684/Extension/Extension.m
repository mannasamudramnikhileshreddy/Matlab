%%First Speech
%movefile('cv-invalid.mat','cv-invalid.csv');
%reading speech audio file
[y,fs1]=audioread('Railway.mp3');
%sound(y,fs1);                           %To hear the speech
y1 = y((100000:120000),1);
%Band Pass Filtering of The Signal
fpass = [fs1/4 fs1/2];
yb = bandpass(y1,fpass,fs1);

%order
p1 = 3;
%Linear Predictive Coding
[a1,g1] = lpc(yb,p1);                   %Output and error of LPC
x1 = yb(end-4096+1:end);
%filter estimate of lpc
est_x1 = filter([0 -a1(2:end)],1,x1);
%[acf,lags] = autocorr(x)

figure(1)
plot(1:100,x1(end-100+1:end),1:100,est_x1(end-100+1:end),'--')
grid
xlabel('Sample Number')
ylabel('Amplitude')
legend('Original signal','LPC estimate')

%error calculation
e1 = x1-est_x1;
[acs1,lags1] = xcorr(e1,'coeff');

figure(2)
plot(lags1,acs1)
grid
xlabel('Lags')
ylabel('Normalized Autocorrelation')
ylim([-0.1 1.1])

%%Second Speech
%movefile('cv-invalid.mat','cv-invalid.csv');
%reading speech audio file
[yy,fs2]=audioread('Airport.mp3');
%sound(yy,fs2);                           %To hear the speech
yy1 = yy((100000:120000),1);
%Band Pass Filtering of The Signal
fpass1 = [fs2/4 fs2/2];
yb1 = bandpass(y1,fpass1,fs2);

%order
p2 = 3;
%Linear Predictive Coding
[a2,g2] = lpc(yb1,p2);                   %Output and error of LPC
x2 = yb1(end-4096+1:end);
%filter estimate of lpc
est_x2 = filter([0 -a2(2:end)],1,x2);
%[acf,lags] = autocorr(x)

figure(3)
plot(1:100,x2(end-100+1:end),1:100,est_x2(end-100+1:end),'--')
grid
xlabel('Sample Number')
ylabel('Amplitude')
legend('Original signal','LPC estimate')

%error calculation
e2 = x2-est_x2;
[acs2,lags2] = xcorr(e2,'coeff');

figure(4)
plot(lags2,acs2)
grid
xlabel('Lags')
ylabel('Normalized Autocorrelation')
ylim([-0.1 1.1])
