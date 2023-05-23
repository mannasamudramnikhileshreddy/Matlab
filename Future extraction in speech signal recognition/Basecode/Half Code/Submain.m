%%First Speech
%movefile('cv-invalid.mat','cv-invalid.csv');
%reading speech audio file
[y,fs1]=audioread('Railway.mp3');
%sound(y,fs1);                           %To hear the speech
y1 = y((100000:120000),1);

%%Feature Extraction
%order
p1 = 3;
%Linear Predictive Coding
[a1,g1] = lpc(y1,p1);                   %Output and error of LPC
x1 = y1(end-4096+1:end);
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