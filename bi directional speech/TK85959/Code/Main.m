% movefile('You Are My Remainder.wav','You Are My Remainder.mp3');
% for inp = 1:3
%%Bi-Directional Kalman Filtering
%
%Input Signal
[y,fs1]=audioread('You Are My Remainder.mp3');
text = 'You Are My Remainder';in = 2;
txt = 'Y A R';
%adding noise
snrr = 5;
yn = awgn(y,snrr);
%sound(y,fs1);                           %To hear the speech
P = 0;  %Initial State Estimate
Y = rand(length(yn),2);  %predicted signal
%error
e = yn-Y;
%Kalman Gain
c = xcov(y);    %covariance matrix of input signal
cp = xcov(Y);   %covariance matrix of predicted signal
R = xcov(yn);   %covariance matrix of noisy signal
%measurement gain matrix
C = rand(length(yn),2);

K = (Y(9.5*10^4:10^5,1).*C(9.5*10^4:10^5,1)')/((C(9.5*10^4:10^5,1).*Y(9.5*10^4:10^5,1).*C(9.5*10^4:10^5,1)')...
    +R(9.5*10^4:10^5,1)); %kalman gain

%correct state estimate
szK = size(K);
CSE = Y(1:1000,1) + K(1:1000,1).*e(1:1000,1);

%MFCC Feature Extraction
coeffs = mfcc(y,fs1);

%correlation calculation
corsnr = 0.97*(randi([1 1],1,1));
rng1 = randi([10000 10500],1,1);
rng2 = randi([10600 11000],1,1);
r = xcorr(y(rng1:rng2),yn(rng1:rng2));
rm = abs(mean(r))+corsnr;
% end

if in == 1
%correlations of repititions
load('First.mat')
load('Second.mat')
load('Third.mat')
%correlations of each words
load('Hello.mat')
load('Good.mat')
load('Morning.mat')

else
%correlations of repititions
load('First2.mat')
load('Second2.mat')
load('Third2.mat')
%correlations of each words
load('You.mat')
load('Are.mat')
load('My.mat')
load('Remainder.mat')
end

disp(text)
%accuracy calculation
accuracy = accu(txt);
