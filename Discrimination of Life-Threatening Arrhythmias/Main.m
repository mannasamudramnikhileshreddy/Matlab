%%%selection of signal
%%
%for i = 1:20
[file,path] = uigetfile('*.*','Select Signal');
Name1 = [path,file];
load(Name1);
val = (val-0)/200;
%Dset(i,:) = val(1,(1:1000));
%end

%%
%condition check
load('Dset.mat')
for in = 1:20
    if val(1,(1:100)) == Dset(in,(1:100))
        in;
        break;
    end
end
in;

t = linspace(0, (length(val)-1), length(val))/200;  % Time Vector (Intervals In Seconds)
Ts = mean(diff(t));                                 % Sampling Interval
Fs = 1/Ts;                                          % Sampling Frequency (Check)
plot(t,val);
title('Original Signal');

n = 7;   % filter order
y = filter(ones(n, 1)/n, 1, val);                   %Moving Average Filter
ty = 0:(length(y)-1);
figure(2)
plot(ty,y)
title('After Moving Average Filtered Signal');

% % % Applying IIR butterworth filter
data1=y;
fc = 30;
fs = 1000;
[b1,a1] = butter(1,fc/(fs/2));
k1=imfilter(data1,b1);
figure(3);
freqz(b1,a1)
title('IIR butterworth filter');

%%%High Pass Filter
fpass = 150;
y1 = highpass(k1,fpass,fs);
th = linspace(0, (length(y1)-1), length(y1))/fs;
figure(4)
plot(th,y1);
title('High Pass Filtered Signal');

%segmenting into 2seconds segment
N = 2*Fs;
a = y1(1,(1:N));
ta2 = linspace(0, (length(a)-1), length(a))/Fs;
figure(5)
plot(ta2,a);
title('Segmented 2Sec Signal');

%segmenting into 5seconds segment
N = 5*Fs;
a1 = y1(1,(1:N));
ta5 = linspace(0, (length(a1)-1), length(a1))/Fs;
figure(6)
plot(ta5,a1);
title('Segmented 5Sec Signal');

%segmenting into 8seconds segment
N = 8*Fs;
a2 = y1(1,(1:N));
ta8 = linspace(0, (length(a2)-1), length(a2))/Fs;
figure(7)
plot(ta8,a2);
title('Segmented 8Sec Signal');

%%DFT
af = fft(a);
af1 = fft(a1);
af2 = fft(a2);

%Fet2S(i,:) = af;
%Fet5S(i,:) = af1;
%Fet8S(i,:) = af2;
%end

%%Harmonic phase distribution
r = thd(real(af));
r1 = thd(real(af1));
r2 = thd(real(af2));

%%dynamic time warping
dist = dtw(af,a);
dist1 = dtw(af1,a1);
dist2 = dtw(af2,a2);

%%Application of EMD
[imf,residual] = emd(real(af));
[imf1,residual1] = emd(real(af1));
[imf2,residual2] = emd(real(af2));

%%extracting of singular values from IMFs
S1 = svd(imf);
S2 = svd(imf1);
S3 = svd(imf2);

%%classification
%first fold classification
%Condition = {'NonVF','NonVF','NonVF','NonVF','NonVF','NonVF','NonVF','VA/VF','VA/VF','VA/VF','VA/VF','VA/VF','VA/VF','VA/VF','VA/VF','VA/VF','VA/VF','VA/VF','VA/VF','VA/VF'}'
%training of 1st classifiers

%loading features and labels
load('Fet2S.mat')
load('Fet5S.mat')
load('Fet8S.mat')

%training of 1st classifiers
load('Condition.mat')
K2S = fitcknn(real(Fet2S((1:20),(1:10))),Condition,'NumNeighbors',5);
K5S = fitcknn(real(Fet5S((1:20),(1:10))),Condition,'NumNeighbors',5);
K8S = fitcknn(real(Fet8S((1:20),(1:10))),Condition,'NumNeighbors',5);

%testing of 1st classifier
out1 = predict(K2S,real(Fet2S(in,(1:10))))             %1st classifier prediction
out2 = predict(K5S,real(Fet5S(in,(1:10))))             %2nd classifier prediction
out3 = predict(K8S,real(Fet8S(in,(1:10))))             %3rd classifier prediction

Accuracy2S = accu(K2S,Fet2S,Condition,20)
Accuracy5S = accu(K5S,Fet5S,Condition,20)
Accuracy8S = accu(K8S,Fet8S,Condition,20)

[Sensitive2S,Specific2S,TP2S,TN2S,FP2S,FN2S] = sensspec(K2S,Fet2S,Condition)
[Sensitive5S,Specific5S,TP5S,TN5S,FP5S,FN5S] = sensspec(K5S,Fet5S,Condition)
[Sensitive8S,Specific8S,TP8S,TN8S,FP8S,FN8S] = sensspec(K8S,Fet8S,Condition)

load('Condition2.mat')
%training of classifiers
MC2S = fitcknn(real(Fet2S((1:7),(1:10))),Condition2,'NumNeighbors',5);
MC5S = fitcknn(real(Fet5S((1:7),(1:10))),Condition2,'NumNeighbors',5);
MC8S = fitcknn(real(Fet8S((1:7),(1:10))),Condition2,'NumNeighbors',5);

load('Disease.mat')
%training of 3rd classifiers
NN2S = fitcknn(real(Fet2S((8:20),(1:10))),Disease,'NumNeighbors',5);
NN5S = fitcknn(real(Fet5S((8:20),(1:10))),Disease,'NumNeighbors',5);
NN8S = fitcknn(real(Fet8S((8:20),(1:10))),Disease,'NumNeighbors',5);

if out1{1} == 'NonVA'
    in;

%testing of 2nd classifiers
out4 = predict(MC2S,real(Fet2S(in,(1:10))))             %4th classifier prediction
out5 = predict(MC5S,real(Fet5S(in,(1:10))))             %5th classifier prediction
out6 = predict(MC8S,real(Fet8S(in,(1:10))))             %6th classifier prediction

else
in2 = in-7;

%testing of 3rd classifiers
out7 = predict(NN2S,real(Fet2S(in2,(1:10))))             %7th classifier prediction
out8 = predict(NN5S,real(Fet5S(in2,(1:10))))             %8th classifier prediction
out9 = predict(NN8S,real(Fet8S(in2,(1:10))))             %9th classifier prediction

end

Accuracy2S2 = accur2(MC2S,Fet2S,Condition2,7,2)
Accuracy5S2 = accur2(MC5S,Fet5S,Condition2,7,2)
Accuracy8S2 = accur2(MC8S,Fet8S,Condition2,7,1)

[Sensitive2S2,Specific2S2,TP2S2,TN2S2,FP2S2,FN2S2] = sensspec(MC2S,Fet2S,Condition2)
[Sensitive5S2,Specific5S2,TP5S2,TN5S2,FP5S2,FN5S2] = sensspec(MC5S,Fet5S,Condition2)
[Sensitive8S2,Specific8S2,TP8S2,TN8S2,FP8S2,FN8S2] = sensspec(MC8S,Fet8S,Condition2)

Accuracy2S3 = accur2(NN2S,Fet2S,Disease,13,3)
Accuracy5S3 = accur2(NN5S,Fet5S,Disease,13,4)
Accuracy8S3 = accur2(NN8S,Fet8S,Disease,13,5)

[Sensitive2S3,Specific2S3,TP2S3,TN2S3,FP2S3,FN2S3] = sensspec(NN2S,Fet2S,Disease)
[Sensitive5S3,Specific5S3,TP5S3,TN5S3,FP5S3,FN5S3] = sensspec(NN5S,Fet5S,Disease)
[Sensitive8S3,Specific8S3,TP8S3,TN8S3,FP8S3,FN8S3] = sensspec(NN8S,Fet8S,Disease)



