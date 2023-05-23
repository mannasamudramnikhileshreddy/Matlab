clc;
clear all;
close all;
warning off;

% for inp = 1:10
[file,path] = uigetfile('*.*','Select Signal');
Name1 = [path,file];
a1 = load(Name1);
a1 = struct2cell(a1);
a1 = cell2mat(a1);
a1 = a1(1,:);
% P(inp,1:1000) = a1(1,1:1000);
% end
fs = 250;
t = 0:(1/fs):(length(a1)-1)/fs;
figure;
plot(t,a1);
xlabel('time');
ylabel('Amplitude');
title('ECG signal');

%Pre-processing
[imf,residual] = emd(a1);
ti = 0:(1/fs):(length(imf)-1)/fs;
figure;
subplot(4,1,1);
plot(ti,imf(:,1));
xlabel('time');
ylabel('Amplitude');
subplot(4,1,2);
plot(ti,imf(:,2));
xlabel('time');
ylabel('Amplitude');
subplot(4,1,3);
plot(ti,imf(:,3));
xlabel('time');
ylabel('Amplitude');
subplot(4,1,4);
plot(ti,imf(:,4));
xlabel('time');
ylabel('Amplitude');
title('IMFs of ECG Signal');

%Application of LPF
fpass = fs/2;
yl = lowpass(a1,fpass,fs);
tl = 0:(1/fs):(length(yl)-1)/fs;
figure;
plot(tl,yl);
xlabel('time');
ylabel('Amplitude');
title('ECG signal after Application of LPF');

%Application of HPF
fpassh = 30;
yh = highpass(a1,fpassh,fs);
th = 0:(1/fs):(length(yh)-1)/fs;
figure;
plot(th,yh);
xlabel('time');
ylabel('Amplitude');
title('ECG signal after Application of HPF');

%Application of DPF
Nf = 50; 
Fpass = 100; 
Fstop = 120;

d = designfilt('differentiatorfir','FilterOrder',Nf, ...
    'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
    'SampleRate',fs);

fvtool(d,'MagnitudeDisplay','zero-phase','fs',fs)

%Feature Extraction
%Peak Detection
[pks,locs] = findpeaks(a1);
findpeaks(a1)
title('Peak Detection');

%ECG Signal Segmentation and Wave Modeling
s1 = a1(1,locs(1):locs(1)+3);s2 = a1(1,locs(2):locs(2)+3);s3 = a1(1,locs(3):locs(3)+3);
s = a1(60:100);netout  = 1;
seg = cat(2,s1,s2,s3);
figure;
plot(th,yh);
xlabel('time');
ylabel('Amplitude');
title('Segmented ECG signal');

ts = 0:(1/fs):(length(s)-1)/fs;
figure;
plot(ts,s);
xlabel('time');
ylabel('Amplitude');
title('Segment of an ECG signal');

%Classification
%Training
load('P.mat')
load('labels.mat')
Pin = P(1,1:100);lab = labels(1);
net = newff([0 99],[5 1],{'tansig' 'purelin'});

Y = sim(net,Pin);

%Testing 
[Ypred,Yout] = predic(netout)

%accuracy
accuracy = accur(netout)