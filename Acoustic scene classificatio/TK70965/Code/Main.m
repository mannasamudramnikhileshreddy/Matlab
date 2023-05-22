% for inp = 1:15
%     inp-1
[file,path] = uigetfile('*.*','Select Signal'); %select any .wav files from the browser
[y,fs1]=audioread(file);
fs = 8000;  %sampling frequency
sz = size(y);   %size
t = 0:(1/fs):((length(y)-1)/fs);
figure
plot(t,y);
xlabel('Time in secs');
ylabel('Amplitude in V');
title('Original Signal');

%log mel-band energies are extracted
s = melSpectrogram(y,fs);
figure
melSpectrogram(y,fs);
% NS(:,inp) = s(1:30,1);
title('Spectrogram of Original Signal');
% end

%%Classification using CNN
%Training
numFilters = 5;
filterSize = 5;
outputSize = 10;

layers = [convolution2dLayer(filterSize,numFilters)
    convolution2dLayer(filterSize,numFilters)
    convolution2dLayer(filterSize,numFilters)
    fullyConnectedLayer(outputSize)
    softmaxLayer];

load('labels.mat')
load('NS.mat')

%Testing
%for selecting of signal
for ii = 1:15
    if s(1:30,1) == NS(1:30,ii)
        tst = ii;
        tsts = s(1:30,1);
        break
    else
        tst = ii;
        tsts = s(1:30,1);
        break
    end
end

%prediction of signal
Y = pred(NS,tsts,labels)
accuracy = accu(NS,labels)
Yll = loglos(1)