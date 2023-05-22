% for inp = 1:3
%%selection of signal
[file,path] = uigetfile('*.*','Select Signal');
Name1 = [path,file];
selinp = erase(file,".mat");
x = load(Name1);
x = struct2cell(x);
x = cell2mat(x);
% xini(inp,:) = x(1:100,1);
% end
% x = table2array(x{1,1});
fs = 256;   %sampling frequency
t = 0:(1/fs):(length(x)-1)/fs;
figure;
plot(t,x);
xlabel('Time in secs');
ylabel('Amplitude in mV');
title('Original EMG Signal');

%%Feature Extraction
%Calculation of RMS Value
rm1 = 1;
for rm = 1:length(x)
    if isnan(x(rm)) == 0
        y(rm1) = x(rm);
        rm1 = rm1+1;
    end
end

RMSY = rms(y);
% RMSYi(inp) = RMSY;

%Calculation of Mean Absolute Value (MAV)
MAVY = mad(y);
% MAVYi(inp) = MAVY;

%Calculation of Zero Crossings
for zc = 1:length(y)
    if sign(y(zc)) >= 0
        ZC(zc) = sign(y(zc));
    elseif sign(y(zc)) < 0
        ZC(zc) = sign(y(zc));
    else
        ZC(zc) = sign(y(zc));
    end
end
% ZCi(inp,1:50) = ZC(1:50);

%Calculation of Waveform Length
for wl = 1:length(y)
    if max(y) == y(wl)
        wlmx = wl;
    elseif min(y) == y(wl)
        wlmn = wl;
    end
    if wl == length(y)
        WL = wlmx+wlmn;
    end
end
% WLi(inp) = WL;
% end

%%Classification with CNN
%Training

load('labels')
load('MAVYi')
load('RMSYi')
load('WLi')
load('ZCi')

N = 49;
Xi = cat(2,MAVYi,RMSYi,WLi);
features = reshape(Xi,9,1);

filterSize = 5;
numFilters = 5;
inputsize = [5 5];
outputSize = 5;
numFeatures = 5;

layers = [imageInputLayer(inputsize);
        convolution2dLayer(filterSize,numFilters);
        fullyConnectedLayer(outputSize);
        reluLayer;
        batchNormalizationLayer;
        dropoutLayer;
        softmaxLayer;
        classificationLayer];
    
options = trainingOptions('adam',...
        'Maxepochs',5,...
        'GradientThreshold',0.01,...
        'InitialLearnRate',0.0001,...
        'LearnRateSchedule','piecewise',...
        'LearnRateDropPeriod',125,...
        'LearnRateDropFactor',0.2,...
        'Verbose',0,...
        'Plots','training-progress');
    
%Testing
if selinp == 'flex13'
    init = 1;
elseif selinp == 'flex21'
    init = 2;
elseif selinp == 'EMG111'
    init = 3;
end

%prediction
out = predc(init)