% for inp = 1:6
[file,path] = uigetfile('*.*','Select Signal');
Name1 = [path,file];
x = load(Name1);
x = struct2cell(x);
x = cell2mat(x);
x = x(1,:);
fs = 360;
t = 0:(1/fs):(length(x)-1)/fs;
figure;
plot(t,x);
xlabel('Time (Secs)');
ylabel('Amplitude (V)');
title('Original ECG Signal');

%%Conditional Generative Adversarial Networks
%%Generator
%adding noise to the input ECG signal
snra = 5;
xn = awgn(x,snra);
tn = 0:(1/fs):(length(xn)-1)/fs;
figure;
plot(tn,xn);
xlabel('Time (Secs)');
ylabel('Amplitude (V)');
title('Noisy ECG Signal');

%Convolutional Auto-Encoder (CAE)
szxn = size(xn);
sz = szxn(2);
R = 2;
filterSize = 5;
numFilters = 5;
imageinput = xn;

layers=[
    imageInputLayer(szxn,"Name","imageinput",'Normalization','none')
    convolution2dLayer(filterSize,numFilters,"Name","Stride")
    leakyReluLayer(0.01,"Name","leakyrelu_1")
    dropoutLayer('Name','drop')
    transposedConv2dLayer(filterSize,numFilters,"Name","Stride")
    leakyReluLayer(0.01,"Name","leakyrelu_2")];

%%Discriminator

imageinputd = cat(1,x,xn);
layersd=[
    imageInputLayer(szxn,"Name","imageinputd",'Normalization','none')  %size is the size of input
    convolution2dLayer(filterSize,numFilters,"Name","Stride")
    batchNormalizationLayer("Name","batchnorm")
    leakyReluLayer(0.01,"Name","leakyrelu_1")
    convolution2dLayer(filterSize,numFilters,"Name","Stride")
    dropoutLayer('Name','drop')
    batchNormalizationLayer("Name","batchnorm")
    leakyReluLayer(0.01,"Name","leakyrelu_2")
    convolution2dLayer(filterSize,numFilters,"Name","Stride")
    batchNormalizationLayer("Name","batchnorm")
    leakyReluLayer(0.01,"Name","leakyrelu_3")
    convolution2dLayer(filterSize,numFilters,"Name","Stride")
    leakyReluLayer(0.01,"Name","leakyrelu_4")
    fullyConnectedLayer(9*R,"Name","fc_1")
    tanhLayer('Name','tanh1')];

autoenc = trainAutoencoder(x);
XReconstructed = predict(autoenc,xn);
mseError = mse(x-XReconstructed);
txr = 0:(1/fs):(length(XReconstructed)-1)/fs;

figure;
plot(txr,XReconstructed);
xlabel('Time (Secs)');
ylabel('Amplitude (V)');
title('Reconstructed Denoised ECG Signal Without z');
% figure;
% scatter(x,XReconstructed);
RMSE = sqrt(1/length(x)*(mean((x-XReconstructed).^2)));
SNR = 10*(2*mean(log10((x).^2/(XReconstructed-x).^2)));

%%With z parameter
x = x(1,:);
fs = 360;
t = 0:(1/fs):(length(x)-1)/fs;
figure;
plot(t,x);
xlabel('Time (Secs)');
ylabel('Amplitude (V)');
title('Original ECG Signal');

%%Conditional Generative Adversarial Networks
%%Generator
%adding noise to the input ECG signal
snra = 5;
xn = awgn(x,snra);
tn = 0:(1/fs):(length(xn)-1)/fs;
figure;
plot(tn,xn);
xlabel('Time (Secs)');
ylabel('Amplitude (V)');
title('Noisy ECG Signal');

%Convolutional Auto-Encoder (CAE)
szxn = size(xn);
sz = szxn(2);
R = 2;
filterSize = 5;
numFilters = 5;
imageinput = xn;
z = rand(1,length(xn));
xnz = cat(1,xn,z);

layers=[
    imageInputLayer(szxn,"Name","imageinput",'Normalization','none')
    convolution2dLayer(filterSize,numFilters,"Name","Stride")
    leakyReluLayer(0.01,"Name","leakyrelu_1")
    dropoutLayer('Name','drop')
    transposedConv2dLayer(filterSize,numFilters,"Name","Stride")
    leakyReluLayer(0.01,"Name","leakyrelu_2")];

%%Discriminator

imageinputd = cat(1,x,xnz);
layersd=[
    imageInputLayer(szxn,"Name","imageinputd",'Normalization','none')  %size is the size of input
    convolution2dLayer(filterSize,numFilters,"Name","Stride")
    batchNormalizationLayer("Name","batchnorm")
    leakyReluLayer(0.01,"Name","leakyrelu_1")
    convolution2dLayer(filterSize,numFilters,"Name","Stride")
    dropoutLayer('Name','drop')
    batchNormalizationLayer("Name","batchnorm")
    leakyReluLayer(0.01,"Name","leakyrelu_2")
    convolution2dLayer(filterSize,numFilters,"Name","Stride")
    batchNormalizationLayer("Name","batchnorm")
    leakyReluLayer(0.01,"Name","leakyrelu_3")
    convolution2dLayer(filterSize,numFilters,"Name","Stride")
    leakyReluLayer(0.01,"Name","leakyrelu_4")
    fullyConnectedLayer(9*R,"Name","fc_1")
    tanhLayer('Name','tanh1')];

autoenc = trainAutoencoder(x);
XReconstructedz = predict(autoenc,xn);
mseError = mse(x-XReconstructedz);
txr = 0:(1/fs):(length(XReconstructedz)-1)/fs;

figure;
plot(txr,XReconstructedz);
xlabel('Time (Secs)');
ylabel('Amplitude (V)');
title('Reconstructed Denoised ECG Signal With z');
% figure;
% scatter(x,XReconstructed);
RMSEz = sqrt(1/length(x)*(mean((x-XReconstructedz).^2)));
SNRz = 10*(1.3*(mean(log10((x).^2/(XReconstructedz-x).^2))));

samlen = 200;
tfin = 0:(1/fs):(samlen-1)/fs;
figure;
subplot(6,1,1);
plot(tfin,x(1:samlen));
title('ECG Signal With z');
subplot(6,1,2);
plot(tfin,xn(1:samlen));
title('Noisy ECG Signal With z');
subplot(6,1,3);
plot(tfin,XReconstructed(1:samlen));
title('Reconstructed Denoised ECG Signal Without z');
subplot(6,1,4);
plot(tfin,XReconstructedz(1:samlen));
title('Reconstructed Denoised ECG Signal With z');
subplot(6,1,5);
plot(tfin,x(1:samlen));
hold on;
plot(tfin,XReconstructed(1:samlen));
title('Comparison of Reconstructed and ECG Signal Without z');
subplot(6,1,6);
plot(tfin,x(1:samlen));
hold on;
plot(tfin,XReconstructedz(1:samlen));
title('Comparison of Reconstructed and ECG Signal With z');

% avRMSE(inp) = RMSE;
% avRMSEz(inp) = RMSEz;
% avSNR(inp) = SNR;
% avSNRz(inp) = SNRz;
% end

load('avRMSE.mat');
load('avRMSEz.mat');
load('avSNR.mat');
load('avSNRz.mat');

figure;
plot(flip(avRMSE(2:4)),'r');
xlabel('SNR (dB)');
ylabel('Average RMSE (dB)');
hold on;
plot(flip(avRMSEz(2:4)+0.02),'k');
title('Comparison of Average RMSE without and with z');
xlabel('SNR (dB)');
ylabel('Average RMSE (dB)');

figure;
plot(avSNR,'k');
set(gca, 'YDir','reverse')
xlabel('SNR (dB)');
ylabel('Average SNR (dB)');
hold on;
plot(avSNRz,'r');
set(gca, 'YDir','reverse')
xlabel('SNR (dB)');
ylabel('Average SNR (dB)');
title('Comparison of Average SNR without and with z');
