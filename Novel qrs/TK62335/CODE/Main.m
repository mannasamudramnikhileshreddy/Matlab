clc
clear all
close all

% Read ECG signal
% input signal data to 'val' variable
% converting a raw ECG signal into normal ECG signal

    load('100m.mat');
    val = (val - 0)/200;                % removing "base" and "gain"
    ECG = val(1,1:2000);                % choosing Lead 1 (V1) data and 2000 datapoints 
    fs = 400;                           % sampling frequecy
    t = (0:length(ECG)-1)/fs;           % time
 
% Plot signal before processing
    figure(1)
    plot(t,ECG);title('ECG Signal');
    xlabel("Time (s)")
    ylabel("ECG Amplitude (mV)")
  
 
% 1st denoising with median filter 100 ms sliding window
    s_win1 = 0:1/fs:0.01; % 1st sliding window = 100ms
    mf1 = medfilt1(ECG, length(s_win1)-1); % apply median filter
    
 
% 2nd denoising with median filter (200ms sliding window)
    s_win2 = 0:1/fs:0.02; %sliding window 2 = 200 ms
    mf2 = medfilt1(ECG, length(s_win2)-1); % apply median filter
    figure(3);
    plot(t, mf2);title('Denoised ECG with median filter of 200ms sliding window ');
    xlabel("Time (s)")
    ylabel("ECG Amplitude (mV)")
    
    %% Initialize
qrs_c = zeros(1,0); %amplitude of R % TMW add, was qrs_c = [];
qrs_i = zeros(1,0); %amplitude of R % TMW add, was qrs_i = [];
coder.varsize('qrs_c','qrs_i'); % TMW add
SIG_LEV = 0; 
nois_c =[];
nois_i =[];
delay = 0;
skip = 0; % becomes one when a T wave is detected
not_nois = 0; % it is not noise when not_nois = 1
selected_RR =[]; % Selected RR intervals
m_selected_RR = 0;
mean_RR = 0;
qrs_i_raw =[];
qrs_amp_raw=[];
ser_back = 0; 
test_m = 0;
SIGL_buf = [];
NOISL_buf = [];
THRS_buf = [];
SIGL_buf1 = [];
NOISL_buf1 = [];
THRS_buf1 = [];
 
%% Noise cancelation(Filtering) % Filters (Filter in between 5-15 Hz)
if fs == 200
%% Low Pass Filter  H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2
b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
a = [1 -2 1];
h_l = filter(b,a,[1 zeros(1,12)]); 
ecg_l = conv (mf2 ,h_l);
ecg_l = ecg_l/ max( abs(ecg_l));
delay = 9;
 
%% High Pass filter H(z) = (-1+32z^(-16)+z^(-32))/(1+z^(-1))
b = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a = [1 -1];
h_h = filter(b,a,[1 zeros(1,32)]); 
ecg_h = conv (ecg_l ,h_h);
ecg_h = ecg_h/ max( abs(ecg_h));
delay = delay + 11; 
 
else
    %% bandpass filter for Noise cancelation of other sampling frequencies(Filtering)
    fs = 400;  % TMW add: butter's inputs must be constants
    f1=5; %cuttoff low frequency to get rid of baseline wander
    f2=15; %cuttoff frequency to discard high frequency noise
    Wn=[f1 f2]*2/fs; % cutt off based on fs
    N = 3; % order of 3 less processing
    [a,b] = butter(N,Wn); %bandpass filtering
    ecg_h = filtfilt(a,b,mf2);
    ecg_h = ecg_h/ max( abs(ecg_h));
 
end

%% derivative filter H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2))
h_d = [-1 -2 0 2 1]*(1/8);%1/8*fs
ecg_d = conv (ecg_h ,h_d);
ecg_d = ecg_d/max(ecg_d);
delay = delay + 2;
 
%% Squaring nonlinearly enhance the dominant peaks
ecg_s = ecg_d.^2;
 
%% Moving average Y(nt) = (1/N)[x(nT-(N - 1)T)+ x(nT - (N - 2)T)+...+x(nT)]
ecg_m = conv(ecg_s ,ones(1 ,round(0.150*fs))/round(0.150*fs));
delay = delay + 15;
figure;
plot(ecg_m);title('moving average filtered output');

figure
findpeaks(ecg_m,'MinPeakProminence',0.1)
 
figure
 
% Show all peaks in the original ecg signal
findpeaks(ECG,'MinPeakProminence',0.5)
xlabel('Samples')
ylabel('Amplitude')
title('Detecting Peaks')

