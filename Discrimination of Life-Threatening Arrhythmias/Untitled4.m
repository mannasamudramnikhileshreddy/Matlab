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
%load('Dset.mat')
%for in = 1:20
    %if val(1,(1:100)) == Dset(in,(1:100))
        %in;
        %break;
    %end
%end
%in;

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

%%Harmonic phase distribution
r = thd(y1(:,1));
%r1 = thd(real(af1));
%r2 = thd(real(af2));