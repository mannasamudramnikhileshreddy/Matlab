a = load('100m.mat');
a = struct2cell(a);
a = cell2mat(a);
a1 = a(1,:);
a2 = a(2,:);
%doublea=double(a)
%a = cast(a(1,:),'single')
Ts = -1;
A = [1.1269   -0.4940    0.1129 
     1.0000         0         0 
          0    1.0000         0];
B = [-0.3832
      0.5919

      0.5191];

C = [1 0 0];

D = 0;
sys = ss(A,[B B],C,D,Ts,'InputName',{'a1' 'a2'},'OutputName','y');
t = 1:length(a);
plot(t,a);
xlabel('Time in Secs');
ylabel('Amplitude');
title('Input ECG Signal');
Q=2;
%doubleQ=double(Q)
R=1;
%doubleR=double(R)
b=kalman(sys,Q,R);
%N = 10;
%Fc = 0.4;
%[c,d] = butter(N,Fc);
iir = dsp.IIRFilter('Numerator',b.C(:,1)',...
    'Denominator',1);
[pks,locs] = findpeaks(a1);
%for i = 1:length(a1)
segment = a1(locs(9):locs(26));
segment2= a1(locs(11):locs(66));
segment3=a1(locs(1):locs(34));
ts = 1:length(segment);
figure
plot(ts,segment)
xlabel('Time in Secs');
ylabel('Amplitude');
title('Segmented ECG Signal Segment1');

ts2 = 1:length(segment2);
figure
plot(ts2,segment2)
xlabel('Time in Secs');
ylabel('Amplitude');
title('Segmented ECG Signal Segment2');

ts3 = 1:length(segment3);
figure
plot(ts3,segment3);
xlabel('Time in Secs');
ylabel('Amplitude');
title('Segmented ECG Signal using fiducial features');

%name = 'segment3';
%f = dbwavf(wname)

xden = wden(segment3,'sqtwolog','s','mln',3,'sym4');
ts4=1:length(xden);
figure
plot(ts4,xden);
xlabel('Time in Secs');
ylabel('Amplitude');
title('Segmented ECG Signal using Symlet Wavelet');

xden1 = wden(segment3,'sqtwolog','s','mln',3,'db4');
ts5=1:length(xden1);
figure
plot(ts5,xden1);
xlabel('Time in Secs');
ylabel('Amplitude');
title('Segmented ECG Signal using Daubechies Wavelet');

x1=segment3(:,1);
x2=xden(:,1);
Euc_distance = norm(x1-x2);
DTW1=dtw(x1,x2);
for i=1:locs(:,1)
    z(i)= locs(i+1)-locs(i);
end    

load('Y.mat');
k=[-30.0690705951149,-31.2530726208355,-32.3960278832253,-33.3927983886654,-34.0692916816571;-30.3606028308014,-31.5711397866414,-33.0178298467484,-33.9303370940205,-34.2699126160638]
Mdl = fitcknn(k,Y);
out = predict(Mdl,k(2,:))
accuracy1 = accur(Mdl,k,Y,1)

%non fiducial symlet4
load('nffeats4.mat');
Mdl2 = fitcknn(nffeats4,Y);
out2 = predict(Mdl2,nffeats4(1,:))
accuracy2 = accur(Mdl2,nffeats4,Y,2)

%non fiducial daubechies4
load('nffeatd4.mat');
Mdl3 = fitcknn(nffeatd4,Y);
out3 = predict(Mdl3,nffeatd4(1,:))
accuracy3 = accur(Mdl3,nffeatd4,Y,3)

load('accuracyf.mat');
load('accuracynfs.mat');
load('accuracynfd.mat');
figure;
seg_wind = 0:0.1:0.6;
plot(seg_wind,accuracyf,'r-o');
hold on;
plot(seg_wind,accuracynfs,'g-o');
hold on;plot(seg_wind,accuracynfd,'k-o');
xlabel('Segmentation Window');
ylabel('Accuracy');
title('Accuracy Comparison');
legend('accuracyf','accuracynfd','accuracynfs','location','southeast');
