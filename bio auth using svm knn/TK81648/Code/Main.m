% for in = 1:10
[file,path] = uigetfile('*.*','Select Signal');
Name1 = [path,file];
x1 = load(Name1);
x1 = struct2cell(x1);
x1 = cell2mat(x1);
x1 = x1(1,:);
fs = 1000;
t1 = 0:(1/fs):(length(x1)-1)/fs;
figure;
plot(t1,x1);
xlabel('Time in Secs');
ylabel('Amplitude in mV');
title('Original ECG signal of sample');

%%Application of EMD
%
[imf,residue] = emd(x1,'MaxNumIMF',8);
% features(in,:) = imf(1:1000,1)';
% end
te1 = 0:(1/fs):(length(imf)-1)/fs;
figure;
subplot(4,1,1);
plot(te1,imf(:,1));
xlabel('Time in Secs');
ylabel('Imf 1');
title('Imfs Extracted from sample');
subplot(4,1,2);
plot(te1,imf(:,2));
xlabel('Time in Secs');
ylabel('Imf 2');
subplot(4,1,3);
plot(te1,imf(:,3));
xlabel('Time in Secs');
ylabel('Imf 3');
subplot(4,1,4);
plot(te1,imf(:,4));
xlabel('Time in Secs');
ylabel('Imf 4');

figure;
subplot(2,1,1);
plot(te1,imf(:,5));
xlabel('Time in Secs');
ylabel('Imf 5');
title('Imfs Extracted from sample');
subplot(2,1,2);
plot(te1,imf(:,6));
xlabel('Time in Secs');
ylabel('Imf 6');

%%Feature Extraction
%
%shannon energy
SE = imf(:,1).*log(abs(imf(:,1)));
%skewness
sk = skewness(imf(:,1));
%variance
mimf = mean(imf(:,1));
N = length(imf(:,1));
for i = 1:N
    V(i) = ((imf(i,1)-mimf)^2)/N;
end
%occupied band width
bw = obw(imf(:,1));
%median frequency
[pks,locs] = findpeaks(imf(:,1));

for ii = 1:length(locs)
    freq(ii) = medfreq(imf(locs(ii)-1:locs(ii)+1,1))*50;
end

figure;
plot(V(1:700)/100,freq(1:700),'ro');
xlabel('Variance');
ylabel('Median Frequency');
title('Frequency Distribution');

%%classification
%
load('labels.mat')
load('features.mat')

%%for selection of signal
%
ssimf = imf(1:1000,1)';
for ss = 1:10
    if features(ss,:) == ssimf
        ssin = ss;
    end
end

%classification using SVM
%Linear SVM
Mdl = fitcsvm(features,labels);
Y = predict(Mdl,features(ssin,:))

%Multi-Class SVM or Quadratic SVM
Mdlm = fitcecoc(features,labels);
Ym = predict(Mdlm,features(ssin,:))

%Binary SVM
load('labelsb.mat')
Mdlr = fitrsvm(features,labelsb);
Yr = predict(Mdlr,features(ssin,:))

if Yr <= 0.1
    Yrout = 'NotEnrolled'
else
    Yrout = 'Enrolled'
end

%Gaussian SVM
load('labelsq.mat');
Mdlg = fitrlinear(features,labelsq);
Yg = predict(Mdlg,features(ssin,:))

if Yg <= 0.1
    Ygout = 'NotEnrolled';
else
    Ygout = 'Enrolled';
end

%classification using KNN
%KNN with weighted
Mdlkw = fitcknn(features,labels);
Ykw = predict(Mdlkw,features(ssin,:))

%KNN with 5 Nearest Neighbors
Mdlk5 = fitcknn(features,labels,'NumNeighbors',5);
Yk5 = predict(Mdlk5,features(ssin,:))

%KNN with 9 Nearest Neighbors
Mdlk9 = fitcknn(features,labels,'NumNeighbors',9);
Yk9 = predict(Mdlk9,features(ssin,:))

%classification using Decision Tree
tree = fitctree(features,labels);
Ytree = predict(tree,features(ssin,:))

accuracyL = accu(Mdl,features,labels,1);
accuracyQ = accu(Mdlm,features,labels,2);
accuracyB = accu(Mdlr,features,labelsb,3);
accuracyG = accu(Mdlg,features,labelsq,4);
accuracyK = accu(Mdlkw,features,labels,5);
accuracyK9 = accu(Mdlkw,features,labels,6);
accuracyK5 = accu(Mdlkw,features,labels,7);
accuracyDT = accu(tree,features,labels,8);