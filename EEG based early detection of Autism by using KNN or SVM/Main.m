
% for in = 1:10
%%Selection of Signal
[file,path] = uigetfile('*.*','Select Signal');
Name1 = [path,file];
x = load(Name1);
x = struct2cell(x);
x = cell2mat(x);
x = x(1,:);
% x = table2array(x{1,1});
fs1 = 20;   %sampling frequency1
fs2 = 50;    %sampling frequency2
fs3 = 30;    %sampling frequency3
t = 0:(1/fs1):(length(x)-1)/fs1;
figure;
plot(t,x);
xlabel('Time in secs');
ylabel('Amplitude in mV');
title('Original EEG Signal');


%%Pre-processing
%Feature Extraction
%Extracting Frequency statistics
ifq1 = instfreq(x,fs1)*10; %extracting freqs
ifq2 = instfreq(x,fs2)*10; %extracting freqs
ifq3 = instfreq(x,fs3)*10; %extracting freqs
ifq = cat(1,ifq1,ifq2,ifq3);
szifq = size(ifq);

%characteristics extraction using band energies
%extracting alpha, beta and theta
iiffa = 1;iiffb = 1;iiffc = 1;
alpha = 0;beta = 0;btheta = 0;
for iiff1 = 1:szifq(1)
    for iiff2 = 1:szifq(2)
    if ifq(iiff1,iiff2) >= 8 && ifq(iiff1,iiff2) < 13
        alpha(iiffa) = ifq(iiff1,iiff2);
        iiffa = iiffa+1;
    elseif ifq(iiff1,iiff2) >= 13 && ifq(iiff1,iiff2) < 30
        beta(iiffb) = ifq(iiff1,iiff2);
        iiffb = iiffb+1;
    elseif ifq(iiff1,iiff2) >= 4 && ifq(iiff1,iiff2) < 8
        btheta(iiffc) = ifq(iiff1,iiff2);
        iiffc = iiffc+1;
    end
    if iiff1 == szifq(1)
        if length(beta) == 1
            beta = alpha(1:3);
        end
        if length(btheta) > 1 && length(btheta) < 3
            btheta(2,:) = alpha(1:3);
        end
    end
    end
end

alphatst = mean(alpha);
betatst = mean(beta);
bthetatst = mean(btheta);
tstfet = cat(2,alphatst,betatst,bthetatst);
% alphain(in) = mean(alpha);
% betain(in) = mean(beta);
% bthetain(in) = mean(btheta);
% if in == 10
% features = cat(1,alphain,betain,bthetain);
% end
% end

%%Classification
load('features.mat')
load('labels.mat')

%Separating the dataset into training and testing with labels

%Training
load('trainfet.mat')
load('trainlabels.mat')

datinp = trainfet;
labinp = trainlabels;

Mdl = fitcecoc(datinp,labinp);

%Testing
ypred = predict(Mdl,tstfet);  %prediction

%Predicted output
if ypred == '0'
    out = 'Autism'
elseif ypred == '1'
    out = 'NotAut'
end

%accuracy
load('testfet.mat')
load('testlabels.mat')

accuracy = accur(labels,features,Mdl);
