% for inp = 1:10
%     inp
%     
%%selection of signal
[file,path] = uigetfile('*.*','Select Signal');
Name1 = [path,file];
x = load(Name1);
x = struct2cell(x);
x = cell2mat(x(1,:));
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

%Pre-Processing of EEG Signal
n = 2;Rs = 60;Ws = [0.5,0.6];
[b,a] = cheby2(n,Rs,Ws,'bandpass');
figure;
freqz(b,a,[],fs1)

subplot(2,1,1)
ylim([-100 20])

%%Different Types of Feature Extraction Processes
%Principal Composition Analysis (PCA)
coeff = pca(x); %Extracting PCA coefficients from the input signal
coeff = coeff(1:1000,1);
% coeffinp(1:1000,inp) = coeff(1:1000,1);
% end

%Extracting Frequency statistics
ifq1 = instfreq(x,fs1); %extracting freqs
ifq2 = instfreq(x,fs2); %extracting freqs
ifq3 = instfreq(x,fs3); %extracting freqs
ifq = cat(1,ifq1,ifq2,ifq3);
szifq = size(ifq);
% ifqi(inp) = ifq(1);

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
alphai = alpha(1);
betai = beta(1);
bthetai = btheta(1);

%band ratios of EEG signal
betalp = beta(1)/alpha(3);
betthe = beta(1)/btheta(1);
alpthe = alpha(1)/btheta(1);
% end

%%Classification
%Training
load('labelbr.mat');
load('bndenrg.mat');
load('bndrat.mat');
load('bndfre.mat');
load('coeffinp.mat');

Mdlbr = fitcsvm(bndrat',labelbr);
Mdlfq = fitcsvm(bndfre',labelbr);
Mdlbe = fitcsvm(bndenrg',labelbr);
Mdlpc = fitcsvm(coeffinp',labelbr);

%Testing
test = input('Enter any number between 1 and 10: ');
tstbr = predict(Mdlbr,bndrat(1:3,test)')
tstfq = predict(Mdlfq,bndfre(1,test)')
tstbe = predict(Mdlbe,bndenrg(1:3,test)')
tstpc = predict(Mdlpc,coeff')

%accuracy
accuracybr = accur(Mdlbr,bndrat',labelbr,0.5);
accuracyfq = accur(Mdlfq,bndfre',labelbr,0.5);
accuracybe = accur(Mdlbe,bndenrg',labelbr,0.5);
accuracypc = accur(Mdlpc,coeffinp',labelbr,0.9);
