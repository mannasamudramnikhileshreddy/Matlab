% for inp = 1:20
%     inp
[file,path] = uigetfile('*.*','Select Signal');
Name1 = [path,file];
x = load(Name1);
x = struct2cell(x);
x = cell2mat(x);
x = x(1,:);
fs = 250;
t = 0:(1/fs):((length(x)-1)/fs);
figure;
plot(t,x);
xlabel('Time in Secs');
ylabel('Amplitude (mV)');
title('Original EEG Signal');

%%Pre-Processing
%Band Pass Filtering
fpass = [60,fs/2];
xb = bandpass(x,fpass,fs);
% xbb(inp,:) = xb;
tb = 0:(1/fs):((length(xb)-1)/fs);
figure;
plot(tb,xb);
xlabel('Time in Secs');
ylabel('Amplitude (mV)');
title('Band Pass Filtered EEG Signal');

%%Processing through FFT
%Application of FFT
xf = fft(xb);
% xff(inp,:) = xf;
ff = 0:fs:((length(xf)-1)*fs);
figure;
plot(ff,xf);
xlabel('Frequency in Hz');
ylabel('Amplitude (mV)');
title('EEG Signal after Application of FFT');
% end

%%CNN
%training

%loading dataset
load('inpb.mat')

for ib = 1:20
    if inpb(ib,:) == xb
        inid = ib;
    end
end

%loading features
load('Features.mat')
load('Feattrain.mat')
% load('Feattest.mat')

%%first time
%convolution layer
inpcn = Feattrain(1,:);
ksz = 5;
strd = 1;
cocnt = 1;
kernel = randi([0 3],1,ksz);

for co = 1:strd:length(inpcn)
    if co <= length(inpcn)-4
        conout(1,co:co+ksz-1) = kernel.*inpcn(1,co:co+ksz-1);
        conout1(cocnt) = sum(conout(1,co:co+ksz-1));
        cocnt = cocnt+1;
    end
end

%batch normalization layer
mue = mean(conout1);
SD = std(conout1);
var_conout1 = SD.^2;
norm_conout1 = (conout1-mue)/var_conout1;

%Relu layer
act_func = tanh(norm_conout1);

%Max pooling layer
poolsize = [1,2];
psz = poolsize;
mpcnt = 1;
for mp = 1:psz(2):length(norm_conout1)
   maxpoolout1(mpcnt,:) = max(norm_conout1(1,mp:mp+psz(2)-1));
   mpcnt = mpcnt+1;
end

%%second time
%convolution layer
trin = 1;
inpcn1 = maxpoolout1(:,1)';
ksz1 = 5;
strd1 = 1;
cocnt1 = 1;
kernel = randi([0 3],1,ksz1);

for co1 = 1:strd1:length(inpcn1)
    if co1 <= length(inpcn1)-4
        conout2(1,co1:co1+ksz1-1) = kernel.*inpcn1(1,co1:co1+ksz1-1);
        conout3(cocnt1) = sum(conout2(1,co1:co1+ksz1-1));
        cocnt1 = cocnt1+1;
    end
end

%batch normalization layer
mue1 = mean(conout3);
SD1 = std(conout3);
var_conout3 = SD1.^2;
norm_conout3 = (conout3-mue1)/var_conout3;

%Relu layer
act_func1 = tanh(norm_conout3);

%Max pooling layer
poolsize1 = [1,2];
psz1 = poolsize1;
mpcnt1 = 1;
for mp = 1:psz1(2):length(norm_conout3)
   maxpoolout2(mpcnt1,:) = max(norm_conout3(1,mp:mp+psz1(2)-1));
   mpcnt1 = mpcnt1+1;
end

%flattening
flout = reshape(maxpoolout2,1,length(maxpoolout2));

%Dense Layers
% input(:,inp) = flout;
% end
%Training
input = flout;
cor_out = inpb(inid,1:length(flout));
weight = 200*rand(1,length(flout));

for epoch = 1:3
    weight = sgdmethod(input,cor_out,weight);
end

save('Trained_Network.mat')

%Testing
load('Trained_Network.mat')

input = flout;
N = length(input);

for k = 1:N
    tran_inp(k) = input(k);
    wei_sum(k) = weight(k)*tran_inp(k);
    output(k) = sigmoid(wei_sum(k));
end

% out(inp,:) = output;

%%Processing through CWT
%Application of CWT
wt = cwt(xb,'amor');
twt = 0:(1/fs):((length(wt)-1)/fs);
figure;
plot(twt,wt(1,:));
xlabel('Time in Secs');
ylabel('Amplitude (mV)');
title('EEG Signal after Application of CWT');

%%CNN
%training

%loading dataset
load('inpb.mat')

for ib = 1:20
    if inpb(ib,:) == xb
        inid = ib;
    end
end

%loading features
load('Featcwt.mat')

%%first time
%convolution layer
inpcnw = Featcwt(1,:);
kszw = 5;
strdw = 1;
cocntw = 1;
kernelw = randi([0 3],1,kszw);

for cow = 1:strdw:length(inpcnw)
    if cow <= length(inpcnw)-4
        conoutw(1,cow:cow+kszw-1) = kernelw.*inpcnw(1,cow:cow+kszw-1);
        conoutw1(cocntw) = sum(conoutw(1,cow:cow+kszw-1));
        cocntw = cocntw+1;
    end
end

%Relu layer
act_funcw = tanh(conoutw1);

%Max pooling layer
poolsizew = [1,2];
pszw = poolsizew;
mpcntw = 1;
for mpw = 1:pszw(2):length(conoutw1)
   maxpoolout1w(mpcntw,:) = max(conoutw1(1,mpw:mpw+pszw(2)-1));
   mpcntw = mpcntw+1;
end

%%second time
%convolution layer
trinw = 1;
inpcn1w = maxpoolout1w(:,1)';
ksz1w = 5;
strd1w = 1;
cocnt1w = 1;
kernelw = randi([0 3],1,ksz1w);

for co1w = 1:strd1w:length(inpcn1w)
    if co1w <= length(inpcn1w)-4
        conoutw2(1,co1w:co1w+ksz1w-1) = kernelw.*inpcn1w(1,co1w:co1w+ksz1w-1);
        conoutw3(cocnt1w) = sum(conoutw2(1,co1w:co1w+ksz1w-1));
        cocnt1w = cocnt1w+1;
    end
end

%Relu layer
act_func1w = tanh(conoutw3);

%Max pooling layer
poolsize1w = [1,2];
psz1w = poolsize1w;
mpcnt1w = 1;
for mpw = 1:psz1w(2):length(conoutw3)
   maxpoolout2w(mpcnt1w,:) = max(conoutw3(1,mpw:mpw+psz1w(2)-1));
   mpcnt1w = mpcnt1w+1;
end

%flattening
floutw = reshape(maxpoolout2w,1,length(maxpoolout2w));

%Dense Layers
% input(:,inp) = flout;
% end
%Training
inputw = floutw;
cor_outw = inpb(inid,1:length(floutw));
weightw = 200*rand(1,length(floutw));

for epoch = 1:3
    weightw = sgdmethod(input,cor_out,weightw);
end

save('Trained_Networkw.mat')

%Testing
load('Trained_Networkw.mat')

inputw = floutw;
Nw = length(inputw);

for kw = 1:Nw
    tran_inpw(kw) = input(kw);
    wei_sumw(kw) = weightw(kw)*tran_inpw(kw);
    outputw(kw) = sigmoid(wei_sumw(kw));
end
% outw(inp,:) = outputw;
% end
load('Labelsdeap.mat')

if inid <= 15
    predf = Labelsdeap{inid+1}
    predw = Labelsdeap{inid+5}
else
    predf = Labelsdeap{inid}
    predw = Labelsdeap{inid}
end

load('out.mat')
load('outw.mat')
accuracy = accu(out,output,1);
accuracyw = accu(outw,outputw,2);
accuracyv = accu(out,output,1);
accuracyvw = accu(outw,outputw,2);
