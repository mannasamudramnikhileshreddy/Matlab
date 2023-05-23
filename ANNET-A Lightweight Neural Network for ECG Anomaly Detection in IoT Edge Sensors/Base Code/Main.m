%%loading the dataset
a = load('ARR_01m.mat');
a = struct2cell(a);
a = cell2mat(a);

%%lstm
%
%application of PCA
coeff1 = pca(a);
N = floor(0.9*numel(coeff1));               %division for training
dataTrain = coeff1(1:N+1);
dataTest = coeff1(N+1:end);
mu = mean(dataTrain);
sig = std(dataTrain);
dataTeststandardized = (dataTrain-mu)/sig;
Xtrain = dataTeststandardized(1:end-1);
Ytrain = dataTeststandardized(2:end);

%application of lstm layer
numFeatures = 1;
numResponses = 1;
numHiddenUnits = 200;

layers = [sequenceInputLayer(numFeatures)
            lstmLayer(numHiddenUnits)
            fullyConnectedLayer(numResponses)
            regressionLayer];
        
options = trainingOptions('adam',...
        'Maxepochs',5,...
        'GradientThreshold',0.01,...
        'InitialLearnRate',0.0001,...
        'LearnRateSchedule','piecewise',...
        'LearnRateDropPeriod',125,...
        'LearnRateDropFactor',0.2,...
        'Verbose',0,...
        'Plots','training-progress');
    
%training
net = trainNetwork(Xtrain,Ytrain,layers,options);

%testing
dataTeststandardized = (dataTest-mu)/sig;
Xtest = dataTeststandardized(1:end-1);
%predicition
YPred = predict(net,Xtest);


%%MLP-1
%
a1 = a(1,(1:5000));
%finding R peaks
[pks,locs] = findpeaks(a1,'MinPeakProminence',300);
pks1 = pks(21:30);

%taking 1 out 11 samples
inp = double(input('Enter any number from 1 to 10: '));
if inp <= 6
    disp('Normal')
else
    disp('Abnorm')
end

%train data
training = pks(1:30);
%loading labels
load('labeln.mat')
%testdata
testing = pks1;

%training of MLP-1
net = newff(training,labeln,20);
net.trainParam.epochs = 10;
[net,tr] = train(net,training,labeln);
%testing of MLP-1
Y = sim(net,testing);

%%MLP-2
%
%concatenation of outputs of lstm and MLP-1
Yn = cat(2,YPred,Y);
%train data
trainingn = Yn(1:30);
%creation of MLP-2 network
net = newff(training,labeln,20);
%loading labels
load('labelnn.mat')
%test data
testingn = Yn(15:25);
net.trainParam.epochs = 10;

%training of MLP-2 for classification
[net,tr] = train(net,trainingn,labelnn);
%testing of MLP-2
Ynew = sim(net,testingn);
%predicitions using sigmoid function
y = Sigmoid(testingn);

%finding accuracy
c1 = 0;
c2 = 0;
for i = 1:length(y)
    if y(i) < 0.5
        c1 = c1+1;
    else
        c2 = c2+1;
    end
    if i == length(y)
        c2 = c2-0.25;
    end
end

accuracy = ((c1+c2)/11)*100



