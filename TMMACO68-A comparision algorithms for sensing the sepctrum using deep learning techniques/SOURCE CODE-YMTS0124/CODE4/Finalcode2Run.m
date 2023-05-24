clc;
clear all;
close all;
warning off all
parameter=64;   
fftlength=64; 
a1=input('Enter the data to be selected 1 for nosie 2 for signal with noise: ');
nloop=1000;  
level=8;                                       
for i=1:nloop    
           Eb_db = -20;
           Eb = 10.^(Eb_db./10);
           Pf = 0:0.05:1; 
           pd=zeros(1,length(Pf));
           pdd=zeros(1,length(Pf));
            for m = 1:length(Pf)
           input_dat =randn(1,parameter);
           magin = sqrt(Eb).*randn(1,parameter); 
           x=3.*sin(input_dat);
           tx_sig=x.*input_dat;
           tx_sig = abs(round(magin + tx_sig));
           H0= magin + input_dat; 
           entropynoi(m) = wentropy(magin,'shannon') ;
           entropy(m) = wentropy(tx_sig,'shannon') ;
            end
            if a1==1
                a2=entropynoi;
            elseif a1==2
                a2=entropy;
            end
end
    load('spectrumdata.mat')  
    XTrain =Xtrain;
    YTrain=categorical(Ytrain);
    inputSize = 1;
    numHiddenUnits = 100;
    numClasses = 2;
    layers = [ ...
        sequenceInputLayer(inputSize)
        lstmLayer(numHiddenUnits,'OutputMode','last')
        reluLayer       
        fullyConnectedLayer(numClasses)
        softmaxLayer
        classificationLayer];
    maxEpochs = 100;
    miniBatchSize = 200;
    options = trainingOptions('sgdm', ...
          'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'GradientThreshold',4, ...
         'LearnRateSchedule', 'piecewise', ...
          'LearnRateDropFactor', 0.2, ...
          'LearnRateDropPeriod', 5, ...
        'Verbose',false, ...
        'Plots','training-progress');
    net = trainNetwork(XTrain,YTrain,layers,options);
    YPred = classify(net,a2,'MiniBatchSize',miniBatchSize);
    Spectrumstatu=double(YPred);
    if Spectrumstatu==1
        msgbox('Primary User is acessing the spectrum');
    elseif Spectrumstatu==2
        msgbox('Primary User is not acessing the spectrum');
    end        
        figure,
        bar(magin*10);
        hold on
        axis([0 20 0 2]);
        xlabel('AWGN noise in amplitude');
        ylabel('Ni-number of iterations of i-th bin');
        title('Histogram of data with user absence');
        hold off    
        figure,bar(H0*10);
        hold on
        axis([0 20 0 5]);
        xlabel('Signal+AWGN noise in amplitude');
        ylabel('Ni-number of iterations of i-th bin');
        title('Histogram of data with user presence');
        hold off   
    