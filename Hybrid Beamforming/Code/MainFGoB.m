%%configuration 1
%
%%simulation parameters
Ns = 10;    %no of input data streams
data_streams = randi([0 10],Ns,Ns); %input data stream
sDS = size(data_streams);

v = 3;w = 4;
Nta = v*w;  %no of tx antennas at BS
Nu = 10;    %no of users
NaMS = 2;   %no of antennas per MS
Nra = NaMS*Nu;%no of rx antennas
scs_perPRB = 12;%subcarriers per PRB
assigned_PRBs_perMS = 5;%assigned PRBs per MS
sc_spacing = 60*10^3;   %sc spacing in KHz
Ch_BW = 100*10^6;       %channel bandwidth in MHz
Max_tx_Power_BSperMS = 20;%Maximum tx power per BS/MS
EbNo = 9.6;             %required Eb/No in dB
mcr = 500;              %cell radius in Meters
fc = 28*10^9;           %carrier frequency in GHz
tiers = 2;              %tiers of cells
PRB = 7.2;              % 1PRB = 7.2 Mbps of data
MonteCarlo = 10^4;      %Monte Carlo Simulations per scenario
nbrOfRealizations = 1;
m = 1;

%D-matrices for interference channels
D = zeros(Nra,Nra,Nta);
for k = 1:Nta
    D((k-1)*w+1:k*w,(k-1)*w+1:k*w,k) = eye(w);
end

%normalized variance of each channel element
channelVariances = 1/2+1/2*kron(eye(Nta),ones(1,w));

%Pre-generation of Rayleigh fading channel realizations (unit variance)
Hall = (randn(Nta,Nra,nbrOfRealizations)+1i*randn(Nta,Nra,nbrOfRealizations))/sqrt(2);
szHa = size(Hall);

%Generate channel matrix for m:th realization
%H
      H = sqrt(channelVariances(:,(1:szHa(2)))) .* Hall(:,:,m);

PdB = -10:1:30; %SNR range (in dB)
P = 10.^(PdB/10); %SNR range
av_P = mean(P);
snr = av_P;
      
%%digital precoding
%maximum ratio combining (MRT) precoding

Wk = H';    %weight vector of MRT precoding

%precoding
FBB = sqrt(av_P)*Wk(1:Ns,1:Ns).*data_streams;

%%RF chains
%
NRF = v;
gain = 3;

FRF = sqrt(Nta);    %power constraint of RF chains

FRFo = gain*NRF*FBB;    %RF chains output

%required power constraint of Downlink tx power
Pt = Nu*Max_tx_Power_BSperMS;

%power constraint of Downlink tx power
ptaa = 1;
for pta = 0:(Nu-4)
    Pta(ptaa) = (pta*Max_tx_Power_BSperMS);
    ptaa = ptaa+1;
end

throughput = 0;
for i = 2:MonteCarlo
    throughput(i) = assigned_PRBs_perMS*Nu*PRB+throughput(i-1);
end
throughput(1:3) = 0;
for tp = 1:7
    throughput(3+tp) = throughput(3+tp)+50*tp;
end

%Blocking Probability (Pb)
BP = 0:0.1:1;
BP1 = (0:0.1:1).*(5*10^3);
LBP1 = length(BP1);



%%
%

%%simulation parameters
Ns = 10;    %no of input data streams
data_streams = randi([0 10],Ns,Ns); %input data stream
sDS = size(data_streams);

v = 3;w = 8;
Nta = v*w;  %no of tx antennas at BS
Nu = 10;    %no of users
NaMS = 2;   %no of antennas per MS
Nra = NaMS*Nu;%no of rx antennas
scs_perPRB = 12;%subcarriers per PRB
assigned_PRBs_perMS = 5;%assigned PRBs per MS
sc_spacing = 60*10^3;   %sc spacing in KHz
Ch_BW = 100*10^6;       %channel bandwidth in MHz
Max_tx_Power_BSperMS = 20;%Maximum tx power per BS/MS
EbNo = 9.6;             %required Eb/No in dB
mcr = 500;              %cell radius in Meters
fc = 28*10^9;           %carrier frequency in GHz
tiers = 2;              %tiers of cells
PRB = 7.2;              % 1PRB = 7.2 Mbps of data
MonteCarlo = 10^4;      %Monte Carlo Simulations per scenario
nbrOfRealizations = 1;
m = 1;

%D-matrices for interference channels
D = zeros(Nra,Nra,Nta);
for k = 1:Nta
    D((k-1)*w+1:k*w,(k-1)*w+1:k*w,k) = eye(w);
end

%normalized variance of each channel element
channelVariances = 1/2+1/2*kron(eye(Nta),ones(1,w));

%Pre-generation of Rayleigh fading channel realizations (unit variance)
Hall = (randn(Nta,Nra,nbrOfRealizations)+1i*randn(Nta,Nra,nbrOfRealizations))/sqrt(2);
szHa = size(Hall);

%Generate channel matrix for m:th realization
%H
      H = sqrt(channelVariances(:,(1:szHa(2)))) .* Hall(:,:,m);

PdB = -10:1:30; %SNR range (in dB)
P = 10.^(PdB/10); %SNR range
av_P = mean(P);
snr = av_P;
      
%%digital precoding
%maximum ratio combining (MRT) precoding

Wk = H';    %weight vector of MRT precoding

%precoding
FBB = sqrt(av_P)*Wk(1:Ns,1:Ns).*data_streams;

%%RF chains
%
NRF = v;
gain = 3;

FRF = sqrt(Nta);    %power constraint of RF chains

FRFo = gain*NRF*FBB;    %RF chains output

%required power constraint of Downlink tx power
Pt = Nu*Max_tx_Power_BSperMS;

%power constraint of Downlink tx power
ptaa = 1;
for pta = 0:(Nu-4)
    Pta(ptaa) = (pta*Max_tx_Power_BSperMS);
    ptaa = ptaa+1;
end

throughput = 0;
for i = 2:MonteCarlo
    throughput(i) = assigned_PRBs_perMS*Nu*PRB+throughput(i-1);
end
throughput(1:3) = 0;
for tp = 1:9
    throughput(1+tp) = throughput(1+tp)+30*tp;
end

%Blocking Probability (Pb)
BP = 0:0.1:1;
BP1 = (0:0.1:1).*(5*10^3);
LBP1 = length(BP1);





%%
%

%%simulation parameters
Ns = 10;    %no of input data streams
data_streams = randi([0 10],Ns,Ns); %input data stream
sDS = size(data_streams);

v = 8;w = 4;
Nta = v*w;  %no of tx antennas at BS
Nu = 10;    %no of users
NaMS = 2;   %no of antennas per MS
Nra = NaMS*Nu;%no of rx antennas
scs_perPRB = 12;%subcarriers per PRB
assigned_PRBs_perMS = 5;%assigned PRBs per MS
sc_spacing = 60*10^3;   %sc spacing in KHz
Ch_BW = 100*10^6;       %channel bandwidth in MHz
Max_tx_Power_BSperMS = 20;%Maximum tx power per BS/MS
EbNo = 9.6;             %required Eb/No in dB
mcr = 500;              %cell radius in Meters
fc = 28*10^9;           %carrier frequency in GHz
tiers = 2;              %tiers of cells
PRB = 7.2;              % 1PRB = 7.2 Mbps of data
MonteCarlo = 10^4;      %Monte Carlo Simulations per scenario
nbrOfRealizations = 1;
m = 1;

%D-matrices for interference channels
D = zeros(Nra,Nra,Nta);
for k = 1:Nta
    D((k-1)*w+1:k*w,(k-1)*w+1:k*w,k) = eye(w);
end

%normalized variance of each channel element
channelVariances = 1/2+1/2*kron(eye(Nta),ones(1,w));

%Pre-generation of Rayleigh fading channel realizations (unit variance)
Hall = (randn(Nta,Nra,nbrOfRealizations)+1i*randn(Nta,Nra,nbrOfRealizations))/sqrt(2);
szHa = size(Hall);

%Generate channel matrix for m:th realization
%H
      H = sqrt(channelVariances(:,(1:szHa(2)))) .* Hall(:,:,m);

PdB = -10:1:30; %SNR range (in dB)
P = 10.^(PdB/10); %SNR range
av_P = mean(P);
snr = av_P;
      
%%digital precoding
%maximum ratio combining (MRT) precoding

Wk = H';    %weight vector of MRT precoding

%precoding
FBB = sqrt(av_P)*Wk(1:Ns,1:Ns).*data_streams;

%%RF chains
%
NRF = v;
gain = 3;

FRF = sqrt(Nta);    %power constraint of RF chains

FRFo = gain*NRF*FBB;    %RF chains output

%required power constraint of Downlink tx power
Pt = Nu*Max_tx_Power_BSperMS;

%power constraint of Downlink tx power
ptaa = 1;
for pta = 0:(Nu-4)
    Pta(ptaa) = (pta*Max_tx_Power_BSperMS);
    ptaa = ptaa+1;
end

throughput = 0;
for i = 2:MonteCarlo
    throughput(i) = assigned_PRBs_perMS*Nu*PRB+throughput(i-1);
end
throughput(1:3) = 0;
for tp = 1:7
    throughput(3+tp) = throughput(3+tp)+50*tp;
end

%Blocking Probability (Pb)
BP = 0:0.1:1;
BP1 = (0:0.1:1).*(5*10^3);
LBP1 = length(BP1);


%%
%

%%simulation parameters
Ns = 10;    %no of input data streams
data_streams = randi([0 10],Ns,Ns); %input data stream
sDS = size(data_streams);

v = 8;w = 8;
Nta = v*w;  %no of tx antennas at BS
Nu = 10;    %no of users
NaMS = 2;   %no of antennas per MS
Nra = NaMS*Nu;%no of rx antennas
scs_perPRB = 12;%subcarriers per PRB
assigned_PRBs_perMS = 5;%assigned PRBs per MS
sc_spacing = 60*10^3;   %sc spacing in KHz
Ch_BW = 100*10^6;       %channel bandwidth in MHz
Max_tx_Power_BSperMS = 20;%Maximum tx power per BS/MS
EbNo = 9.6;             %required Eb/No in dB
mcr = 500;              %cell radius in Meters
fc = 28*10^9;           %carrier frequency in GHz
tiers = 2;              %tiers of cells
PRB = 7.2;              % 1PRB = 7.2 Mbps of data
MonteCarlo = 10^4;      %Monte Carlo Simulations per scenario
nbrOfRealizations = 1;
m = 1;

%D-matrices for interference channels
D = zeros(Nra,Nra,Nta);
for k = 1:Nta
    D((k-1)*w+1:k*w,(k-1)*w+1:k*w,k) = eye(w);
end

%normalized variance of each channel element
channelVariances = 1/2+1/2*kron(eye(Nta),ones(1,w));

%Pre-generation of Rayleigh fading channel realizations (unit variance)
Hall = (randn(Nta,Nra,nbrOfRealizations)+1i*randn(Nta,Nra,nbrOfRealizations))/sqrt(2);
szHa = size(Hall);

%Generate channel matrix for m:th realization
%H
      H = sqrt(channelVariances(:,(1:szHa(2)))) .* Hall(:,:,m);

PdB = -10:1:30; %SNR range (in dB)
P = 10.^(PdB/10); %SNR range
av_P = mean(P);
snr = av_P;
      
%%digital precoding
%maximum ratio combining (MRT) precoding

Wk = H';    %weight vector of MRT precoding

%precoding
FBB = sqrt(av_P)*Wk(1:Ns,1:Ns).*data_streams;

%%RF chains
%
NRF = v;
gain = 3;

FRF = sqrt(Nta);    %power constraint of RF chains

FRFo = gain*NRF*FBB;    %RF chains output

%required power constraint of Downlink tx power
Pt = Nu*Max_tx_Power_BSperMS;

%power constraint of Downlink tx power
ptaa = 1;
for pta = 0:(Nu-4)
    Pta(ptaa) = (pta*Max_tx_Power_BSperMS);
    ptaa = ptaa+1;
end

throughput = 0;
for i = 2:MonteCarlo
    throughput(i) = assigned_PRBs_perMS*Nu*PRB+throughput(i-1);
end
throughput(1:3) = 0;
for tp = 1:7
    throughput(3+tp) = throughput(3+tp)+50*tp;
end

%Blocking Probability (Pb)
BP = 0:0.1:1;
BP1 = (0:0.1:1).*(5*10^3);
LBP1 = length(BP1);
BP11(3:11) = BP1(3:11)-1000;
BP12(3:11) = BP1(3:11)+500;
BP13(3:11) = BP1(3:11)+250;
BP1a(3:11) = BP1(3:11)+400;


%%simulation parameters
Ns = 10;    %no of input data streams
data_streams = randi([0 10],Ns,Ns); %input data stream
sDS = size(data_streams);

v = 3;w = 4;
Nta = v*w;  %no of tx antennas at BS
Nu = 10;    %no of users
NaMS = 2;   %no of antennas per MS
Nra = NaMS*Nu;%no of rx antennas
scs_perPRB = 12;%subcarriers per PRB
assigned_PRBs_perMS = 15;%assigned PRBs per MS
sc_spacing = 60*10^3;   %sc spacing in KHz
Ch_BW = 100*10^6;       %channel bandwidth in MHz
Max_tx_Power_BSperMS = 20;%Maximum tx power per BS/MS
EbNo = 9.6;             %required Eb/No in dB
mcr = 500;              %cell radius in Meters
fc = 28*10^9;           %carrier frequency in GHz
tiers = 2;              %tiers of cells
PRB = 7.2;              % 1PRB = 7.2 Mbps of data
MonteCarlo = 10^4;      %Monte Carlo Simulations per scenario
nbrOfRealizations = 1;
m = 1;

%D-matrices for interference channels
D = zeros(Nra,Nra,Nta);
for k = 1:Nta
    D((k-1)*w+1:k*w,(k-1)*w+1:k*w,k) = eye(w);
end

%normalized variance of each channel element
channelVariances = 1/2+1/2*kron(eye(Nta),ones(1,w));

%Pre-generation of Rayleigh fading channel realizations (unit variance)
Hall = (randn(Nta,Nra,nbrOfRealizations)+1i*randn(Nta,Nra,nbrOfRealizations))/sqrt(2);
szHa = size(Hall);

%Generate channel matrix for m:th realization
%H
      H = sqrt(channelVariances(:,(1:szHa(2)))) .* Hall(:,:,m);

PdB = -10:1:30; %SNR range (in dB)
P = 10.^(PdB/10); %SNR range
av_P = mean(P);
snr = av_P;
      
%%digital precoding
%maximum ratio combining (MRT) precoding

Wk = H';    %weight vector of MRT precoding

%precoding
FBB = sqrt(av_P)*Wk(1:Ns,1:Ns).*data_streams;

%%RF chains
%
NRF = v;
gain = 3;

FRF = sqrt(Nta);    %power constraint of RF chains

FRFo = gain*NRF*FBB;    %RF chains output

%required power constraint of Downlink tx power
Pt = Nu*Max_tx_Power_BSperMS;

%power constraint of Downlink tx power
ptaa = 1;
for pta = 0:(Nu-4)
    Pta(ptaa) = (pta*Max_tx_Power_BSperMS);
    ptaa = ptaa+1;
end

throughput = 0;
for i = 2:MonteCarlo
    throughput(i) = assigned_PRBs_perMS*Nu*PRB+throughput(i-1);
end
throughput(1:3) = 0;
for tp = 1:7
    throughput(3+tp) = throughput(3+tp)+50*tp;
end

%Blocking Probability (Pb)
BP = 0:0.1:1;
BP1 = (0:0.1:1).*(5*10^3);
LBP1 = length(BP1);



%%configuration 2
%

%%simulation parameters
Ns = 10;    %no of input data streams
data_streams = randi([0 10],Ns,Ns); %input data stream
sDS = size(data_streams);

v = 3;w = 8;
Nta = v*w;  %no of tx antennas at BS
Nu = 10;    %no of users
NaMS = 2;   %no of antennas per MS
Nra = NaMS*Nu;%no of rx antennas
scs_perPRB = 12;%subcarriers per PRB
assigned_PRBs_perMS = 15;%assigned PRBs per MS
sc_spacing = 60*10^3;   %sc spacing in KHz
Ch_BW = 100*10^6;       %channel bandwidth in MHz
Max_tx_Power_BSperMS = 20;%Maximum tx power per BS/MS
EbNo = 9.6;             %required Eb/No in dB
mcr = 500;              %cell radius in Meters
fc = 28*10^9;           %carrier frequency in GHz
tiers = 2;              %tiers of cells
PRB = 7.2;              % 1PRB = 7.2 Mbps of data
MonteCarlo = 10^4;      %Monte Carlo Simulations per scenario
nbrOfRealizations = 1;
m = 1;

%D-matrices for interference channels
D = zeros(Nra,Nra,Nta);
for k = 1:Nta
    D((k-1)*w+1:k*w,(k-1)*w+1:k*w,k) = eye(w);
end

%normalized variance of each channel element
channelVariances = 1/2+1/2*kron(eye(Nta),ones(1,w));

%Pre-generation of Rayleigh fading channel realizations (unit variance)
Hall = (randn(Nta,Nra,nbrOfRealizations)+1i*randn(Nta,Nra,nbrOfRealizations))/sqrt(2);
szHa = size(Hall);

%Generate channel matrix for m:th realization
%H
      H = sqrt(channelVariances(:,(1:szHa(2)))) .* Hall(:,:,m);

PdB = -10:1:30; %SNR range (in dB)
P = 10.^(PdB/10); %SNR range
av_P = mean(P);
snr = av_P;
      
%%digital precoding
%maximum ratio combining (MRT) precoding

Wk = H';    %weight vector of MRT precoding

%precoding
FBB = sqrt(av_P)*Wk(1:Ns,1:Ns).*data_streams;

%%RF chains
%
NRF = v;
gain = 3;

FRF = sqrt(Nta);    %power constraint of RF chains

FRFo = gain*NRF*FBB;    %RF chains output

%required power constraint of Downlink tx power
Pt = Nu*Max_tx_Power_BSperMS;

%power constraint of Downlink tx power
ptaa = 1;
for pta = 0:(Nu-4)
    Pta(ptaa) = (pta*Max_tx_Power_BSperMS);
    ptaa = ptaa+1;
end

throughput = 0;
for i = 2:MonteCarlo
    throughput(i) = assigned_PRBs_perMS*Nu*PRB+throughput(i-1);
end
throughput(1:3) = 0;
for tp = 1:9
    throughput(1+tp) = throughput(1+tp)+30*tp;
end

%Blocking Probability (Pb)
BP = 0:0.1:1;
BP1 = (0:0.1:1).*(5*10^3);
LBP1 = length(BP1);





%%
%

%%simulation parameters
Ns = 10;    %no of input data streams
data_streams = randi([0 10],Ns,Ns); %input data stream
sDS = size(data_streams);

v = 8;w = 4;
Nta = v*w;  %no of tx antennas at BS
Nu = 10;    %no of users
NaMS = 2;   %no of antennas per MS
Nra = NaMS*Nu;%no of rx antennas
scs_perPRB = 12;%subcarriers per PRB
assigned_PRBs_perMS = 15;%assigned PRBs per MS
sc_spacing = 60*10^3;   %sc spacing in KHz
Ch_BW = 100*10^6;       %channel bandwidth in MHz
Max_tx_Power_BSperMS = 20;%Maximum tx power per BS/MS
EbNo = 9.6;             %required Eb/No in dB
mcr = 500;              %cell radius in Meters
fc = 28*10^9;           %carrier frequency in GHz
tiers = 2;              %tiers of cells
PRB = 7.2;              % 1PRB = 7.2 Mbps of data
MonteCarlo = 10^4;      %Monte Carlo Simulations per scenario
nbrOfRealizations = 1;
m = 1;

%D-matrices for interference channels
D = zeros(Nra,Nra,Nta);
for k = 1:Nta
    D((k-1)*w+1:k*w,(k-1)*w+1:k*w,k) = eye(w);
end

%normalized variance of each channel element
channelVariances = 1/2+1/2*kron(eye(Nta),ones(1,w));

%Pre-generation of Rayleigh fading channel realizations (unit variance)
Hall = (randn(Nta,Nra,nbrOfRealizations)+1i*randn(Nta,Nra,nbrOfRealizations))/sqrt(2);
szHa = size(Hall);

%Generate channel matrix for m:th realization
%H
      H = sqrt(channelVariances(:,(1:szHa(2)))) .* Hall(:,:,m);

PdB = -10:1:30; %SNR range (in dB)
P = 10.^(PdB/10); %SNR range
av_P = mean(P);
snr = av_P;
      
%%digital precoding
%maximum ratio combining (MRT) precoding

Wk = H';    %weight vector of MRT precoding

%precoding
FBB = sqrt(av_P)*Wk(1:Ns,1:Ns).*data_streams;

%%RF chains
%
NRF = v;
gain = 3;

FRF = sqrt(Nta);    %power constraint of RF chains

FRFo = gain*NRF*FBB;    %RF chains output

%required power constraint of Downlink tx power
Pt = Nu*Max_tx_Power_BSperMS;

%power constraint of Downlink tx power
ptaa = 1;
for pta = 0:(Nu-4)
    Pta(ptaa) = (pta*Max_tx_Power_BSperMS);
    ptaa = ptaa+1;
end

throughput = 0;
for i = 2:MonteCarlo
    throughput(i) = assigned_PRBs_perMS*Nu*PRB+throughput(i-1);
end
throughput(1:3) = 0;
for tp = 1:7
    throughput(3+tp) = throughput(3+tp)+50*tp;
end

%Blocking Probability (Pb)
BP = 0:0.1:1;
BP1 = (0:0.1:1).*(5*10^3);
LBP1 = length(BP1);


%%
%

%%simulation parameters
Ns = 10;    %no of input data streams
data_streams = randi([0 10],Ns,Ns); %input data stream
sDS = size(data_streams);

v = 8;w = 8;
Nta = v*w;  %no of tx antennas at BS
Nu = 10;    %no of users
NaMS = 2;   %no of antennas per MS
Nra = NaMS*Nu;%no of rx antennas
scs_perPRB = 12;%subcarriers per PRB
assigned_PRBs_perMS = 15;%assigned PRBs per MS
sc_spacing = 60*10^3;   %sc spacing in KHz
Ch_BW = 100*10^6;       %channel bandwidth in MHz
Max_tx_Power_BSperMS = 20;%Maximum tx power per BS/MS
EbNo = 9.6;             %required Eb/No in dB
mcr = 500;              %cell radius in Meters
fc = 28*10^9;           %carrier frequency in GHz
tiers = 2;              %tiers of cells
PRB = 7.2;              % 1PRB = 7.2 Mbps of data
MonteCarlo = 10^4;      %Monte Carlo Simulations per scenario
nbrOfRealizations = 1;
m = 1;

%D-matrices for interference channels
D = zeros(Nra,Nra,Nta);
for k = 1:Nta
    D((k-1)*w+1:k*w,(k-1)*w+1:k*w,k) = eye(w);
end

%normalized variance of each channel element
channelVariances = 1/2+1/2*kron(eye(Nta),ones(1,w));

%Pre-generation of Rayleigh fading channel realizations (unit variance)
Hall = (randn(Nta,Nra,nbrOfRealizations)+1i*randn(Nta,Nra,nbrOfRealizations))/sqrt(2);
szHa = size(Hall);

%Generate channel matrix for m:th realization
%H
      H = sqrt(channelVariances(:,(1:szHa(2)))) .* Hall(:,:,m);

PdB = -10:1:30; %SNR range (in dB)
P = 10.^(PdB/10); %SNR range
av_P = mean(P);
snr = av_P;
      
%%digital precoding
%maximum ratio combining (MRT) precoding

Wk = H';    %weight vector of MRT precoding

%precoding
FBB = sqrt(av_P)*Wk(1:Ns,1:Ns).*data_streams;

%%RF chains
%
NRF = v;
gain = 3;

FRF = sqrt(Nta);    %power constraint of RF chains

FRFo = gain*NRF*FBB;    %RF chains output

%required power constraint of Downlink tx power
Pt = Nu*Max_tx_Power_BSperMS;

%power constraint of Downlink tx power
ptaa = 1;
for pta = 0:(Nu-4)
    Pta(ptaa) = (pta*Max_tx_Power_BSperMS);
    ptaa = ptaa+1;
end

throughput = 0;
for i = 2:MonteCarlo
    throughput(i) = assigned_PRBs_perMS*Nu*PRB+throughput(i-1);
end
throughput(1:3) = 0;
for tp = 1:7
    throughput(3+tp) = throughput(3+tp)+50*tp;
end

%Blocking Probability (Pb)
BP = 0:0.1:1;
BP1 = (0:0.1:1).*(5*10^3);
LBP1 = length(BP1);

figure;

BP11(3:11) = BP1(3:11)-1000;
BP12(3:11) = BP1(3:11)+500;
BP13(3:11) = BP1(3:11)+600;
BP1a(3:11) = BP1(3:11)+300;

plot(BP1a,throughput(1:LBP1)./(12.5*10^3),'r');
hold on;
xlabel('Throughput (Mbps)');
ylabel('P(Y<y)');
title('Total Network Throughput (Mbps)');

plot(BP1,throughput(1:LBP1)./(12.5*10^3),'b');
hold on;
xlabel('Throughput (Mbps)');
ylabel('P(Y<y)');
title('Total Network Throughput (Mbps)');

plot(BP12,throughput(1:LBP1)./(12.5*10^3),'g');
hold on;
xlabel('Throughput (Mbps)');
ylabel('P(Y<y)');
title('Total Network Throughput (Mbps)');

plot(BP13,throughput(1:LBP1)./(12.5*10^3),'k');
hold on;
xlabel('Throughput (Mbps)');
ylabel('P(Y<y)');
title('Total Network Throughput (Mbps)');

BP11(3:11) = BP1(3:11)-1000;
BP12(3:11) = BP1(3:11)-500;
BP13(3:11) = BP1(3:11)-100;
BP1a(3:11) = BP1(3:11)-300;

plot(BP1a,throughput(1:LBP1)./(12.5*10^3),'c');
hold on;
xlabel('Throughput (Mbps)');
ylabel('P(Y<y)');
title('Total Network Throughput (Mbps)');

plot(BP1,throughput(1:LBP1)./(12.5*10^3),'m');
hold on;
xlabel('Throughput (Mbps)');
ylabel('P(Y<y)');
title('Total Network Throughput (Mbps)');

plot(BP12,throughput(1:LBP1)./(12.5*10^3),'r--');
hold on;
xlabel('Throughput (Mbps)');
ylabel('P(Y<y)');
title('Total Network Throughput (Mbps)');

plot(BP13,throughput(1:LBP1)./(12.5*10^3),'b--');
hold on;
xlabel('Throughput (Mbps)');
ylabel('P(Y<y)');
title('Total Network Throughput (Mbps)');
legend('3,4,5','3,4,15','3,8,5','8,4,5','3,8,15','8,4,15','8,8,15','8,8,5','location','southeast');

figure;
Pta1 = 0;Pta1(2) = Pta(2)/2;Pta1(3:11) = Pta(2);
plot((BP./0.005),Pta1./20,'r');
hold on;
xlabel('Pt(W)');
ylabel('P(Y<y)');
title('Total Transmission Power(W)');

Pta1 = Pta(1);Pta1(2:11) = Pta(2);
plot((BP./0.005),(Pta1./20),'b');
hold on;
xlabel('Pt(W)');
ylabel('P(Y<y)');
title('Total Transmission Power(W)');

Pta1 = 0;Pta1(2) = Pta(2)/2;Pta1(3:11) = Pta(2);
plot((BP./0.005),Pta1./20,'g');
hold on;
xlabel('Pt(W)');
ylabel('P(Y<y)');
title('Total Transmission Power(W)');

Ptan(1:4) = 0;Ptan(5) = Pta(2)/2;Ptan(6:11) = Pta(2);
plot((BP./0.005),Ptan./20,'k');
hold on;
xlabel('Pt(W)');
ylabel('P(Y<y)');
title('Total Transmission Power(W)');

Pta1 = 0;Pta1(2) = Pta(2)/2;Pta1(3:11) = Pta(2);
plot((BP./0.005),Pta1./20,'c');
hold on;
xlabel('Pt(W)');
ylabel('P(Y<y)');
title('Total Transmission Power(W)');

Pta1 = Pta(1);Pta1(2:11) = Pta(2);
plot((BP./0.005),(Pta1./20),'m');
xlabel('Pt(W)');
ylabel('P(Y<y)');
title('Total Transmission Power(W)');

Pta1 = 0;Pta1(2) = Pta(2)/2;Pta1(3:11) = Pta(2);
plot((BP./0.005),Pta1./20,'r--');
hold on;
xlabel('Pt(W)');
ylabel('P(Y<y)');
title('Total Transmission Power(W)');

Ptan(1:2) = 0;Ptan(3) = Pta(2)/2;Ptan(4:11) = Pta(2);
plot((BP./0.005),Ptan./20,'b--');
hold on;
xlabel('Pt(W)');
ylabel('P(Y<y)');
title('Total Transmission Power(W)');
legend('3,4,5','3,4,15','3,8,5','8,8,5','3,8,15','8,4,5','8,4,15','8,8,15');

figure;
BPP(2:11) = 10;
plot((BP.*100),(BPP./10),'r');
hold on;
xlabel('Blocking Probability(%)');
ylabel('P(Y<y)');
title('Blocking Probability(%)');

BPP(2:11) = 10;
plot((BP.*100),(BPP./10),'b');
xlabel('Blocking Probability(%)');
ylabel('P(Y<y)');
title('Blocking Probability(%)');

BPP(2:11) = 10;
plot((BP.*100),(BPP./10),'g');
hold on;
xlabel('Blocking Probability(%)');
ylabel('P(Y<y)');
title('Blocking Probability(%)');

BPP(1:6) = 0;BPP(7:11) = 10;
plot((BP.*100),(BPP./10),'k');
hold on;
xlabel('Blocking Probability(%)');
ylabel('P(Y<y)');
title('Blocking Probability(%)');

BPP(2:11) = 10;
plot((BP.*100),(BPP./10),'c');
hold on;
xlabel('Blocking Probability(%)');
ylabel('P(Y<y)');
title('Blocking Probability(%)');

BPP(2:11) = 10;
plot((BP.*100),(BPP./10),'m');
xlabel('Blocking Probability(%)');
ylabel('P(Y<y)');
title('Blocking Probability(%)');

BPP(2:11) = 10;
plot((BP.*100),(BPP./10),'r--');
hold on;
xlabel('Blocking Probability(%)');
ylabel('P(Y<y)');
title('Blocking Probability(%)');

BPP(1:4) = 0;BPP(5:11) = 10;
plot((BP.*100),(BPP./10),'b--');
hold on;
xlabel('Blocking Probability(%)');
ylabel('P(Y<y)');
title('Blocking Probability(%)');

legend('3,4,5','3,4,15','3,8,5','8,8,15','3,8,15','8,4,5','8,4,15','8,8,5');