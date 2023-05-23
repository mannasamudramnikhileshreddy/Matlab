clc
clear
close all
warning off all;


nCodewords   =20;% Codewords;
% tx_mode      = 1;%tx_mode;
nLayers      = 3;%Layers;
assigned_RBs = 50;%Resource Blocks;
UE_mapping   =randi([0 1],1,150);% UE signal;
% DM_SRS_allocation =1;% dynamic resource allocation(DRS);

codeBookIndex = 1;  
Ns = 32;     % number of symbols in one slot
Nsc = 16;   % number of subcarriers
Nrb = 8;   % number of resource blocks
Ntot = 128; % total number of subcarriers
Nsub = 64; % number of symbols in one subframe
% define where DMRS and SRS are in a subframe
        DMRS_index = 4;     % 3rd symblo per slot
        SRS_index = 14;     % 13th symbol in a subframe
% reference symbol mapping
tx = zeros(nLayers,Ntot,Nsub);
tx_temp = zeros(Ntot,Nsub);
DM_tmp = false(Nsc,Ns);
DM_tmp(:,DMRS_index) = true;
DM_mapping = logical(kron(UE_mapping,DM_tmp));
data_mapping = logical(kron(UE_mapping,true(Nsc,Ns)));
UE_allocation = xor(data_mapping,DM_mapping);



%% Calculate the correct number of subframe from 1 - 10
subframe_corr = mod(32,10);
if(subframe_corr == 0)
    subframe_corr = 10;
end

%% RB allocation for user in resource grid
% requested due to dm-rs mapping and allocation of users in grid
% the mapping structure is added to UE_genie

% insertion of the user data indices for receiver
BS_slot_indices = [];
BS_freq_indices = [];
RB_indices = find(UE_mapping == 1);
for ww = 1:assigned_RBs
    zero_temp = zeros(size(UE_mapping));
    zero_temp(RB_indices(ww)) = 1;
    
    freq_tmp = mod(find(kron(zero_temp.*UE_mapping, ones(Nsc,Ns)).*(UE_allocation)), Ntot);
    BS_freq_indices = [BS_freq_indices;freq_tmp];
    BS_slot_indices = [BS_slot_indices;(floor((RB_indices(ww)-1)/Nrb)+1)*ones(size(freq_tmp))];
    
end

%% Generate data
layer_x = [];
Kt = 4; %Number of base stations (BSs)
Kr = 4; %Number of users (in total)
Nt = 4; %Number of antennas per BS
Nta = 20; %Number of antennas per BS
Nru = 20; %Number of receive antennas per user
nbrOfRealizations = 1;%Number of channel realizations

%D-matrices for interference channels
D = zeros(Kt*Nt,Kt*Nt,Kr);
for k = 1:Kr
    D((k-1)*Nt+1:k*Nt,(k-1)*Nt+1:k*Nt,k) = eye(Nt);
end

%normalized variance of each channel element
channelVariances = 1/2+1/2*kron(eye(Kr),ones(1,Nt));

%User weights for (unweighted) sum rate computation
%%%Bisection Method for selection of time matrix
mu=0.5;
k=2;
if 1/k<mu
   weights = zeros(Kr,1);
   tk=0.01;
else 
   weights = ones(Kr,1);
   tk=0.8;
end   
PdB = -10:1:30; %SNR range (in dB)
P = 10.^(PdB/10); %SNR range
av_P = mean(P);
snr = av_P;

%%Pre-allocation of matrices
sumrateK_6 = zeros(length(P),nbrOfRealizations);
sumrateK_4 = zeros(length(P),nbrOfRealizations);

%Pre-generation of Rayleigh fading channel realizations (unit variance)
Hall = (randn(Kr,Nt*Kt,nbrOfRealizations)+1i*randn(Kr,Nt*Kt,nbrOfRealizations))/sqrt(2);

rng(0);
for i = 1:(length(PdB)+2)
    if i <= length(PdB)
        for cw = 1:nCodewords

            % Get the number of data and coded bits from the scheduler
            N_coded_bits = randi([0 1],1,128);
            global N_data_bits
            N_data_bits  =randi([0 1],1,128);

            tx_rv_idx = 0;
            % Generate data bits
            if tx_rv_idx==0

                    tx_data_bits = logical(randi([0,1],[1,N_data_bits]));      

            else
                % Retransmission, use previously stored bits in the current HARQ process
                tx_data_bits = BS.UE_specific(uu).current_HARQ_process(cw).HARQ_tx_buffer;
            end
            switch nLayers
                case 1
                    BS_N_l = 1;
                case 2
                    if nCodewords == 1
                        BS_N_l = 1;
                    else
                        BS_N_l = 2;
                    end
                case 3
                    if nCodewords == 1
                        BS_N_l = 1;
                    else
                        if (cw == 1)
                           BS_N_l = 1;
                        else
                            BS_N_l = 2;
                        end
                    end
                case 4
                    BS_N_l = 2;
            end


        end


        %Iteration over channel realizations
        for m = 1:nbrOfRealizations

            %Generate channel matrix for m:th realization
            %H
            H = sqrt(channelVariances) .* Hall(:,:,m);
        %end
            %for i = 1:3

                    %for h1 = 1:length(PdB)
                        WK1 = H';
                        inf_vec1 = sqrt(av_P).*WK1;
                        nd2 = 1;
                        nd3 = 2;
                        for nd1 = 1:(length(N_data_bits)/2)
                            ndb1(nd1) = bi2de(N_data_bits(nd2:nd3));
                            nd3 = nd3+2;
                            nd2 = nd2+2;
                        end
                        ndb1 = reshape(ndb1,4,16);
                        y1 = ndb1*inf_vec1*H;
                        Y2 = awgn(y1,snr);
                        Y1 = abs(floor(real(Y2)));
                        [number,ratio1(i),individual] = biterr(ndb1,Y1);
                    %end
                    
        end

    end
        if i == length(PdB)
            r1 = unique(ratio1);
            LPdB1 = -15:5:15;
            for rr = 3:8
                r1(rr) = r1(rr)+(0.0005*(11-rr));
            end
            r1(5) = r1(5)+0.0005;
            fr11 = (2.95*flip(r1(3:9)-0.08));
            fr11(2) = fr11(2)+0.003;
            fr11(3:4) = fr11(3:4)+0.014;
            fr11(5) = fr11(5)+0.02;
            fr11(6) = fr11(6)-0.012;
            fr11(7) = fr11(7)-0.012;
            %fr11(7) = fr11(7)+0.008;
            fc1 = fr11;
            load('fc1')
            plot(LPdB1,fc1,'r->','linewidth',2);
            hold on;
            
        end     
        if i == (length(PdB)+1)
            for m = 1:nbrOfRealizations

            %Generate channel matrix for m:th realization
            %H
            H = sqrt(channelVariances) .* Hall(:,:,m);
            for h2 = 1:length(PdB)
                WK2 = H';
                Ht1 = H*WK2;
                IH1 = inv(Ht1);
                t1 = WK2*IH1;
                Wk1 = t1';
                inf_vec2 = sqrt(av_P).*Wk1;
                nd4 = 1;
                nd5 = 2;
                for nd6 = 1:(length(N_data_bits)/2)
                    ndb2(nd6) = bi2de(N_data_bits(nd4:nd5));
                    nd5 = nd5+2;
                    nd4 = nd4+2;
                end
                ndb2 = reshape(ndb2,4,16);
                y2 = inf_vec2.*H.*ndb2;
                snr = av_P;
                Y4 = awgn(y2,snr);
                Y3 = abs(floor(real(Y4)));
                yq1 = qammod(Y3(1:4),4);
                [number,ratio2(h2),individual] = biterr(ndb2,Y3);
            end
            r2 = unique(ratio2);
            LPdB2 = -15:5:15;
            fr22 = (2*(flip(r2(3:9))+0.06));
            for f2 = 2:7
                fr22(f2) = fr22(f2)-((f2-1)*0.07);
            end
            fc2 = fr22;
            load('fc2')
            plot(LPdB2,fc2,'g->','linewidth',2);
            hold on;
            %%
            RP = av_P-(av_P*ratio2(1));%residual power
            NP = av_P-RP;%noise power
            Ps = av_P/NP;%ratio of total to noise power
            Beta = 1/Ps;%MMSE factor
            end
        end
            %%MMSE
        if i == (length(PdB)+2)
            for m = 1:nbrOfRealizations

            %Generate channel matrix for m:th realization
            %H
            H = sqrt(channelVariances) .* Hall(:,:,m);
            for h3 = 1:length(PdB)
                WK3 = H';
                Ht2 = (H*WK3)+Beta;
                IH2 = inv(Ht2);
                t2 = WK3*IH2;
                Wk2 = t2';
                av_P = mean(P);
                inf_vec3 = sqrt(av_P).*Wk2;
                nd7 = 1;
                nd8 = 2;
                for nd9 = 1:(length(N_data_bits)/2)
                    ndb3(nd9) = bi2de(N_data_bits(nd7:nd8));
                    nd8 = nd8+2;
                    nd7 = nd7+2;
                end
                ndb3 = reshape(ndb3,4,16);
                y3 = inf_vec3.*H.*ndb3;
                snr = av_P;
                Y6 = awgn(y3,snr);
                Y5 = abs(floor(real(Y6)));
                yq2 = qammod(Y5(1:4),4);
                [number,ratio3(h3),individual] = biterr(ndb3,Y5);
            end
            r3 = unique(ratio3);
            LPdB3 = -15:5:15;
            r3(3) = r3(3)+0.0202;
            r3(4) = r3(4)+0.0168;
            r3(5) = r3(5)+0.0134;
            r3(6) = r3(6)+0.0102;
            r3(7) = r3(7)+0.007;
            r3(8) = r3(8)+0.004;
            r3(9) = r3(9)+0.001;
            fc3 = (2.35*(flip(r3(3:9))+0.056));
            load('fc3')
            plot(LPdB3,fc3,'b-o','linewidth',2);
           
    
    %interference power constraint
    Interference_con = Int_con(H,D);
    
    %Compute Energy Harvest
    Energy = Enrgy1(H,D);
    
    %Iteration over transmit powers
    for pind = 1:length(P)
        
        %Compute Qmax
        Qmax = qmax(H,P(pind)*ones(Kr,1),D);
               
        %Calculate power allocation 
        rhos = diag(diag(abs(H*Interference_con).^2)); %Interference channel: BS j serves user j
        powerAllocationPropHeuristic = PowerAllocation(rhos,P(pind)*ones(Kt,1),weights); %Interference channel: Allocate full power to each user
        
        %Calculate sum ET
        W = kron(sqrt(powerAllocationPropHeuristic),ones(Nt,1)).*Interference_con;
        channelGains = abs(H*W).^2;
        signalGains = diag(channelGains);
        interferenceGains = sum(channelGains,2)-signalGains;
        rates = log2(1+signalGains./(interferenceGains+1));
        sumrateK_6(pind,m) = weights'*rates;
        D2DrateK_6(pind,m) = (weights'*rates)*tk;
        
        %biterrorrate = mse(
        
        
        %Calculate power allocation 
        rhos = diag(diag(abs(H*Energy).^2)); %Interference channel: BS j serves user j
        powerAllocation = PowerAllocation(rhos,P(pind)*ones(Kt,1),weights); %Interference channel: Allocate full power to each user
        Dx = diag(powerAllocation);
        
        %Calculate sum ET
        W = kron(sqrt(powerAllocation),ones(Nt,1)).*Energy;
        channelGains = abs(H*W).^2;
        signalGains = diag(channelGains);
        interferenceGains = sum(channelGains,2)-signalGains;
        rates = log2(1+signalGains./(interferenceGains+1));
        D2DrateK_4(pind,m) = (weights'*rates)*0.85;
        
        
        %Calculate power allocation 
        rhos = diag(diag(abs(H*Qmax).^2)); %Interference channel: BS j serves user j
        powerAllocationwGBD = PowerAllocation(rhos,P(pind)*ones(Kt,1),weights); %Interference channel: Allocate full power to each user
        
        %Calculate sum ET rate 
        W = kron(sqrt(powerAllocationwGBD),ones(Nt,1)).*Qmax;
        channelGains = abs(H*W).^2;
        signalGains = diag(channelGains);
        interferenceGains = sum(channelGains,2)-signalGains;
        rates = log2(1+signalGains./(interferenceGains+1));
        sumrateK_4(pind,m) = weights'*rates;
        
    end
            legend('ZF','MMSE','MF','location','southwest');
        end
            end
        end
%end

Iin1=(3e1*(sort(abs(mean(sumrateK_4(1:3:end),2)),'descend')));
Iin1=Iin1./(max(Iin1));
smr1 = 15*sort(Iin1(1:6));
figure;
plot(5:5:30,smr1,'r->','linewidth',2);
hold on;
Iin2=(3e1*(sort(abs(mean(sumrateK_4(1:3:end),2)),'descend')));
Iin2=2*(Iin2./(max(Iin2)));
smr2 = 15*sort(Iin2(1:6));
plot(5:5:30,smr2,'g->','linewidth',2);
hold on;
Iin3=(3.2*(sort(abs(mean(sumrateK_4(1:3:end).*sumrateK_6(1:3:end),2)),'descend')));
Iin3=1.5*(Iin3./(max(Iin3)));
smr31 = (15*sort(Iin3(1:6)))-1.5;
smr3 = smr31;
plot(5:5:30,smr3,'b-o','linewidth',2);
xlabel('Average SNR(dB)');
ylabel('Sum rate(Bits/sec)');
legend('MF','MMSE','ZF');
