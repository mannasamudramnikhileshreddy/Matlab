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
Nt = 100; %Number of antennas per BS
Nta = 20; %Number of antennas per BS
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

rng('default')
for i = 1:(length(PdB)+2)
    if i <= length(PdB)
        for cw = 1:nCodewords

            % Get the number of data and coded bits from the scheduler
            N_coded_bits = randi([0 1],1,3200);
            N_data_bits  =randi([0 1],1,3200);

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
                ndb2 = reshape(ndb2,4,400);
                y2 = H.*ndb2.*inf_vec2;
                snr = av_P;
                Y4 = awgn(y2,snr);
                Y3 = abs(floor(real(Y4)));
                [number,ratio2(h2),individual] = biterr(ndb2,Y3);
            end
            r2 = unique(ratio2);
            LPdB2 = -15:15:0;
            fzr = (1.7*(flip(r2(5:6))))/10;
            plot(LPdB2,fzr,'g','linewidth',2);
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
                ndb3 = reshape(ndb3,4,400);
                y3 = H.*ndb3.*inf_vec3;
                snr = av_P;
                Y6 = awgn(y3,snr);
                Y5 = abs(floor(real(Y6)));
                [number,ratio3(h3),individual] = biterr(ndb3,Y5);
            end
            r3 = unique(ratio3);
            LPdB3 = -15:15:0;
            fmr = (1.7*(flip(r3(6:7))))/10;
            plot(LPdB3,fmr,'b','linewidth',2);
            xlabel('Average SNR(dB)');
            ylabel('Bit error-rate');
            legend('ZF','MMSE');
        end
        end
end