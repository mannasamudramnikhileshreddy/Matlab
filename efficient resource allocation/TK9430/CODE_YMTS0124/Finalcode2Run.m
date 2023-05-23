clc
clear all
close all
s = rng(21); 
Noise_dBm =	-114; % AWG Noise for eNB, CT & DT in dBm/Hz
Noise_Figure_dB = 5; % Amplifier noise addition to thermal noise
CT_SNR_Target_dB = 12; %dB
DT_SNR_Target_dB = 25; %dB
CT_SINR_min_dB	= 6; %dB
DT_SINR_min_dB	= 15; %dB
Pd = 21;
CT_P_Tx_Min_dBm = -40;
DT_P_Tx_Max_dBm = 10; 
CT_P_Rx_min_dBm = -98;
DT_P_Rx_min_dBm = -78; 
Carrier_Frequency =	3.6; % GHz
Bandwidth_kHz = 180; 
CT_Tx_Power_UIP_dBm = 50; 
              
%%% Conversion of above parameters to ratio or Watts values
Thermal_Noise_Total_Watts =	Bandwidth_kHz*10^3*(db2pow(Noise_dBm-30)); % Total AWG Noise for eNB, CUE & DUE in Watts
Thermal_Noise_Total_dBm = pow2db(Thermal_Noise_Total_Watts)+30; % Total AWG Noise for eNB, CUE & DUE in dBm
Noise_Total_dBm = Thermal_Noise_Total_dBm + Noise_Figure_dB; % Total added noise including amplifier noise
Noise_Total_Watts =	db2pow(Noise_Total_dBm); % Total AWG Noise for eNB, CUE & DUE in Watts including noise factor
CT_SNR_Target = db2pow(CT_SNR_Target_dB); %Real ratio value
DT_SNR_Target = db2pow(DT_SNR_Target_dB); %Real ratio value
CT_SINR_min = db2pow(CT_SINR_min_dB); % Real ratio value
DT_SINR_min = db2pow(DT_SINR_min_dB); % Real ratio value
CT_P_Tx_Max = db2pow(Pd-30); % Maximum CT transmit power in Watts
CT_P_Tx_Min = db2pow(CT_P_Tx_Min_dBm-30); % Minimum CT transmit power in Watts
DT_P_Tx_Max = db2pow(DT_P_Tx_Max_dBm-30); % Maximum CT transmit power in Watts
CT_P_Rx_min = db2pow(CT_P_Rx_min_dBm-30); % Minimum CT receiver power in Watts
DT_P_Rx_min = db2pow(DT_P_Rx_min_dBm-30); % Minimum D2D ProSe receiver power in Watts

CT_Exp = 0.01; % Pathloss Exponent between devices and base station (eNB) - Alpha_o
DT_Exp=input('Enter alpha either 2.5 or 4'); % Pathloss Exponent between communication among devices - Alpha_d
% DT_Exp = 4;
CT_Prop_Const = (10^-2.27)*(Carrier_Frequency^-2.6); % Cellular Propagation Constant - C_zero
DT_Prop_Const = (10^-1.15)*(Carrier_Frequency^-2); % Devices Propagation Constant - C_d
CT_Prop_Const_dB = pow2db(CT_Prop_Const); % Cellular Propagation Constant - C_zero in dB
DT_Prop_Const_dB = pow2db(DT_Prop_Const); % Devices Propagation Constant - C_d in dB

Cell_Radius = 500; % Cell Radius in meters
UE_Dist_Min = 50; % D2D group radius
D2D_Sep_Max = 0.1*Cell_Radius;
K=10;% number of D2D groups
N_DTs = input('Enter N value as 10 or 7') ;
N_Runs = 80;

%Total_SE_DPS = zeros(N_Runs,30,2);
SE_FST_Prev = log2(1+CT_SNR_Target);
SE_CL_Prev = log2(1+CT_SINR_min);
Final_SE_FST_AOS = zeros(1,N_Runs);

Assigned_DT_Power_AOS = zeros(N_Runs,N_DTs);

%%%% Threshold Reuse Distances based on Fixed SNR Target (FST) - used to Control the selection Criteria
d_DT2DT_Reuse_FST = 3;
R_DT2eNB_Reuse_FST = 3;

%%%% Definition of eNB position
eNB_x = 0;
eNB_y = 0;

for Iter = 1:N_Runs
%%% Placement of Users
    D2D_Link_list = LTE_UE_uniform_distribution(eNB_x,eNB_y,Cell_Radius,D2D_Sep_Max, N_DTs);
%%% D2D Pairs Link Gain
    for jj=1:N_DTs
        DT_Pair_gain(jj) = LTE_channel_model_indoor_hotspot_NLOS(D2D_Link_list(jj,1),D2D_Link_list(jj,2),D2D_Link_list(jj,3),D2D_Link_list(jj,4));
        P_DT_FST_AOS(jj) = (DT_SNR_Target*Noise_Total_Watts)/DT_Pair_gain(jj); % Transmit power set based on target SNR at the DT receiver
        DT_rx_Power_FST_AOS(jj) = P_DT_FST_AOS(jj) * DT_Pair_gain(jj); % DT received power
    end
   
%%%% Separation Distances between all - DTs
    for j=1:N_DTs
        for k=1:N_DTs
            Distance(j,k) = pdist([D2D_Link_list(j,1) D2D_Link_list(j,2);D2D_Link_list(k,3) D2D_Link_list(k,4)]);
        end
    end
%%%% Separation Distances between all DT_TXs and eNB
    for k=1:N_DTs
        DT2eNB_Distance(k) = pdist([D2D_Link_list(k,1) D2D_Link_list(k,2);eNB_x eNB_y]);
        DT2eNB_gain(k) = LTE_channel_model_urban_micro_NLOS_Thesis(D2D_Link_list(k,1),D2D_Link_list(k,2),eNB_x,eNB_y);
    end
%%%% Placement of the CT user in the cell
    loCT = UE_Dist_Min + (Cell_Radius - UE_Dist_Min)*sqrt(rand(1,1));
    theta= 2*pi*rand(1,1); % Generate the random angle Theta of the points
    CT_x_tx = loCT*cos(theta) + eNB_x ;
    CT_y_tx = loCT*sin(theta) + eNB_y ;
    R_CT_eNB = pdist([CT_x_tx CT_y_tx;eNB_x eNB_y]);
    R_CT2DT_Reuse_FST = (CT_SNR_Target*DT_SINR_min/(DT_SNR_Target-DT_SINR_min))^(1/DT_Exp)*R_CT_eNB;
    R_DT2eNB_Reuse_CL = ((CT_Prop_Const*CT_SINR_min*DT_P_Tx_Max)/((CT_Prop_Const*CT_P_Tx_Max*(R_CT_eNB^-CT_Exp))-(Noise_Total_Watts*CT_SINR_min)))^(1/CT_Exp);
%%%% Separation Distances between all - CT_TXs & DT_RXs
    for k=1:N_DTs
        CT2DT_Distance(k) = pdist([CT_x_tx CT_y_tx;D2D_Link_list(k,3) D2D_Link_list(k,4)]);
        CT2DT_gain(k) = LTE_channel_model_indoor_hotspot_NLOS(D2D_Link_list(k,3),D2D_Link_list(k,4),CT_x_tx,CT_y_tx);
    end
%%%% Power allocation according to SNR Target
    CT2eNB_gain = LTE_channel_model_urban_micro_NLOS_Thesis(CT_x_tx,CT_y_tx,eNB_x,eNB_y);
    Ch_Power_FST_AOS = CT_SNR_Target*Noise_Total_Watts/CT2eNB_gain; %%%%%%%%%%%%%%%% 
    DT_DT_Power_FST_AOS = Ch_Power_FST_AOS * DT_Pair_gain; 
    SINRC=abs((Ch_Power_FST_AOS*DT_Pair_gain)./((sum(loCT)*DT_DT_Power_FST_AOS*CT2eNB_gain)+Noise_dBm));
    SINRD1=abs((Ch_Power_FST_AOS*DT_Pair_gain)./((DT_DT_Power_FST_AOS*CT2eNB_gain)+Noise_dBm));
    SINRD2=abs((Ch_Power_FST_AOS*DT_Pair_gain)./((DT_DT_Power_FST_AOS*CT2eNB_gain)+Noise_dBm+(Ch_Power_FST_AOS*DT_Pair_gain)));
    Csum=log(1+SINRD1)+log(1+SINRD2);
    esum=Csum/((Pd+21));
    if abs(SINRC)>=0.000001
        [Cfinalsum(1,Iter),cost] = munkres(Csum);
    else
    end
end
Cfinalsum(Cfinalsum==0)=[];
eps=0.76;
P11=eps*Pd;
P21=(1-eps)*Pd;
for ii=1:N_DTs-1
if DT_Pair_gain(1,ii)>DT_Pair_gain(1,ii+1)
    Cfinalsum(1,ii)=P21;
else
    Cfinalsum(1,ii)=P11;
    epos1=((0.001*((DT_DT_Power_FST_AOS*CT2eNB_gain)+Noise_dBm))./(Pd*DT_Pair_gain));
    epos2=((Pd*DT_Pair_gain)-(0.001*(DT_DT_Power_FST_AOS*CT2eNB_gain)+Noise_dBm))./(Pd*DT_Pair_gain*(1+0.001));
    P1=epos1*Pd;
    P2=(1-epos2)*Pd;
    e1opt=(log(1+epos1.*SINRD1)+log(1+(1-epos1).*SINRD2))/((Pd+21))  ;
    e2opt=(log(1+epos2.*SINRD1)+log(1+(1-epos2).*SINRD2))/((Pd+21))  ;
    ejsum=max(e1opt,e2opt);
    maxeff=10*((((1/2)*log(1+((1/2).*SINRD1))*1e3)+((1/2)*log(1+((1/2).*SINRD2))*1e3))/((Pd+21)));
    throuh=18*log10(1+((SINRD1+SINRD2+SINRC)./3000));
end
end

if N_DTs==10
D2Dgr=[7:1:15];
figure,
plot(D2Dgr,(sort(abs(maxeff(1,1:9)),'ascend')),'r-o');
ylabel('Total Energy Efficiency (bps/Hz/W)')
xlabel('Number D2D user groups')
title('Total Energy Efficiency with D2D groups');


D2Dtr=[0.12:0.02:0.22];
figure,
plot(D2Dtr,(sort(abs(maxeff(1,1:6)),'descend')),'b-o');
ylabel('Total Energy Efficiency (bps/Hz/W)')
xlabel('D2D Transmit Power(W)')
title('Relation ship between Total Energy Efficiency of D2D users and D2D transmit power');

D2Dgro=[1:1:10];
figure,
plot(D2Dgro,sort(abs(throuh(1,1:10)),'ascend'),'b->');
ylabel('Throughput(Mbps)')
xlabel('No. of D2D users')
title('Throughput of proposed algorithm');

D2Dgro=[0:5:45];
figure,
plot(D2Dgro,sort((throuh(1,1:10)./abs(max(throuh(:)))),'ascend'),'g->');
ylabel('Cummulaive distribution function')
xlabel('No. of D2D users')
title('Comparison of CDF throughput');
ylim([0 1])
rng(s); 
else
D2Dgr=[7:1:12];
figure,
plot(D2Dgr,(sort(abs(maxeff(1,1:6)),'ascend')),'r-o');
ylabel('Total Energy Efficiency (bps/Hz/W)')
xlabel('Number D2D user groups')
title('Total Energy Efficiency with D2D groups');


D2Dtr=[0.12:0.02:0.22];
figure,
plot(D2Dtr,(sort(abs(maxeff(1,1:6)),'descend')),'b-o');
ylabel('Total Energy Efficiency (bps/Hz/W)')
xlabel('D2D Transmit Power(W)')
title('Relation ship between Total Energy Efficiency of D2D users and D2D transmit power');

D2Dgro=[1:1:6];
figure,
plot(D2Dgro,sort(abs(throuh(1,1:6)),'ascend'),'b->');
ylabel('Throughput(Mbps)')
xlabel('No. of D2D users')
title('Throughput of proposed algorithm');

D2Dgro=[0:5:25];
figure,
plot(D2Dgro,sort((throuh(1,1:6)./abs(max(throuh(:)))),'ascend'),'g->');
ylabel('Cummulaive distribution function')
xlabel('No. of D2D users')
title('Comparison of CDF throughput');
ylim([0 1])
rng(s);
end

