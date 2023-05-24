clc
close all;
clear all;
%%% INITIALIZATION %%
xm=100;
ym=100;
sinkx=50;
sinky=125;              
numNodes = 100; % number of nodes
ETX=50*0.000000001;  %tx energy
ERX=50*0.000000001;  %rx energy
Efs=10*0.000000000001;  %free space loss
Emp=0.0013*0.000000000001;   %multipath loss
L=3200;
%Data Aggregation Energy
EDA=5*0.000000001;  %compression energy
a=1.5;   %fraction of energy enhancment of advance nodes
rmax=2500 ; %maximum number of rounds
Eelec=50*0.000000001;
Eini=1;
do=sqrt(Efs/Emp);
 % Aux Parameters
t = 0;
currentError = 1;

sample1 = [ ( 75.2 - 32 ).*rand( 1, 20 ) + 32, ...
           ( 60 - 20 ).*rand( 1,20 ) + 20, ...
           ( 15 -( -3.4 ) ).*rand( 1, 20 ) + ( -3.4 ),... 
           (50 - 32 ).*rand( 1, 20 ) + (-15)....
           ( 35.2 - 22 ).*rand( 1, 20 ) + (-2.4)... 
          ]';
sample2 = [ ( 19.2 - 1 ).*rand( 1,20 ) + 1, ...
           ( -10 - ( -42.5 ) ).*rand( 1,20 ) + ( -42.5 ), ...
           ( 20.4 - ( -5 ) ).*rand( 1,20 ) + ( -5 ),...
           ( 85.2 - 12 ).*rand( 1, 20 ) + (-3)....
           ( 95.2 - 22 ).*rand( 1, 20 ) + (-1.4)...
          ]';
Inputs = [ sample1, sample2 ];
Bs=[sinkx,sinky];
[nInputs, nInputSamples ] = size( Inputs );
    % Number of Clusters
    nC = 5;
    % fuzzification parameter
    m = 2; 
    % Generate the random association of inputs to the clusters, values [ 0, 1 ]
    n3 = rand( nInputs, nC );
    nTotal = sum( n3, 2 );
    RandomAssociationValues = ( n3./nTotal );
    % Initial Membership Matrix
    U = RandomAssociationValues;
    % Cluster's centroids
    Cf = zeros( nC, nInputSamples );
    Imax=10;
    P=0.85;
    % Aux Parameters
    t = 0;
    currentError = 1;
    error=1e-9;
    imp=0.001;
    % plot parameters
    minAxis = min( Inputs );
    maxAxis = max( Inputs );
     d1=0.765*xm/2;  %distance between cluster head and base station
K=5;
 numClusters =sqrt(1.262*n3/2*pi)*xm/d1; %optimal no. of clusters
    while currentError > error 
        U0 = U;

% Calculate the cluster's centroids
for ia=1:100
        for i = 1 : 1 : nC
            for j = 1 : 1 : nInputSamples
                Cf( i, j ) = ( sum ( Inputs( :, j ) .* ( U( :, i ).^m ) ) )/( sum ( U( :, i ).^m ) );
            end
        end
        % calculate dissimilarly between the inputs and centroids using
        % euclidean distance   
        distanceFromCluster = zeros( nInputs, nC );
   for k = 1 :1: nC       
         distance = sum( ( ( Inputs - Cf( k, : ) ).^2 ), 2 );
            distanceFromCluster( : , k) = sqrt( distance );
   end
        % update membership matrix values
              den = sum( ( ( 1./distanceFromCluster ).^( 1/( m-1 ) ) ), 2 );
        for z = 1 : 1 : nC
            num = ( ( 1./distanceFromCluster( :, z ) ).^( 1/( m-1 ) ) ) ./ den;
            U( :, z ) = num';
        end 
       cluthre= (numNodes*P)/ nC ;
      for i1=1:Cf
         B = sort(K);
         min_cluster=min(B);
      end
      if min_cluster > cluthre
            break;
       else
     for i11=1:Cf
         for j1=1:numNodes
       distanceFromCentroid = zeros(numNodes, nInputSamples );
              distance1(:,:) = sum( ( ( [Cf(1,1),Cf( 1, 2)]- Inputs).^2 ), 2 );
              distanceFromCentroid = (sqrt( distance1));
               B1 = sort( distanceFromCentroid,'descend');     
         end
     end    
       for i2=1:Cf
           Idx = knnsearch(Inputs,  Cf( 5, : ));

   Inputs1 =cluthre*Inputs;

 [idx, dist] = knnsearch(Inputs1,  Cf( 5, : ),'dist','cityblock','k',2);           
       end
        
 currentError = ( sqrt( ( U - U0 ).^2 ) );
 t = t + 1;
      end
end
    end  
Econ=100000;
dBs=50;
dF=33;
 for i=1:K
      for j=1 numClusters(i)
        Er= Eini-Econ;
          for i=1:100
             for j=1:K
                  if dBs > dF
                      dFCH=min( ( ( Inputs - Cf(j) ).^2 ), 2 );
                  elseif dBs < do
%                       dFCH =((Inputs - BS ).^2 ), 2 );
                  end
             end
          end
                   if j>1
                     dBCH=min ( ( ( Inputs - Cf(j-1) ).^2 ), 2 );
                   else
                     dBCH=0;
                   end
                     ACD=(min( dFCH, dBCH))/(max(dFCH, dBCH));
                     F= ((Er)./(dFCH + dBCH))+ACD(100,2);
      end
             Tb=1./F;
             ETH=0.1;
             n=100;
             Echp=50*0.000000001;
             En=2*0.000000001;
             Rchs=(n-1)*(ETH/Echp);
             Erth=Eini-ETH;
             Rn=Er/En;
             T=intersect( Rchs,Rn);
            if (Econ+T*Eini)<=Eini
                 ETH= Econ+(Eini-T)*Eini;
            else
                  ETH=Econ+Er;
            end
          end   

xm=100;      %diameters of sensor network
ym=100;
sink.x=100;  %distance of base station from the network
sink.y=100;
n = 100;         %no of nodes
p=0.1;          %probibilty of a node to become cluster head
Eo=0.5;          %energy supplied to each node
ETX=50*0.000000001;     %transmiter energy per node
ERX=50*0.000000001;        %reciever energy per mode
Efs=10*0.0000000000001;     %amplification energy when d is less than d0
Emp=0.0013*0.000000000001;      %amplification energy  when d is greater than d0
%Data Aggregation Energy
EDA=5*0.000000001;
rmax=5000;           %no of rounds
do=sqrt(Efs/Emp);
h=100;
sv=0;                                  %%%%%%previously Sensed value S(v)
tempi=50;
tempf=200;
min_dis_cluster=0;
m=0;
a=0;

for i=1:1:n
    S3(i).xd=rand(1,1)*xm;         %it will distribute the nodes in 1 dimension in x axis randomly.
    S3(i).yd=rand(1,1)*ym;           %it will distribute the nodes in 1 dimension in y axis randoml.
    S3(i).G=0;                        % as the no of node that have been cluster head is zero 0
    S3(i).E=Eo;%%*(1+rand*a);                
    S3(i).type='N';
end

S3(n+1).xd=sink.x;
S3(n+1).yd=sink.y;
Emp=0.0013*0.000000000001;

%counter for CHs
countCHs3=0;
%counter for CHs per round
cluster3=1;

allive3=200;
flag_first_dead3=0;
packets_TO_BS3=0;
packets_TO_CH3=0;
s3=0;
th=0.00000000000001;
for r=0:1:rmax
    cv = tempi + (tempf-tempi).*rand(1,1);  %%%%%%Current sensing value C(v)
    %Operation for epoch
    if(mod(r, round(1/p) )==0)
        for i=1:1:n
            S3(i).G=0;
            S3(i).cl=0;
        end
    end   
    
    %Number of dead nodes
    dead3=0;
    %Number of dead Advanced Nodes
    dead_a3=0;
    %Number of dead Normal Nodes
    dead_n3=0;
       
    for i=1:1:n
        %checking if there is a dead node
        if (S3(i).E<=0)
            dead3=dead3+1;          
        
        end        
        if (S3(i).E<=th && S3(i).E>0)
            s3=s3+1;
        end
        
        if S3(i).E>0
            S3(i).type='N';
        end     
        
    end
    
    DEAD3(r+1)=dead3;
    STATISTICS.ALLIVE3(r+1)=allive3-dead3;
    %When the first node dies
    if (dead3==1)
        if(flag_first_dead3==0)
            first_dead3=r
            flag_first_dead3=1
        end
    end
    
    countCHs3=0;
    cluster3=1;
    for i=1:1:n
        if(S3(i).E>=th && s3<10)  %Checking threshold and nmbr of sleep nodes
            temp_rand=rand;
            if ( (S3(i).G)<=0)                
                %Election of Cluster Heads
                if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
                    countCHs3=countCHs3+1;
                    S3(i).type='C';
                    S3(i).G=round(1/p)-1;
                    C(cluster3).xd=S3(i).xd;
                    C(cluster3).yd=S3(i).yd;
                
                    distance=sqrt( (S3(i).xd-(S3(n+1).xd) )^2 + (S3(i).yd-(S3(n+1).yd) )^2 );
                    C(cluster3).distance=distance;
                    C(cluster3).id=i;
                    X(cluster3)=S3(i).xd;
                    Y(cluster3)=S3(i).yd;
                    cluster3=cluster3+1;
                    
                    %Calculation of Energy dissipated
                    distance;
                    
                    if (cv >= h)
                                              
                        if (distance>do)
                            S3(i).E=S3(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
                            packets_TO_BS3=packets_TO_BS3+1;
                            countCHs3=countCHs3+1;
                                                   
                        end
                        if (distance<=do)
                            S3(i).E=S3(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
                            packets_TO_BS3=packets_TO_BS3+1;
                            countCHs3=countCHs3+1;
                        end
                        
                         
                    end
                end
                
            end
        end
    end
    COUNTCHS3(r+1)=countCHs3;
    

    
    %Election of Associated Cluster Head for Normal Nodes
    for i=1:1:n
        if ( S3(i).type=='N' && S3(i).E>=th && s3<10)
            if(cluster3-1>=1)
                min_dis=sqrt( (S3(i).xd-S3(n+1).xd)^2 + (S3(i).yd-S3(n+1).yd)^2 );
                min_dis_cluster3=1;
                for c=1:1:cluster3-1
                    temp=min(min_dis,sqrt( (S3(i).xd-C(c).xd)^2 + (S3(i).yd-C(c).yd)^2 ) );
                    if ( temp<min_dis )
                        min_dis=temp;
                        min_dis_cluster3=c;
                    end
                end
                
                %Energy dissipated by associated Cluster Head
                if (cv >= h)
                    
                    
                    min_dis;
                    if (min_dis>do)
                        S3(i).E=S3(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis));
                    end
                    if (min_dis<=do)
                        S3(i).E=S3(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
                    end
                    
                    if(min_dis>0)
                        S3(C(min_dis_cluster3).id).E = S3(C(min_dis_cluster3).id).E- ( (ERX + EDA)*4000 );
%                         PACKETS_TO_CH(r+1)=n-dead-cluster+1;
                    end
                   
                end           
                              
            end
        end
        
         %   Activation of sleep nodes when their nmbr exceed 9
        if(s3>=10 && S3(i).E>=0)
            if(cluster-1>=1)
                min_dis=Inf;
                min_dis_cluster3=0;
                for c=1:1:cluster-1
                    temp=min(min_dis,sqrt( (S3(i).xd-C(c).xd)^2 + (S3(i).yd-C(c).yd)^2 ) );
                    if ( temp<min_dis )
                        min_dis=temp;
                        min_dis_cluster3=c;
                    end
                end               
                
                min_dis;
                if (cv >= h)
                if (min_dis>do)
                    S3(i).E=S3(i).E- ( ETX*(4000) + Emp*4000*( min_dis *min_dis * min_dis * min_dis));
                end
                if (min_dis<=do)
                    S3(i).E=S3(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
                end
                S3(C(min_dis_cluster3).id).E =S3(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 );
                packets_TO_CH3=packets_TO_CH3+1;
                S3(i).min_dis=min_dis;
                S3(i).min_dis_cluster3=min_dis_cluster3;
                end
            end
        end
        
    end   
    
   PACKETS_TO_BS3(r+1)=packets_TO_BS3;    

end
figure;
r=0:rmax;
plot(r,STATISTICS.ALLIVE3(r+1),'-b');
xlim([1500 5000])
xlabel('No. of Rounds (r)');
ylabel('No. of Nodes Allive')

xm=100;
ym=100;
sink.x=0.5*xm;  %location of sink on x-axis
sink.y=0.5*ym;  %location of sink on y-axis
r1=2500; 
P=0.1;  %probability of cluster heads
Eo=0.5;%initial energy
Echeck=Eo;
ETX=50*0.000000001;  %tx energy
ERX=50*0.000000001;  %rx energy
Efs=10*0.000000000001;  %free space loss
Emp=0.0013*0.000000000001;   %multipath loss
%Data Aggregation Energy
EDA=5*0.000000001;  %compression energy
a=1.5;   %fraction of energy enhancment of advance nodes
rmax=2500 ; %maximum number of rounds
do=sqrt(Efs/Emp);  %distance do is measured
Et=0;  %variable just use below 
m=0.5;
mo=0.4;
b=3;%fraction of extra energy in super nodes
normal=r1*(1-m);
advance=r1*m*(1-mo);
super=r1*m*mo;
m7=0;
c=0.02;
mony=0;
A=0;
for i=1:1:r1
    S1(i).xd=rand(1,1)*xm;

    S4(i).xd=S1(i).xd;
    XR4(i)=S4(i).xd;

    S1(i).yd=rand(1,1)*ym;

    S4(i).yd=S1(i).yd;
    YR4(i)=S4(i).yd;
    S4(i).G=0;

    talha=rand*a;
    S1(i).E=Eo*(1+talha);
    S2(i).E=S1(i).E;
    S3(i).E=S1(i).E;
    S4(i).E=S3(i).E;
    E3(i)= S3(i).E;
    E4(i)= S4(i).E;
    Et=Et+E3(i);
    E(i)= S1(i).E;
    A=A+talha;             
    S1(i).A=talha;
    S2(i).A=talha;
    if (S2(i).E>Echeck)
        mony=mony+1;
    end
    
    Et=Et+E3(i);

    %initially there are no cluster heads only nodes

    S4(i).type='N';
end
figure,
plot(sort(E3,'ascend'),'rd-');
hold on
xlabel('rounds');
ylabel('Energy consumption for nodes');
title('Energy consumption for nodes');
NT=4; %%%%%%%%%%transmitting antennas
N=25; %%%%%%%%%%%%%%%%Resource blocks
BW=180; %%%total bandwidth
Subc=12; %%%%%%%%%%%%%%%total subcarriers
t=0.5; %%%%%%%%time in ms
UEs=40:5:90; %%%%%%%%%%no.of UEs
Tp=30; %%%Transmitting power
Power_Symbol=5;
Kt = 4; %Number of base stations (BSs)
Noise_dBm =	-174;
InfDAC_Flag=1;
Symbols = [ -1-1i,-1+1i,1-1i,1+1i ];  
for m = 1:length(UEs)
D = zeros(Kt*NT,Kt*NT,UEs(m));
for k = 1:UEs(m)
    D((k-1)*NT+1:k*NT,(k-1)*NT+1:k*NT,k) = eye(NT);
end
Bits = randi([0, 1], UEs(m), 1000);
index_s = bi2de(Bits, 'left-msb') + 1; % convert binary to decimal
s =(index_s)';  

%Combined channel matrix will be (Kr x Kt*Nt). This matrix gives the
channelVariances = 1/2+1/2*kron(eye(UEs(m)),ones(1,NT));

%User weights for (unweighted) sum rate computation
weights = ones(UEs(m),1);

%Pre-generation of Rayleigh fading channel realizations (unit variance)
Hall = (randn(UEs(m),NT*UEs(m),length(UEs))+1i*randn(UEs(m),NT*UEs(m),length(UEs)))/sqrt(2);
%Iteration over channel realizations

%Generate channel matrix for m:th realization
H = sqrt(channelVariances) .* Hall(:,:,m);
[W, ~, ~] = ZF_Precoder(H, s, Power_Symbol, Tp,InfDAC_Flag,UEs, Kt);
for k=1:UEs(m)-1
Spatcor(k)=(H(k,:)*(H(k+1,:)'))/(( norm(H(k,:), 'fro'))*( norm(H(k+1,:), 'fro')));
alphamin=min(Spatcor(:));
alphamax=max(Spatcor(:));
alpha_n=(alphamin+alphamax)/2;
muldt=H(k,:).*W(k,:)';
g(k)=(norm(muldt, 'fro').^2)/(Noise_dBm.^2);
d1=(sum(g(k)./(1+g(k).*(alpha_n-(1/g(k)))))).*(1/log(2));
d2=((sum(g(k)).^2./((1+g(k).*(alpha_n-(1/g(k)))))).^2).*(1/log(2));
d3=sum(alpha_n-(1/g(k)))+Tp;
B=BW/N;
devf=B.*(((d1.*d3)-(NT.*sum(log(1+g(k)*(alpha_n-(1/g(k))./(d3.^2)))))))*alpha_n;
if devf<=0
   alpha_n=alphamax;
else
   alpha_n=alphamin;
end
p(k)=(alpha_n-(1/g(k)));
gamman(k)=(p(k)*g(k));
end

pow_all=(max((sum(B.*log(1+gamman)))./(sum(p)+Tp)))/N;
energyef=(B.*log(1+gamman))./(sum(p)+Tp);
specteff=log(1+gamman);
datarate=B.*log(1+gamman);
throughpu=(1000-100)*datarate*Subc*8;
finalpow_all(1,m)= abs(pow_all(:));
finalenergyef(1,m)=abs( max(energyef(:)));
finalspecteff(1,m)= abs(max(specteff(:)));
finalthroughpu(1,m)=abs(max(throughpu));
Jain=((sum(throughpu)./datarate).^2)./(UEs(m)*(sum(throughpu./datarate).^2));
Jainin(1,m)=abs(max(Jain(:)));
end
UEs=60;
pow=20:2:40;
for m = 1:length(pow)
D = zeros(Kt*NT,Kt*NT,UEs);
for k = 1:UEs
    D((k-1)*NT+1:k*NT,(k-1)*NT+1:k*NT,k) = eye(NT);
end
Bits = randi([0, 1], UEs, 1000);
index_s = bi2de(Bits, 'left-msb') + 1; % convert binary to decimal
s =(index_s)';  

%Combined channel matrix will be (Kr x Kt*Nt). This matrix gives the
channelVariances = 1/2+1/2*kron(eye(UEs),ones(1,NT));

%User weights for (unweighted) sum rate computation
weights = ones(UEs,1);

%Pre-generation of Rayleigh fading channel realizations (unit variance)
Hall = (randn(UEs,NT*UEs,length(UEs))+1i*randn(UEs,NT*UEs,length(UEs)))/sqrt(2);
%Iteration over channel realizations

%Generate channel matrix for m:th realization
H = sqrt(channelVariances) .* Hall(:,:,1);
[W, ~, ~] = ZF_Precoder(H, s, Power_Symbol, Tp,InfDAC_Flag,UEs, Kt);
for k=1:UEs-1
Spatcor(k)=(H(k,:)*(H(k+1,:)'))/(( norm(H(k,:), 'fro'))*( norm(H(k+1,:), 'fro')));
alphamin=min(Spatcor(:));
alphamax=max(Spatcor(:));
alpha_n=(alphamin+alphamax)/2;
muldt=H(k,:).*W(k,:)';
g(k)=(norm(muldt, 'fro').^2)/(Noise_dBm.^2);
d1=(sum(g(k)./(1+g(k).*(alpha_n-(1/g(k)))))).*(1/log(2));
d2=((sum(g(k)).^2./((1+g(k).*(alpha_n-(1/g(k)))))).^2).*(1/log(2));
d3=sum(alpha_n-(1/g(k)))+Tp;
B=BW/N;
devf=B.*(((d1.*d3)-(NT.*sum(log(1+g(k)*(alpha_n-(1/g(k))./(d3.^2)))))))*alpha_n;
if devf<=0
   alpha_n=alphamax;
else
   alpha_n=alphamin;
end
p(k)=(alpha_n-(1/g(k)));
gamman(k)=(p(k)*g(k));
end

pow_all1=(max((sum(B.*log(1+gamman)))./(sum(p)+pow(m))))/N;
energyef1=(B.*log(1+gamman))./(sum(p)+pow(m));
specteff1=log(1+gamman);
throghpu1=B.*log(1+gamman);
finalpow_all1(1,m)= abs(pow_all1(:));
finalenergyef1(1,m)=abs( max(energyef1(:)));
finalspecteff1(1,m)= abs(max(specteff1(:)));
finalthroughpu1(1,m)=abs(max(throghpu1));
end

figure,
plot(pow,sort(finalenergyef1*1e4,'descend'),'r-o');
ylabel('Energy Efficiency')
xlabel('Maximum Transmit power(W)')
title('Energy effeciency vs. HAP transmit power');

figure,
plot(pow,sort((finalthroughpu1)/7*1e-5,'ascend'),'r-o');
ylabel('Throughput')
xlabel('Maximum Transmit power(W)')
title('Throughput vs. HAP transmit power');