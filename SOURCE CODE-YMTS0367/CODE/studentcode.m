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
     %--- visualize data when it's possible
        fig = figure( 'Name', 'Fuzzy c Means Ilustration', 'NumberTitle', 'off' );
       
            plot2DData( fig, 1, Bs,Inputs, 0, minAxis, maxAxis, 'Input Data' );
%  plot2DData( fig, 1, Bs,Inputs, 0, 'Input Data' );
%   axis([1 xm 1 ym]);
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
 
 plot2DData( fig, 2, Bs,Inputs, Cf, minAxis, maxAxis, [ 'Centroids at ', 'Epoch: ', num2str( t ) ] );;
      end
end
    end
  plotClusterData( Inputs, U, nC );   
% %   %CHSRA:

%backoff timer mechanism
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
rmax=2500;           %no of rounds
do=sqrt(Efs/Emp);
h=100;
sv=0;                                  %%%%%%previously Sensed value S(v)
tempi=50;
tempf=200;
min_dis_cluster=0;
%Values for Hetereogeneity
%Percentage of nodes than are advanced
m=0;
%\alpha
a=0;

for i=1:1:n
    S3(i).xd=rand(1,1)*xm;         %it will distribute the nodes in 1 dimension in x axis randomly.
    S3(i).yd=rand(1,1)*ym;           %it will distribute the nodes in 1 dimension in y axis randoml.
    S3(i).G=0;                        % as the no of node that have been cluster head is zero 0
    S3(i).E=Eo;%%*(1+rand*a);                
    %initially there are no cluster heads only nodes
    S3(i).type='N';
end

S3(n+1).xd=sink.x;
S3(n+1).yd=sink.y;
Emp=0.0013*0.000000000001;

%counter for CHs
countCHs3=0;
%counter for CHs per round
cluster3=1;

allive3=n;
flag_first_dead3=0;
packets_TO_BS3=0;
packets_TO_CH3=0;
s3=0;
th=0.00000000000001;
for r=0:1:rmax
    r
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
        
        %         check for sleep nodes
        if (S3(i).E<=th && S3(i).E>0)
            s3=s3+1;
        end
        
        if S3(i).E>0
            S3(i).type='N';
        end
        
        
    end
    
    DEAD3(r+1)=dead3;

%     ALLIVE3(r+1)=allive3-dead3;
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
    % for i=1:n
    %          if S(i).E==0
    %              AliveNode(rmax)= AliveNode(rmax)-1;
    %          end
    % end
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
                        %test = cv-sv;
                        %if (test >= s)
                        
                        if (distance>do)
                            S3(i).E=S3(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
                            packets_TO_BS3=packets_TO_BS3+1;
                            countCHs3=countCHs3+1;
%                             PACKETS_TO_BS(i).j=packets_TO_BS;
                            
                        end
                        if (distance<=do)
                            S3(i).E=S3(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
                            packets_TO_BS3=packets_TO_BS3+1;
                            countCHs3=countCHs3+1;
%                             PACKETS_TO_BS(i).j=packets_TO_BS;
                        end
                        
                        %end
                        
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
xlabel('No. of Rounds (r)');
ylabel('No. of Nodes Allive')

%---------------------------

xm=100;
ym=100;
sink.x=0.5*xm;  %location of sink on x-axis
sink.y=0.5*ym;  %location of sink on y-axis
r1=2500; 
P=0.1;  %probability of cluster heads
Eo=0.5;%initial energy
%
Echeck=Eo;
%
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

d1=0.765*xm/2;  %distance between cluster head and base station
K=sqrt(0.5*r1*do/pi)*xm/d1^2; %optimal no. of cluster heads
d2=xm/sqrt(2*pi*K);  %distance between cluster members and cluster head
Er=4000*(2*r1*ETX+r1*EDA+K*Emp*d1^4+r1*Efs*d2^2);  %energy desipated in a round

S4(r1+1).xd=sink.x; %sink is a n+1 node, x-axis postion of a node
S4(r1+1).yd=sink.y; %sink is a n+1 node, y-axis postion of a node

countCHs4=0;  %variable, counts the cluster head

cluster4=1;  %cluster is initialized as 1

flag_first_dead4=0; %flag tells the first node dead

flag_teenth_dead4=0;  %flag tells the 10th node dead

flag_all_dead4=0;  %flag tells all nodes dead

dead4=0;  %dead nodes count initialized to 0

first_dead4=0;

teenth_dead4=0;

all_dead4=0;

allive4=n;
%counter for bit transmitted to Bases Station and to Cluster Heads

packets_TO_BS4=0;

packets_TO_CH4=0;


for i=1:1:super
  
    S4(i).E=Eo*(1+b);
    

end
talha1=super+advance;
for i=super:1:talha1
    
    S4(i).E=Eo*(1+a);
    
end
for i=talha1:1:r1
    
    S4(i).E=Eo;
    
end

for r=0:1:rmax     
    r
  if(mod(r, round(1/P) )==0)
    for i=1:1:r1
        S4(i).G=0;
        S4(i).cl=0;
    end
  end
  
 Ea=Et*(1-r/rmax)/r1;
 dead4=0;
for i=1:1:r1
   
    if (S4(i).E<=0)
        dead4=dead4+1; 
        if (dead4==1)
           if(flag_first_dead4==0)
              first_dead4=r;
              flag_first_dead4=1;
           end
        end
        if(dead4==0.1*n)
           if(flag_teenth_dead4==0)
              teenth_dead4=r;
              flag_teenth_dead4=1;
           end
        end
        if(dead4==n)
           if(flag_all_dead4==0)
              all_dead4=r;
              flag_all_dead4=1;
           end
        end
    end
    if S4(i).E>0
        S4(i).type='N';
    end
end

STATISTICS.DEAD4(r+1)=dead4;
STATISTICS.ALLIVE4(r+1)=allive4-dead4;
countCHs4=0;
cluster4=1;
for i=1:1:r1
    
 
 if Ea>0
     if (S4(i).E<= Eo)
         p(i)=P*E(i)/(1+m*(a+mo*b))*Ea;
     end
         
     else if (S4(i).E<=Eo*(1+a))
         p(i)=P*(1+a)*E(i)/(1+m*(a+mo*b))*Ea;
         end
    if (S4(i).E<=Eo*(1+b))
         p(i)=P*(1+b)*E(i)/(1+m*(a+mo*b))*Ea;
     end
     

 if(S4(i).E>0)
   temp_rand=rand;     
   if ( (S4(i).G)<=0)  
        if((temp_rand<= (p(i)*E(i)*K/(1-p(i)*r*mod(1,p(i))*Ea))))
            countCHs4=countCHs4+1;
            packets_TO_BS4=packets_TO_BS4+1;
            PACKETS_TO_BS4(r+1)=packets_TO_BS4;
             S4(i).type='C';
            S4(i).G=round(1/p(i))-1;
            C4(cluster4).xd=S4(i).xd;
            C4(cluster4).yd=S4(i).yd;
           distance=sqrt( (S4(i).xd-(S4(r1+1).xd) )^2 + (S4(i).yd-(S4(r1+1).yd) )^2 );
            C4(cluster4).distance=distance;
            C4(cluster4).id=i;
            X4(cluster4)=S4(i).xd;
            Y4(cluster4)=S4(i).yd;
            cluster4=cluster4+1;
           distance;
            if (distance>do)
                S4(i).E=S4(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                S4(i).E=S4(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
        end     
    
   end
   
 end 
 end
 end
STATISTICS.COUNTCHS4(r+1)=countCHs4;


for i=1:1:r1
   if ( S4(i).type=='N' && S4(i).E>0 )
     if(cluster4-1>=1)
       min_dis=sqrt( (S4(i).xd-S4(r1+1).xd)^2 + (S4(i).yd-S4(r1+1).yd)^2 );
       min_dis_cluster=0;
       for c=1:1:cluster4-1
           temp=min(min_dis,sqrt( (S4(i).xd-C4(c).xd)^2 + (S4(i).yd-C4(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       
       
       if(min_dis_cluster~=0)    
            min_dis;
            if (min_dis>do) 
                S4(i).E=S4(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S4(i).E=S4(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
      
            S4(C4(min_dis_cluster).id).E = S4(C4(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
            packets_TO_CH4=packets_TO_CH4+1;
       else 
            min_dis;
            if (min_dis>do)
                S4(i).E=S4(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S4(i).E=S4(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
            packets_TO_BS4=packets_TO_BS4+1;
            
       end
        S4(i).min_dis=min_dis;
       S4(i).min_dis_cluster=min_dis_cluster;
   else
            min_dis=sqrt( (S4(i).xd-S4(r1+1).xd)^2 + (S4(i).yd-S4(r1+1).yd)^2 );
            if (min_dis>do)
                S4(i).E=S4(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S4(i).E=S4(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
            packets_TO_BS4=packets_TO_BS4+1;
   end
  end
end
STATISTICS.PACKETS_TO_CH4(r+1)=packets_TO_CH4;
STATISTICS.PACKETS_TO_BS4(r+1)=packets_TO_BS4;
end
% first_dead4
% teenth_dead4
% all_dead4
STATISTICS.DEAD4(r+1)
STATISTICS.ALLIVE4(r+1)
STATISTICS.PACKETS_TO_CH4(r+1)
STATISTICS.PACKETS_TO_BS4(r+1)
STATISTICS.COUNTCHS4(r+1)
r=0:2500;
figure,
plot(r,STATISTICS.PACKETS_TO_BS4,'-.m');
xlabel('x(rounds)');
ylabel('y(packets sent)');
title('Packets sent to the base station');
figure,
plot(sort(E3*100,'ascend'),'rd-');
%% 
 hold on
% plot(sort(totalenergy*100,'ascend'),'rd-');
% plot(sort(Totalnetwor_energy,'ascend'),'rd-');
% plot(1:2501,sort(E3*100,'ascend')); 
% plot(sort(E4/100,'ascend'),'b^');
xlabel('rounds');
ylabel('Energy consumption for nodes');
title('Energy consumption for nodes');


