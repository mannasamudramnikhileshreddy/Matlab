    clc
clear
close all
rng('default')
%%%%%%%%%%% Input parameters for formation of network

nodescou.nodeNo = 30;                                 %Number of nodes
nodescou.nodePosition(1,:) = [1 20 20];     %(Sender node fixed position)
nodescou.nodePosition(2,:) = [2 160 160];             %(Receiver node fixed position)
nodescou.NeighborsNo = 5;
nodescou.range = 50;                %Tolerance distance to became neighbor of one node (Euclidean distance based)
nodescou.atenuationFactor = 1.8;    %Atenuation factor in freespace - ranges from 1.8 to 4 due environment
nodescou.minEnergy = 0;            % Minimum threshold
nodescou.maxEnergy = 100;           % Max threshold
nodescou.energyconsumptionperCicle = 3*1e-14;
nodescou.energyrecoveryperCicle =  2*1e-14;
nodescou.energyfactor = 0.001;
iterationcounter=1;

% Node position sortition

for a = 3 : nodescou.nodeNo
    
   nodescou.nodeId = a; 
   nodespos.x = randi([1 140]); %Xpos sortition
   nodespos.y = randi([1 140]); %Ypos sortition
   nodescou.nodePosition(a,:) = [nodescou.nodeId nodespos.x nodespos.y]; %NodeID, X and Y position into nodePosition table
   
end

% Euclidean Distance calc from one node to all others

for i = 1 : nodescou.nodeNo
    for j = 1: nodescou.nodeNo
        nodespos.x1 = nodescou.nodePosition(i,2); 
        nodespos.x2 = nodescou.nodePosition(j,2); 
        nodespos.y1 = nodescou.nodePosition(i,3); 
        nodespos.y2 = nodescou.nodePosition(j,3);        
        nodescou.euclidiana(i,j) = sqrt(  (nodespos.x1 - nodespos.x2) ^2 + (nodespos.y1 - nodespos.y2)^2  ); 
        
    end
end
nodescou.weights = lt(nodescou.euclidiana,nodescou.range);
G=graph(nodescou.weights,'upper');
for a = 1 : height(G.Edges)
    nodespos.s = G.Edges.EndNodes(a,1);
    nodespos.t = G.Edges.EndNodes(a,2);
    nodespos.Z(a,:) = nodescou.euclidiana(nodespos.s,nodespos.t);
 end
G.Edges.Euclidiana = nodespos.Z(:,1);

[nodescou.nodePosition(:,4)] = nodescou.maxEnergy -(nodescou.maxEnergy-nodescou.minEnergy)*rand(nodescou.nodeNo,1);
for a = 1: length(nodescou.nodePosition(:,1))
   
    nodescou.nodePosition(a,5) = degree(G,nodescou.nodePosition(a,1));
    
end

[G.Edges.Pathloss] = (10*nodescou.atenuationFactor)*log10(G.Edges.Euclidiana);

for a = 1 : height(G.Edges)
	nodespos.Sourcenode = G.Edges.EndNodes(a,1);
	nodespos.Targetnode = G.Edges.EndNodes(a,2);
	G.Edges.SourcenodeXpos(a) = nodescou.nodePosition(nodespos.Sourcenode,2);
	G.Edges.SourcenodeYpos(a) = nodescou.nodePosition(nodespos.Sourcenode,3);
	G.Edges.TargetnodeXpos(a) = nodescou.nodePosition(nodespos.Targetnode,2);
	G.Edges.TargetnodeYpos(a) = nodescou.nodePosition(nodespos.Targetnode,3);
    G.Edges.ActiveEdge(a) = 1;
end

% Graph objects plot
figure('units','normalized','innerposition',[0 0 1 1],'MenuBar','none')
hold on
% Plot the points.
scatter(nodescou.nodePosition(:, 2), nodescou.nodePosition(:, 3), 'LineWidth', 3);
hold on;
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
grid on;
title('Nodes location');
hold off

%%%%%%%%%%%%%%SAR Calculation
dr=sort(nodescou.nodePosition(:,4),'ascend');
fmin=880*1e6;
fmax=915*1e6;
BW=fmax-fmin;
k=6*1e6;
sigma=0.02; %%%%%%%%%%%skull tissue conductivity(S/m)
row=1.1; %%%%%%%%%%density of human tissue(kg/m^3)-https://www.emf-portal.org/en/cms/page/home/technology/radio-frequency/mobile-communication
c=3*1e8;
D=1000;
f=(fmax+fmin)/2;
G=20;
No=70;
num=(sigma*((4*3.14*G*D*f)^2)*30*No*(2^(k/BW)-1));
den=(row*(dr*c).^2);
SARmin=num./den;

figure,
plot(sort((SARmin/1e10),'ascend'),'r-');
xlabel('distance between transmitter and reciever');
ylabel('SAR Value');
legend('Safe level');
title('safe tranmission range');
figure(3)
k=[0:500:3500];
fmax=2*1e9;
fmax1=3*1e9;
fmax2=6*1e9;
fmin= 900*1e6;
for j=1:numel(k)
    n=(BW./k).*(log2((fmax/fmin).^2*(2.^(k/BW)-1))+1);
    n1=(BW./k).*(log2((fmax1/fmin).^2*(2.^(k/BW)-1))+1);
    n2=(BW./k).*(log2((fmax2/fmin).^2*(2.^(k/BW)-1))+1);
end
plot(k,sort(abs(n./1e5),'descend'),'r-');
hold on
plot(k,sort(abs(n1./1e5),'descend'),'b-');
hold on
plot(k,sort(abs(n2./1e5),'descend'),'g-');
legend('fmax=6GHz','fmax=3GHz','fmax=2GHz')
xlabel('Transmission rate in Kpbs');
ylabel('No.of balanced nodes');
title('No. of balanced nodes');

figure(4)
k=[0:100];
fmax=2*1e9;
fmax1=3*1e9;
fmax2=6*1e9;
fmin= 900*1e6;
for j=1:numel(k)
    n=(BW./k).*(log2((fmax/fmin).^2*(2.^(k/BW)-1))+1);
    n1=(BW./k).*(log2((fmax1/fmin).^2*(2.^(k/BW)-1))+1);
    n2=(BW./k).*(log2((fmax2/fmin).^2*(2.^(k/BW)-1))+1);
end
color= ['r','g','b'];
h=bar(-(n./1e8));
set(h,'FaceColor',color(1));
h1=bar(-(n1./1e8));
set(h1,'FaceColor',color(2));
h2=bar(-(n2./1e8));
set(h2,'FaceColor',color(3));
xlabel('NodeID');
ylabel('SAR');
title('SAR for LTE 8 band');




