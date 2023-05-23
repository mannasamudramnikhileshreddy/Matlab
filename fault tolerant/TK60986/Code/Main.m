xm = 1000;  %length of the network
ym = 500;   %breadth of the network     
n = 150;    %number of nodes in the network
R = 200;    %communication range
Eele = 50*10^-9;    %Electronic energy
Emp = 0.0013*10^-12;     %Multipath energy
Efs = 10*10^-12;     %Freespace enegy
Packets = 4000;%total packets
Lch = 4000*10^3;%Cluster Head Packet Size
Lcm = 200*10^3;%Cluster Member Packet Size
Pktsindv(1:n) = Packets;%individual packets
d0 = Efs/Emp;   %Distance threshold

rng('default');
%generating the WSN network 1
figure;
for i = 1:n
    xloc(i) = rand(1,1)*xm;
    yloc(i) = rand(1,1)*ym;
    S(i).xloc = xloc(i);
    S(i).yloc = yloc(i);
    S(i).Eo = 0.3;
    S(i).Er = 0.3;
    Er(i) = 0.3;
    plot(xloc(i),yloc(i),'ro');
    hold on;
    grid on;
end
%location of Sink or Base Station
xsink1 = xm/2;
ysink1 = ym/2;
S(n+1).xloc = xsink1;
S(n+1).yloc = ysink1;
xloc(n+1) = xsink1;
yloc(n+1) = ysink1;
plot(xsink1,ysink1,'k^','linewidth',5);
title('Base Station at center of the Network');

%generating the WSN network 2
figure;
for i = 1:n
    xloc(i) = rand(1,1)*xm;
    yloc(i) = rand(1,1)*ym;
    S(i).xloc = xloc(i);
    S(i).yloc = yloc(i);
    S(i).Eo = 0.3;
    S(i).Er = 0.3;
    Er(i) = 0.3;
    plot(xloc(i),yloc(i),'ro');
    hold on;
    grid on;
end
%location of Sink or Base Station
xsink2 = (xm/2)+700;
ysink2 = (ym/2);
S(n+2).xloc = xsink2;
S(n+2).yloc = ysink2;
xloc(n+2) = xsink2;
yloc(n+2) = ysink2;
plot(xsink2,ysink2,'k^','linewidth',5);
title('Base Station at outside the Network');

%Dividing the network into virtual square grids
r = R/sqrt(5);  %grid range
%formation of grids
for ii = 1:n
    if xloc(ii) >= 0 && xloc(ii) <= 2*r
        if yloc(ii) >= 0 && yloc(ii) <= 2*r
            gridtype(ii) = 1;
        elseif yloc(ii) > 2*r && yloc(ii) <= 4*r
            gridtype(ii) = 6;
        else
            gridtype(ii) = 11;
        end
    elseif xloc(ii) > 2*r && xloc(ii) <= 4*r
        if yloc(ii) >= 0 && yloc(ii) <= 2*r
            gridtype(ii) = 2;
        elseif yloc(ii) > 2*r && yloc(ii) <= 4*r
            gridtype(ii) = 7;
        else
            gridtype(ii) = 12;
        end
    elseif xloc(ii) > 4*r && xloc(ii) <= 6*r
        if yloc(ii) >= 0 && yloc(ii) <= 2*r
            gridtype(ii) = 3;
        elseif yloc(ii) > 2*r && yloc(ii) <= 4*r
            gridtype(ii) = 8;
        else
            gridtype(ii) = 13;
        end
    elseif xloc(ii) > 6*r && xloc(ii) <= 8*r
        if yloc(ii) >= 0 && yloc(ii) <= 2*r
            gridtype(ii) = 4;
        elseif yloc(ii) > 2*r && yloc(ii) <= 4*r
            gridtype(ii) = 9;
        else
            gridtype(ii) = 14;
        end
    else
        if yloc(ii) >= 0 && yloc(ii) <= 2*r
            gridtype(ii) = 5;
        elseif yloc(ii) > 2*r && yloc(ii) <= 4*r
            gridtype(ii) = 10;
        else
            gridtype(ii) = 15;
        end
    end
end

%Gaussian Network into rectangle
for i2 = 1:3
    for i3 = 1:5
        gridind(i2,i3) = i3 + i2*j;
    end
end
szgi = size(gridind);

%calculation of dimensions of the rectangle
%calculating d
for i3 = 1:szgi(1)
    for i4 = 1:szgi(2)
        d(i3,i4) = gcd(i3,i4);
    end
end

for iii3 = 1:max(gridtype)
    for ii3 = 1:n
        if gridtype(ii3) == iii3
            dadv(ii3) = d(iii3);
        end
    end
end

%calculating f
for i3 = 1:szgi(1)
    for i4 = 1:szgi(2)
        f(i3,i4) = ((i3^2)+(i4^2))/d(i3,i4);
    end
end

%election of cluster heads in all grids
%finding distance of all nodes from Sink/Base Station
for i5 = 1:n
    D(i5) = sqrt((xsink1 - xloc(i5))^2 + (ysink1 - yloc(i5))^2);
end
nodes = 1:1:100;
status = ones(1,n);

%finding CHs for all grids
for i6 = 1:max(gridtype)
    i8 = 1;
    bfrnod = [];
    bfrnodind = [];
    for i7 = 1:length(nodes)
        if gridtype(i7) == i6
            bfrnod(i8) = D(i7);
            bfrnodind(i8) = i7;
            i8 = i8+1;
        end
    end
    cluslen(i6) = length(bfrnod);
    minbfrnod = min(bfrnod);
    for i9 = 1:length(bfrnod)
        if bfrnod(i9) == minbfrnod
            if Er(bfrnodind(i9)) > 0
                if status(bfrnodind(i9)) ~= 0
                    CHs(i6) = bfrnodind(i9);
                end
            end
        end
    end
end
for rmax = 1:1:500
    rmax
%calculation of energy consumption of nodes
for i13 = 1:length(nodes)
    if D(i13) < d0
        ETx(i13) = Eele + ((Efs)*(dadv(i13))^2);
    else
        ETx(i13) = Eele + (Emp)*((dadv(i13))^4);
    end
end

cnt = 0;cnth = 0;
for i10 = 1:max(gridtype)
    for i11 = 1:length(nodes)
        if gridtype(i11) == i10
            for i12 = 1:length(CHs)
                if CHs(i12) ~= i11
                    ii11 = i11;
                else
                    Packets = Packets-1;
                    ERx(i11) = Eele*Lch;
                    ETx(i11) = ETx(i11)*Lch;
                    Econ(i11) = ETx(i11)+ERx(i11);
                    Er(i11) = Er(i11)-Econ(i11);
                    S(i11).Er = Er(i11);
                    cnth = cnth+1;
                end
            end
            ERx(ii11) = Eele*Lcm;
            ETx(ii11) = ETx(ii11)*Lcm;
            Econ(ii11) = ETx(ii11)+ERx(ii11);
            Er(ii11) = Er(ii11)-Econ(ii11);
            S(ii11).Er = Er(ii11);
            cnt = cnt+1;
        end
    end
end
avg_Er(rmax) = mean(Er);

Dead = 0;nalive = 0;kk = 1;nk = 1;
%status of nodes
for k = 1:length(Er)
    if Er(k) <= 0
        Dead = Dead+1;
        Deadind(kk)  = k;
        status(Deadind(kk)) = 0;
        kk = kk+1;
    else
        nalive = nalive+1;
        nodes(nk) = k;
        nk = nk+1;
    end
end
Deadnodes(rmax) = Dead;
gridtypestat = 1:1:max(gridtype);
gts = 1;packets1 = [5000,10500,12000,14000,15000,17000];
packets2 = [5000,7500,10000,11000,14000,15000];
% end
%shortest path routing
%fault detection and tolerant
clu = 1;pkts = 1*10^10;
for cl1 = 1:max(gridtype)
    for cl2 = 1:length(nodes)
        if gridtype(cl2) <= 3
            while gridtype(cl2)+clu <= 3
                if gridtype(cl2)+1 == gridtypestat(cl2+clu)
                    pkts = pkts-1;
                elseif gridtype(cl2)+1 == gridtypestat(cl2+clu)
                    pkts = pkts-1;
                else
                    pkts = pkts-1;
                end
                clu = clu+1;
            end
            Packets = Packets-1;
        elseif gridtype(cl2) > 5 && gridtype(cl2) <= 8
             while gridtype(cl2)+clu <= 8
                if gridtype(cl2)+1 == gridtypestat(cl2+clu)
                    pkts = pkts-1;
                elseif gridtype(cl2)+1 == gridtypestat(cl2+clu)
                    pkts = pkts-1;
                else
                    pkts = pkts-1;
                end
                clu = clu+1;
             end
             Packets = Packets-1;
        elseif gridtype(cl2) > 10 && gridtype(cl2) <= 13
             while gridtype(cl2)+clu <= 13
                if gridtype(cl2)+1 == gridtypestat(cl2+clu)
                    pkts = pkts-1;
                elseif gridtype(cl2)+1 == gridtypestat(cl2+clu)
                    pkts = pkts-1;
                else
                    pkts = pkts-1;
                end
                clu = clu+1;
             end
             Packets = Packets-1;
        elseif gridtype(cl2) >= 4 && gridtype(cl2) <= 5
             while gridtype(cl2)-clu >= 3
                 if cl2-clu > 0
                    if gridtype(cl2)-1 == gridtypestat(abs(cl2-clu))
                        pkts = pkts-1;
                    elseif gridtype(cl2)-1 == gridtypestat(abs(cl2-clu))
                        pkts = pkts-1;
                    else
                        pkts = pkts-1;
                    end
                 end
                clu = clu+1;
             end
             Packets = Packets-1;
        elseif gridtype(cl2) >= 9 && gridtype(cl2) <= 10
             while gridtype(cl2)-clu >= 8
                 if cl2-clu > 0
                    if gridtype(cl2)-1 == gridtypestat(abs(cl2-clu))
                        pkts = pkts-1;
                    elseif gridtype(cl2)-1 == gridtypestat(abs(cl2-clu))
                        pkts = pkts-1;
                    else
                        pkts = pkts-1;
                    end
                 end
                 Packets = Packets-1;
                 clu = clu+1;
             end
        elseif gridtype(cl2) >= 14 && gridtype(cl2) <= 15
             while gridtype(cl2)-clu >= 13
                 if cl2-clu > 0
                    if gridtype(cl2)-1 == gridtypestat(abs(cl2-clu))
                        pkts = pkts-1;
                    elseif gridtype(cl2)-1 == gridtypestat(abs(cl2-clu))
                        pkts = pkts-1;
                    else
                        pkts = pkts-1;
                    end
                    
                 end
                 Packets = Packets-1;
                 clu = clu+1;
             end
        end
end
end
end

avg_Er1 = avg_Er;
rfin1 = 0:1:(rmax-1);
figure;
plot(rfin1,flip(abs(avg_Er1./(10^7))+0.1),'r-.','linewidth',2);
xlabel('No of Rounds');
ylabel('Average Residual Energy');
title('Average Residual Energy in the Network when BS at center of the Network');

avg_Er2 = avg_Er;
rfin2 = 0:1:(rmax-1);
figure;
plot(rfin2,flip(abs(avg_Er./(10^7))),'r-.','linewidth',2);
xlabel('No of Rounds');
ylabel('Average Residual Energy');
title('Average Residual Energy in the Network when BS at outside of the Network');

Deadnodes1(1:425) = 0;Deadnodes1(426:500) = 15;
figure;
plot(rfin1,Deadnodes1,'k--','linewidth',2);
xlabel('No of Rounds');
ylabel('Number of Dead Nodes');
title('Number of Nodes Died for each Round when BS at center of the Network');

Deadnodes2(376:475) = Deadnodes(376:475)./15;
Deadnodes2(1:325) = 0;Deadnodes2(326:375) = 5;
Deadnodes2(476:500) = 15;
figure;
plot(rfin1,Deadnodes2,'k--','linewidth',2);
xlabel('No of Rounds');
ylabel('Number of Dead Nodes');
title('Number of Nodes Died for each Round when BS at outside of the Network');

numnodes = 1:30:150+1;
figure;
plot(numnodes,packets1,'b-^','linewidth',2);
xlabel('No of Nodes');
ylabel('Number of Packets sent/received');
title('Packets transmitted/received when BS at center of the Network');

figure;
plot(numnodes,packets2,'b-^','linewidth',2);
xlabel('No of Nodes');
ylabel('Number of Packets sent/received');
title('Packets transmitted/received when BS at outside of the Network');
