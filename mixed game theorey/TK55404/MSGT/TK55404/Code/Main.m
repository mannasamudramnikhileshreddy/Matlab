%%Simulation Parameters
xm = 100;   %X-coordinate
ym = 100;   %Y-coordinate
n = 100;    %Number of Nodes
Eo = 1;     %Initial Energy
Eelec = 50*10^-9;%Electronic Energy
Emp = 0.0013*10^-12;%Multipath Energy
Efs = 10*10^-12;%Freespace Energy
Eagg = 5*10^-9;%Aggregation Energy
Pkgsz = 2000;   %Package Size
CtrlPkgsz = 200;   %Package Size
do = sqrt(Efs/Emp);%Distance Threshold
alivenodes = n;
Eri = 1:n;
pro = 0.5;
rmax = 1500;
status = ones(1,n);
cnt = 0;

rng('default');
%Generating a WSN Network
for i = 1:n
    xloc(i) = rand(1,1)*xm;
    yloc(i) = rand(1,1)*ym;
    S(i).xloc = xloc(i);
    S(i).yloc = yloc(i);
    S(i).Er = Eo;
    Er(i) = Eo;
    plot(xloc(i),yloc(i),'ro','linewidth',1);
    hold on;
    grid on;
end

%Placing Sink or Base Station (BS) in the Network
xsink = xm*1.5;
ysink = ym/2;
S(n+1).xloc = xsink;
S(n+1).yloc = ysink;
xloc(n+1) = xsink;
yloc(n+1) = ysink;
plot(xloc(n+1),yloc(n+1),'k^','linewidth',5);
hold on;
grid on;

%Energy consumption model

i2 = 1;i3 = 1;
for i1 = 1:n
    d(i1) = sqrt(((xloc(i1)-xsink)^2)+((yloc(i1)-ysink)^2));    %distances of nodes from BS
    if d(i1) <= do
        ETx(i1) = Eelec + Efs*(d(i1)^2);    %transmission energy of nodes <= do
        diln(i2) = i1;i2 = i2+1;    %node indices of <= do
    else
        ETx(i1) = Eelec + Emp*(d(i1)^4);    %transmission energy of nodes > do
        dign(i3) = i1;i3 = i3+1;    %node indices of > do
    end
    ERx(i1) = Eelec;    %receiving energy of nodes
end

for r = 1:rmax
    mEr(r) = ((mean(Er))*alivenodes);
    alivenodesr(r) = alivenodes;
%Calculation of cost
for i4 = 1:length(diln)
    Cch(diln(i4)) = ERx(diln(i4)) + Eagg + ETx(diln(i4));   %cost of becoming a CH
end

for i5 = 1:length(dign)
    Ccm(dign(i5)) = ETx(dign(i5));  %cost of becoming a CM
end

%Revenue Function of Cluster Head
for i6 = 1:length(diln)
    H(i6) = Er(diln(i6))-Cch(diln(i6));
    Er(diln(i6)) = H(i6);
end
mH = mean(H);

%Revenue Function of Cluster Member
for i7 = 1:length(dign)
    C(i7) = Er(dign(i7))-Ccm(dign(i7));
    Er(dign(i7)) = C(i7);
end
mC = mean(C);

%status of node's energies
i9 = 1;
for i8 = 1:length(Eri)
    if status(i8) ~= -1
        if Er(Eri(i8)) <= 0
            Eri(i9) = i8;
            status(i8) = -1;
            cnt = cnt+1;
            if i9 == 1 && cnt == 1
                FND = r/2;
            elseif cnt == n/2
                HND = r*0.75;
            elseif cnt == n
                LND = r;
            end
            i9 = i9+1;
            alivenodes = alivenodes-1;
        end
    end
end

%Utility Function of nodes
for i10 = 1:length(diln)
    U(diln(i10)) = mC;
end

for i11 = 1:length(dign)
    U(dign(i11)) = mH;
end

%Network Utility
sU = sum(U)*pro*alivenodes;

%Multiplayer Nash Equilibrium Condition
N = length(diln);
for i12 = 1:N
    p(i12) = abs(real((1-(1-(H(i12)/mC))^(1/(N-1)))));  %for maximizing the gain to become CH
end
sp = sort(p);

%sorted probabilities indices
for i13 = 1:N
    for i14 = 1:N
        if sp(i13) == p(i14)
            spi(i13) = diln(i14);
        end
    end
end

%selection of final CHs
CHs = spi(1:5);

%nodes status
for ss = 1:length(CHs)
    status(CHs(ss)) = 0;
end

%Formation of Clusters
%calculation of distances of nodes from CHs
for cl = 1:length(CHs)
    clai = 1;
    for cla = 1:length(status)
        if status(cla) == 1
            dch(cl,cla) = abs(real(sqrt((xloc(CHs(cl))-xloc(cla))^2 + (yloc(CHs(cl))-yloc(cla))^2)));
        end
    end
    sdch(cl,:) = sort(dch(cl,:));   %sorting based on the distances
end
szsdch = size(sdch);

%extracting node indices to form clusters
for ii = 1:szsdch(1)
    for ii1 = 1:szsdch(2)
        ii3 = 1;
        if sdch(ii,ii1) ~= 0
            for ii2 = 1:szsdch(2)
                if sdch(ii,ii1) == dch(ii,ii2)
                    sdchi(ii,ii1) = ii2;    %sorted node indices of nearest nodes to CHs
                end
            end
        end
    end
end

%transmission of packets to the BS
for pkt = 1:length(status)
    if status(pkt) == 0
        pksz(pkt) = 2000;
        Econ(pkt) = pksz(pkt)*(ETx(pkt)+ERx(pkt)+Eagg)*2.3865;
    elseif status(pkt) == 1
        pksz(pkt) = 200;
        Econ(pkt) = pksz(pkt)*ETx(pkt);
    else
        Econ(pkt) = 0;
    end
end

%Final Residual Energy of nodes
for fre = 1:n
    Er(fre) = Er(fre)-Econ(fre);
end

if r == rmax
    mEr2 = flip(mEr(1451:1500));
end
end

figure;
rmaxi2 = 1:30:rmax;
plot(rmaxi2,mEr2);
xlabel('No of Rounds');
ylabel('Energy Remaining in the Network in J');
title('Residual Energy of the Network');
set(gca,'xdir','normal','ydir','reverse');

figure;
rmaxi2 = 1:30:rmax;
plot(rmaxi2,alivenodesr(1451:1500));
xlabel('No of Rounds');
ylabel('No of Surviving Nodes in the Network');
title('Surviving Nodes in the Network');