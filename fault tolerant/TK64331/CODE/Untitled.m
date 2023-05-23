%%creation of a wsn
xm = 1000;
ym = 500;
n = 100;
initial = 0;
yfinal = 500;
rmax = 500;
Eele = 50*(10^-9);
Emp = 0.0013*(10^-12);
Efs = 10*(10^-12);
CHPS = 4000*(10^3);
CMPS = 200*(10^3);
R = 200;
rng('default');
xlocs = randi([1 1000],1,n);
ylocs = randi([1 500],1,n);
rnum = rand(1,n);
xsink = 500;
ysink = 250;
packets = 6000;
Packets(1:18) = 6000;
PKT2BS = 0;
TPKT2BS = 0;
TPKTS = 0;
TPKTS1 = 0;
TPKTS2 = 0;

figure
for i = 1:n
    plot(xlocs(i),ylocs(i),'ro','linewidth',0.2);
    S(i).xlocations = xlocs(i);
    S(i).ylocations = ylocs(i);
    S(i).Eo = 0.3;
    S(i).Er = S(i).Eo;
    S(i).type = 'N';
    S(i).G = 0;
    D(i) = sqrt((xsink - xlocs(i))^2 + (ysink - ylocs(i))^2);
    hold on;
    grid on;
end

plot(xsink,ysink,'k^','linewidth',1.7);
grid on;

S(n+1).xlocations = xsink;
S(n+1).ylocations = ysink;

%making grids
xlocs = cat(2,xlocs,xsink);
ylocs = cat(2,ylocs,ysink);

r = R/sqrt(5);
gridtype = grids(r,n,xlocs,ylocs,initial,yfinal);
Eo = 0.3;
for i2 = 1:n
    S(i2).gridtype = gridtype(i2);
    S(i2).rnum = rnum(i2);
    Er(i2) = Eo;
end

ugt = unique(gridtype);

%gt = 1;
%selection of cluster heads
for rmax = 1:2
 %   rmax    
for gt = 1:length(ugt)
    CH(gt) = clusterhead(gridtype,n,rnum,Eo,Er,gt);
end

for ii = 1:length(CH)
    S(CH(ii)).type = 'CH';
end

%gaussian integer
ii3 = 1;
for ii1 = 1:6
    for ii2  = 1:3
        GI(ii3) = ii1+ii2*j;
        d(ii3) = gcd(ii1,ii2);
        f(ii3) = ((ii1^2)+(ii2^2))/d(ii3);
        ii3 = ii3+1;
    end
end

%Routing Paths
sink = 101;
CHBS = cat(2,CH,sink);
[Xmin,Ymin,XminI,YminI,XminID,YminID] = Paths(xlocs,ylocs,CHBS);

%%Fault Detection
[Packets,PKT2BS,PKTS,PKTS2,PKTS1,NPackets] = FaultDetection(Packets,CH,gridtype,n);
end
