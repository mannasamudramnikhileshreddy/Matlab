%%Simulation Parameters
xm = 150;
ym = 100;
xsink = 150;
ysink = 50;
n = 100;
k = 10;
rmax = 1000;
Packagesize  = 2000;
Packets = 400000;
Ctrlpackagesize  = 200;
Eelec  = 50*(10^-9);
Emp  = 0.0013*(10^-12);
Efs  = 10*(10^-12);
Eaggr  = 5*(10^-9);
Eo  = 1;
do = sqrt(Efs/Emp);
%%%END%%%

%%Network Generation
%
rng('default')
for i = 1:n
    xlocs(i) = randi([1,100],1,1);
    ylocs(i) = randi([1,100],1,1);
    S(i).xlocs = xlocs(i);
    S(i).ylocs = ylocs(i);
    S(i).Eo = Eo;
    S(i).type = 'N';
    Er(i) = Eo;
    Econ(i) = 0;
    status(i) = 1;
    alivenode(i) = i;
    plot(xlocs(i),ylocs(i),'bo','linewidth',0.5);
    hold on;
    grid on;
end

xlocs(n+1) = 150;
ylocs(n+1) = 50;
S(n+1).xlocs = xlocs(n+1);
S(n+1).ylocs = ylocs(n+1);
plot(xlocs(n+1),ylocs(n+1),'k^','linewidth',2);
grid on;

%%Energy Consumption Model
%

deadnodes = 0;
alivenodes = n;


i2 = 1;
i3 = 1;
for i1 = 1:alivenodes
    d(i1) = sqrt((xlocs(n+1)-xlocs(i1))^2 + (ylocs(n+1)-ylocs(i1))^2);  %distances
    if d(i1) < do
        Etx1(i2) = (k*Eelec) + ((k*Efs)*(d(i1))^2);                      %transmission energy
        Etx(i1) = (k*Eelec) + ((k*Efs)*(d(i1))^2);
        Et1I(i2) = i1;
        i2 = i2+1;
    else
        Etx2(i3) = (k*Eelec) + ((k*Emp)*(d(i1))^4);                      %transmission energy
        Etx(i1) = (k*Eelec) + ((k*Emp)*(d(i1))^4);
        Et2I(i3) = i1;
        i3 = i3+1;
    end
    Erx(i1) = k*Eelec;                                                   %receiving energy
end

%%Revenue Function
%
for r = 1:rmax
    r
%cost of common node
Ccm = Etx1;
%i42 = 1;
for i41 = 1:length(Etx1)
    for i61 = 1:length(Etx)
        if Etx1(i41) == Etx(i61)
            I1(i41) = i61;
            I(i61) = i61;
            Econ(i61) = Etx(i61);
        end
    end
end
for i42 = 1:length(Etx2)
    for i62 = 1:length(Etx)
        if Etx2(i42) == Etx(i62)
           I2(i42) = i62;
           I(i62) = i62;
           Econ(i62) = Etx(i62);
        end
    end
end

%cost of cluster head
for i4 = 1:length(Etx2)
    for i6 = 1:length(Etx)
        if Etx2(i4) == Etx(i6)
            Cch(i4) = Erx(i6) + Eaggr + Etx(i6);
            Econ(i6) = Erx(i6) + Eaggr + Etx(i6);
        end
    end
end

%residual energy of nodes

for re1 = 1:length(I1)
    Er1(re1) = Er(I1(re1)) - Econ(I1(re1));
    Er(I1(re1)) = Er(I1(re1)) - Econ(I1(re1));
end
for re2 = 1:length(I2)
    Er2(re2) = Er(I2(re2)) - Econ(I2(re2));
    Er(I2(re2)) = Er(I2(re2)) - Econ(I2(re2));
end

%calculation of revenue function for common node

for i7 = 1:length(Ccm)
    C(i7) = Er1(i7)-Ccm(i7);
end

%calculation of revenue function for common node

for i8 = 1:length(Cch)
    H(i8) = Er2(i8)-Cch(i8);
end

%sending the packets to the Base Station

for i9 = 1:length(C)
    %Packets = Packets-1;
    Er(i9) = Er(I1(i9))-(Ctrlpackagesize*Econ(I1(i9)));
end

for i10 = 1:length(H)
    Packets = Packets-1;
    Er(i10) = Er(I2(i10))-((Ctrlpackagesize*Packagesize*Econ(I2(i10)))/length(H));
end

for i11 = 1:(length(C)+length(H))
    S(i11).Er = Er(i11);
end

%if r == rmax
i13 = 1;
%n1 = 
for i12 = 1:length(alivenode)
    if Er(i12) <= 0
        deadnodes = deadnodes+1;
        dead(i13) = i12;
        status(dead(i13)) = 0;
        i13 = i13+1;
    end
    if status(i12) == 0
        alivenode(i12) = [];
        alivenodes = alivenodes-1;
    end
end
%end
end