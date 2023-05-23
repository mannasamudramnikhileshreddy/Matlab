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

deadnodes = 0;
alivenodes = n;

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
    ri(i) = rand(1,1);
    %status(i) = 1;
    %alivenode(i) = i;
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

avgE = mean(Er);

i2 = 1;
i3 = 1;
for i1 = 1:n
    d(i1) = sqrt((xlocs(n+1)-xlocs(i1))^2 + (ylocs(n+1)-ylocs(i1))^2);  %distances
    if d(i1) < do
        Etx1(i2) = (k*Eelec) + ((k*Efs)*(d(i1))^2);                      %transmission energy
        Etx(i1) = (k*Eelec) + ((k*Efs)*(d(i1))^2);
        Et1I(i2) = i1;
        S(Et1I(i2)).type = 'Ych';            %selection of strategy for CH
        i2 = i2+1;
    else
        Etx2(i3) = (k*Eelec) + ((k*Emp)*(d(i1))^4);                      %transmission energy
        Etx(i1) = (k*Eelec) + ((k*Emp)*(d(i1))^4);
        Et2I(i3) = i1;
        S(Et2I(i3)).type = 'Nch';            %selection of strategy for CM
        i3 = i3+1;
    end
    Erx(i1) = k*Eelec;                                                   %receiving energy
end

%cost of becoming a CH
for i4 = 1:length(Et1I)
    Cch(i4) = Erx(Et1I(i4)) + Eaggr + Etx(Et1I(i4));
    Econ(Et1I(i4)) = Cch(i4);
end

%cost of becoming a CM
for i5 = 1:length(Et2I)
    Ccm(i5) = Etx(Et2I(i5));
    Econ(Et2I(i5)) = Ccm(i5);
end

%%calculation of Revenue Function
%For becoming a CH i.e.,'Ych'
for i6 = 1:length(Cch)
    H(i6) = Er(Et1I(i6))-Cch(i6);
end

%For becoming a CM i.e.,'Nch'
for i7 = 1:length(Ccm)
    C(i7) = Er(Et2I(i7))-Ccm(i7);
end


