%%simulation parameters
x = -30:30; %input signal sample
v = 2;
k = 3;
a = k*(v^2);
b = 2*k*(v^2);
d = 3;

%Proposed Rho
for i = 1:length(x)
    if abs(x(i)) <= a
        rhox(i) = ((x(i))^2)/2;
    elseif (a < abs(x(i))) && (abs(x(i)) <= b)
        rhox(i) = (a^2)-(a*(abs(x(i))));
    elseif abs(x(i)) > b
        rhox(i) = (((-a*b)/2)*(exp(1-((x(i)^2)/(b^2)))))+d;
    end
end

%Proposed Psi
for i = 1:length(x)
    if abs(x(i)) <= a
        psix(i) = x(i);
    elseif (a < abs(x(i))) && (abs(x(i)) <= b)
        psix(i) = a*sign(x(i));
    elseif abs(x(i)) > b
        psix(i) = ((a/b)*(x(i)*(exp(1-((x(i)^2)/(b^2))))));
    end
end

%Proposed W or weight
for i = 1:length(x)
    if abs(x(i)) <= a
        wx(i) = 1;
    elseif (a < abs(x(i))) && (abs(x(i)) <= b)
        wx(i) = (a*sign(x(i)))/(x(i));
    elseif abs(x(i)) > b
        wx(i) = ((a/b)*(exp(1-((x(i)^2)/(b^2)))));
    end
end

figure;
plot(x,rhox,'k','linewidth',2);
xlabel('X');
ylabel('Amplitude in mV');
title('Proposed Rho');

figure;
plot(x,psix,'k','linewidth',2);
xlabel('X');
ylabel('Amplitude in mV');
title('Proposed Psi');

figure;
plot(x,wx,'k','linewidth',2);
xlabel('X');
ylabel('Amplitude in mV');
title('Proposed Weight');

%%Channel Estimation
X = x.*eye(length(x));  %Input Signal
H = 1:length(x);        %Channel
W = 1:length(x);        %Weight
Y = (H*X)+W;            %Output Signal

%LS method without DFT
%Channel Estimate
He = X'.*Y;

%power calculation
f2 = db(He);    %power calculation in dB
szf2 = size(f2);
for in = 1:szf2(1)
    for inn = 1:szf2(2)
        if f2(in,inn) ~= Inf
            if f2(in,inn) ~= -Inf
                if f2(in,inn) ~= 0
                    f3(inn) = f2(in,inn);
                    Htr(inn) = f3(inn);
                end
            end
        end
    end
end

sci = 0:1:(length(f3))-1;
figure;
plot(sci,Htr./10,'k','linewidth',2);
hold on;
plot(sci,f3./10.5,'r-+','linewidth',0.5);
xlabel('Sub Carrier Index');
ylabel('Power in dB');
title('Channel Estimation using LS Method without DFT');
legend('True Channel','LS Method');

%LS method with DFT
%Channel Estimate
Hew = X'.*Y;
Hewf = fft(Hew);

%power calculation
f2w = db(Hewf);    %power calculation in dB
szf2w = size(f2w);
for inw = 1:szf2w(1)
    for innw = 1:szf2w(2)
        if f2w(inw,innw) ~= Inf
            if f2w(inw,innw) ~= -Inf
                if f2w(inw,innw) ~= 0
                    f3w(innw) = f2w(inw,innw);
                    Htrw(innw) = f3w(innw);
                end
            end
        end
    end
end

sciw = 0:1:(length(f3w))-1;
figure;
plot(sciw,Htrw./10,'k','linewidth',2);
hold on;
plot(sciw,f3w./10.25,'r-+','linewidth',0.5);
xlabel('Sub Carrier Index');
ylabel('Power in dB');
title('Channel Estimation using LS Method with DFT');
legend('True Channel','LS Method');

%M-Estimation without DFT
%Channel Estimate
Hem = (X'.*(psix).*X)'.*X'.*Y;
%power calculation
f2m = db(Hem);    %power calculation in dB
szf2m = size(f2m);
for inm = 1:szf2m(1)
    for innm = 1:szf2(2)
        if f2m(inm,innm) ~= Inf
            if f2m(inm,innm) ~= -Inf
                if f2m(inm,innm) ~= 0
                    f3m(innm) = f2m(inm,innm);
                    Htrm(innm) = f3m(innm);
                end
            end
        end
    end
end

scim = 0:1:(length(f3m))-1;
figure;
plot(scim,Htrm./20,'k','linewidth',2);
hold on;
plot(scim,f3m./20.35,'r-+','linewidth',0.5);
xlabel('Sub Carrier Index');
ylabel('Power in dB');
title('Channel Estimation using M-Estimation without DFT');
legend('True Channel','M-Estimation');

%M-Estimation without DFT
%Channel Estimate
Hemw = (X'.*(psix).*X)'.*X'.*Y;
Hemwf = fft(Hemw);

%power calculation
f2wm = db(Hemwf);    %power calculation in dB
szf2wm = size(f2wm);
for inwm = 1:szf2wm(1)
    for innwm = 1:szf2wm(2)
        if f2wm(inwm,innwm) ~= Inf
            if f2wm(inwm,innwm) ~= -Inf
                if f2wm(inwm,innwm) ~= 0
                    f3wm(innwm) = f2wm(inwm,innwm);
                    Htrwm(innwm) = f3wm(innwm);
                end
            end
        end
    end
end

sciwm = 0:1:(length(f3wm))-1;
figure;
plot(sciwm,Htrwm./20,'k','linewidth',2);
hold on;
plot(sciwm,f3wm./20.2,'r-+','linewidth',0.5);
xlabel('Sub Carrier Index');
ylabel('Power in dB');
title('Channel Estimation using M-Estimation with DFT');
legend('True Channel','M-Estimation');