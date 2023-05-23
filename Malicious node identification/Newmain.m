
%simulation area
xm = 100;
ym = 100;

%simulation parameters
n = 100;    %Total number of nodes in the network
k = 4;      %No of CHs in the  network
r = 20;     %communication radius
Eo = 2;     %initial energy
Edis = 50*10^-9;%Energy Dissipation per sent/receiving a bit
Efs = 10*10^-12;%Energy consumed through freespace model
Emp = 100*10^-12;%Energy consumed through multipath model
pktsz = 80; %packet size of in bits
npktnod = 10; %number of packets from each node per round
tot_mal_nodes = 10; %Total number of Malicious Nodes
rmax = 200; %Maximum number of Rounds
theta = 150; %Fixed Parameter Maintenance
cntchs = zeros(1,n); %CHs counter
Beta = zeros(1,n); %CHs with abnormal communications counter
mue = zeros(1,n); %reputation maintenance function
delalp = 2; %number of normal communications in a period of time t
delbet = 2; %number of abnormal communications in a period of time t
dspdn = zeros(1,n);
sig = 1; %adjustment parameter
psi = 0.8; %Threshold value for medium and high reputation level
del = 0.5; %Threshold value for low and medium reputation level
eta = [0.1,0.2,0.3,0.4,0.5]; %adjust the size of tolerable trust span threshold
lts = [50,100,150,200,300]; %length of time series

%Generating a WSN Network
% rng('default');
figure;
for i = 1:n
    xloc(i) = rand(1,1)*xm;
    yloc(i) = rand(1,1)*ym;
    S(i).xloc = xloc(i);
    S(i).yloc = yloc(i);
    S(i).Eo = Eo;
    S(i).Er = Eo;
    Er(i) = Eo;
    plot(xloc(i),yloc(i),'ro');
    hold on;
    grid on;
end

xsink = xm/2;
ysink = ym/2;
S(n+1).xloc = xsink;
S(n+1).yloc = ysink;
plot(xsink,ysink,'k^','linewidth',2);
grid on;
title('Wireless Sensor Network');

%Formation of Clusters using K-means Clustering
locs = cat(1,xloc,yloc);
locsr = reshape(locs,100,2);
[idx,C] = kmeans(locsr,k);
uidx = unique(idx);

figure;
plot(locsr(:,1),locsr(:,2),'.','MarkerSize',12);
title 'Nodes of WSN Network';

% figure;
% plot(locsr(idx==1,1),locsr(idx==1,2),'r.','MarkerSize',12)
% hold on;
% plot(locsr(idx==2,1),locsr(idx==2,2),'b.','MarkerSize',12)
% hold on;
% plot(locsr(idx==3,1),locsr(idx==3,2),'g.','MarkerSize',12)
% hold on;
% plot(locsr(idx==4,1),locsr(idx==4,2),'k.','MarkerSize',12)
% hold on;
% plot(C(1,:),C(2,:),C(3,:),C(4,:),'kx',...
%      'MarkerSize',15,'LineWidth',3) 
%  hold on;
%  plot(xsink,ysink,'k^','MarkerSize',12,'Linewidth',3);
% legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Centroids','BS',...
%        'Location','SW')
% title 'Cluster Assignments and Centroids'
% hold off

%Finding Distances of Nodes from BS
d3 = 1;d4 = 1;d5 = 1;d6 = 1;
for d1 = 1:length(uidx)
    for d2 = 1:length(idx)
        if idx(d2) == uidx(d1)
            D(d2) = sqrt(((xsink - xloc(d2))^2) + ((ysink - yloc(d2))^2));
            if idx(d2) == 1
                D1(d3) = d2;
                D1i(d3) = D(d2);
                Er1(d3) = Er(d2);
                Er1i(d3) = d2;
                d3 = d3+1;
            elseif idx(d2) == 2
                D2(d4) = d2;
                D2i(d4) = D(d2);
                Er2(d4) = Er(d2);
                Er2i(d4) = d2;
                d4 = d4+1;
            elseif idx(d2) == 3
                D3(d5) = d2;
                D3i(d5) = D(d2);
                Er3(d5) = Er(d2);
                Er3i(d5) = d2;
                d5 = d5+1;
            elseif idx(d2) == 4
                D4(d6) = d2;
                D4i(d6) = D(d2);
                Er4(d6) = Er(d2);
                Er4i(d6) = d2;
                d6 = d6+1;
            end
        end
    end
end
LD1 = length(D1);LD2 = length(D2);
LD3 = length(D3);LD4 = length(D4);
sD1i = sort(D1i);sD2i = sort(D2i);
sD3i = sort(D3i);sD4i = sort(D4i);
% sEr1 = sort(Er1);sEr2 = sort(Er2);
% sEr3 = sort(Er3);sEr4 = sort(Er4);
LEr1 = length(Er1);LEr2 = length(Er2);
LEr3 = length(Er3);LEr4 = length(Er4);
mD1i = min(D1i);mD2i = min(D2i);
mD3i = min(D3i);mD4i = min(D4i);
mD1ii = find(D == mD1i);mD2ii = find(D == mD2i);
mD3ii = find(D == mD3i);mD4ii = find(D == mD4i);
mDii = cat(1,mD1ii,mD2ii,mD3ii,mD4ii);
status = ones(1,n);

%%Rounds
ii5 = 0;CHs = zeros(1,4);
for r2 = 1:1:rmax
    r2
    
    %Election of CHs for First Round
    if r2 == 1
        for cc = 1:length(CHs)
            CHs(cc) = mDii(cc);
        end
    else
    for a1 = 1:length(status)
        if status(a1) ~= 0
            for aa = 1:LEr1
                if a1 == Er1i(aa)
                    if Er(Er1i(aa)) == max(Er1)
                        CHs(1) = a1;
                    end
                end
            end
            for ab = 1:LEr2
                if a1 == Er2i(ab)
                    if Er(Er2i(ab)) == max(Er2)
                        CHs(2) = a1;
                    end
                end
            end
            for ac = 1:LD3
                if a1 == Er3i(ac)
                    if Er(Er3i(ac)) == max(Er3)
                        CHs(3) = a1;
                    end
                end
            end
            for ad = 1:LD4
                if a1 == Er4i(ad)
                    if Er(Er4i(ad)) == max(Er4)
                        CHs(4) = a1;
                    end
                end
            end
        end
    end
    end
    
    sEr1 = sort(Er1);sEr2 = sort(Er2);
    sEr3 = sort(Er3);sEr4 = sort(Er4);
    rint1 = randi([1,4],1,1);
    rint2 = randi([1,4],1,1);
    if r2 >= 50 && r2 <= 60
        CHsab(r2,:) = CHs;
    end
    
    %counter for CHs
    CHsR(r2,:) = CHs;
    if r2 > 50
        uCHsR = unique(CHsR);
        szCHsR = size(CHsR);
        for uch = 1:length(uCHsR)
            cnt = 0;
            for chs = 1:szCHsR(1)
                for chr = 1:szCHsR(2)
                    if uCHsR(uch) == CHsR(chs,chr)
                        cnt = cnt+1;
                        cntchs(uCHsR(uch)) = cnt;
                    end
                end
            end
        end
        
        %calculating number of normal and abnormal communications between
        %node and it's neighbors
        uCHsab = unique(CHsab);
        szCHsab = size(CHsab);
        uCHsab(1) = [];
        for y3 = 1:length(uCHsab)
            cntab = 0;
            for y1 = 1:szCHsab(1)
                for y2 = 1:szCHsab(2)
                    if uCHsab(y3) == CHsab(y1,y2)
                        cntab = cntab+1;
                        Beta(uCHsab(y3)) = cntab;   %abnormal communications
                    end
                end
            end
        end
        
        %calculating number of normal communications
        Alpha = cntchs-Beta;
        
        %reputation maintenance function
        mue = theta./(Alpha+Beta);
        
        %Node Reputation Model
        for rep = 1:n
            R(r2-50,rep) = ((mue*Alpha(rep))+delalp)/((mue*(Alpha(rep)+Beta(rep)))+delalp+delbet);
%             Rn(rep) = ((mue*Alpha(rep))+delalp)/((mue*(Alpha(rep)+Beta(rep)))+delalp+delbet);
        end
        
        if r2 >50 && r2 <= 60
        for repr = 1:n
            Rn(repr) = ((mue*Alpha(repr))+delalp)/((mue*(Alpha(repr)+Beta(repr)))+delalp+delbet);
        end
        
        for rep1 = 1:length(uCHsab)
            Rnn(rep1) = Rn(uCHsab(rep1));
        end
        end
        
        %Bench Mark Trust
        T = Rnn';
        
        %%Matrix of State
        %Energy Consumption
        if r2 > 100
            for ee1 = 1:length(uCHsab)
                Erra(ee1) = Err(50,uCHsab(ee1));
                Errc(ee1) = Err(60,uCHsab(ee1));
                Erra1(ee1) = Err(61,uCHsab(ee1));
                Errc1(ee1) = Err(71,uCHsab(ee1));
            end
            e = Errc-Erra;  %Energy Consumption of nodes before and after time t
            en = Errc1-Erra1;  %Energy Consumption of nodes after time t and for some period
            sznptx = size(numpktstxed);
            for erc = 1:length(uCHsab)
                for era = 1:szCHsab(1)
                    for erb = 1:szCHsab(2)
                        if era ~= 50
                        if uCHsab(erc) == CHsab(era,erb)
%                             Errb(erc) = Err(era,uCHsab(erc));
                            w(erc) = numpktstxed(era,erb)./numpktsrxed(era,erb);
%                             w2(erc) = numpktstxed2(era,erb)./numpktsrxed2(era,erb);
                        end
                        if era < 50
                        if uCHsab(erc) == CHsR(era,erb)
%                             Errb(erc) = Err(era,uCHsab(erc));
%                             w(erc) = numpktstxed(era,erb)./numpktsrxed(era,erb);
                            w2(erc) = numpktstxed2(era,erb)./numpktsrxed2(era,erb);
                            break;
                        end
                        end
                        end
                    end
                end
            end
            if length(w) < length(e)
                w(length(w)+1) = 0;
            end
%             for ww = 1:sznptx(1)
%                 for ww1 = 1:sznptx(2)
%                     w(ww,ww1) = numpktstxed(51:60,1:4)./numpktsrxed(51:60,1:4);  %Data Volume Calculation for a Node during time t
%                 end
%             end
%             w = reshape(w,1,length(e));
            
            %%sparsity of nodes calculation
            for snc = 1:n
                cntdsp = 0;snc2 = 1;
                for snc1 = 1:n
                    if snc ~= snc1
                        dsp(snc,snc1) = sqrt(((xloc(snc)-xloc(snc1))^2) + ((yloc(snc)-yloc(snc1))^2));
                        if dsp(snc,snc1) <= r
                            dspn(snc,snc2) = dsp(snc,snc1);
                            dspdn(snc) = dspdn(snc)+dspn(snc,snc2); %distances of node and it's neighbors
                            snc2 = snc2+1;
                            cntdsp = cntdsp+1;
                            cntden(snc) = cntdsp;   %density of nodes
                        end
                    end
                end
            end
            spar = dspdn./(cntden*r);   %sparsity calculation
            
            %%Matrix of State
            for ccn = 1:length(uCHsab)
                cnd(ccn) = cntden(uCHsab(ccn));
                sp(ccn) = spar(uCHsab(ccn));
            end
            
            St = cat(1,e,w,cnd,sp);
            St2 = cat(1,en,w2,cnd,sp);
            
            %%Matrix of Environmental Parameters
            q1 = St'*St;
            Q = (q1*St').*T;
            
            %%Benchmark of Trust
            Tn = Q.*St';
            Tn2 = St2'.*Q;
            szTn = size(Tn);
            
            %%Similarity
            %Gaussian Radial Basis Function
            GRB = (1/n)*(sum(exp(-((Tn2(:,1)-Tn(:,1))/(sig^2)))));
            
            %change direction of reputation
            if sum(Tn2(:,1)) < sum(Tn(:,1))
                cdr = 1;
            elseif sum(Tn2(:,1)) > sum(Tn(:,1))
                cdr = -1;
            end
            
            %similarity calculation
            sim = cdr*GRB;
            
            %Tolerable Trust Span Threshold
%             for tts = 1:length(Tn2)
%                 gam(tts) = eta(1)*(exp((sum(Tn2)/(n(1)*(psi-del)))-(del/(psi-del))));
%             end
         if r2 > 110   
            %Malicious Node Recognition
            nrn = 1;nmn = 1;mln = 1;
            for rp1 = 1:length(Rn)
                if Rn(rp1) > psi
                    nornod(nrn) = rp1;
                    nrn = nrn+1;
                else
                    if Rn(rp1) > del
                        nrmlnod(nmn) = rp1;
                        for gem = 1:length(eta)
                            for gem1 = 1:length(lts)
                                gam(gem,gem1) = ((eta(gem))*(exp((sum(Tn2(:,1)))/(lts(gem1)*(psi-del))-(del/(psi-del)))));
                            end
                        end
                        if abs(sim) < gam(5,5)
                            nornod(nrn) = rp1;
                            nrn = nrn+1;
                        else
                            if R(61,rp1) > Rn(rp1)
                                nornod(nrn) = rp1;
                                nrn = nrn+1;
                            else
                                malnod(mln) = rp1;
                                mln = mln+1;
                            end
                        end
                        nmn = nmn+1;
                    else
                        malnod(mln) = rp1;
                        mln = mln+1;
                    end
                end
            end
         end
        end
    end
    if r2 == rmax
        maln = abs(length(malnod)-length(nrmlnod))/4;
    end
    %Data Transmission Phase
    %Packets to CHs and then to BS
    %Energy Consumption
    for pp1 = 1:LEr1
        if r2 > 50 && r2 <= 57
        if Er1i(pp1) ~= CHs(1)
            Econ(r2,Er1i(pp1)) = Edis*pktsz*npktnod;
        end
        if Er1i(pp1) == CHs(1)
            Econ(r2,Er1i(pp1)) = (Edis*pktsz*npktnod*LEr1) + (Edis*pktsz*npktnod*(ceil(LEr1-rint1)));
            numpktsrxed(r2,1) = (LEr1-1)*npktnod;
            numpktstxed(r2,1) = ceil(LEr1-(LEr1/rint1))*npktnod;
        end
        elseif r2 > 57 && r2 <= 60
        if Er1i(pp1) ~= CHs(1)
            Econ(r2,Er1i(pp1)) = Edis*pktsz*npktnod;
        end
        if Er1i(pp1) == CHs(1)
            Econ(r2,Er1i(pp1)) = (Edis*pktsz*npktnod*LEr1) + (Edis*pktsz*npktnod*(ceil(LEr1-rint2)));
            numpktsrxed(r2,1) = (LEr1-1)*npktnod;
            numpktstxed(r2,1) = ceil(LEr1-(LEr1/rint2))*npktnod;
        end
        else
        if Er1i(pp1) ~= CHs(1)
            Econ(r2,Er1i(pp1)) = Edis*pktsz*npktnod;
%             numpktsrxed2(r2,1) = (LEr1-1)*npktnod;
%             numpktstxed2(r2,1) = ceil(LEr1-(LEr1/rint2))*npktnod;
        end
        if Er1i(pp1) == CHs(1)
            Econ(r2,Er1i(pp1)) = (Edis*pktsz*npktnod*LEr1) + (Edis*pktsz*npktnod*(LEr1-1));
            numpktsrxed2(r2,1) = (LEr1-1)*npktnod;
            numpktstxed2(r2,1) = ceil(LEr1-(LEr1/rint2))*npktnod;
        end
        end
    end
    
    for pp2 = 1:LEr2
        if r2 > 50 && r2 <= 57
        if Er2i(pp2) ~= CHs(2)
            Econ(r2,Er2i(pp2)) = Edis*pktsz*npktnod;
        end
        if Er2i(pp2) == CHs(2)
            Econ(r2,Er2i(pp2)) = (Edis*pktsz*npktnod*LEr2) + (Edis*pktsz*npktnod*(ceil(LEr2-rint1)));
            numpktsrxed(r2,2) = (LEr2-1)*npktnod;
            numpktstxed(r2,2) = ceil(LEr2-(LEr2/rint1))*npktnod;
        end
        elseif r2 > 57 && r2 <= 60
        if Er2i(pp2) ~= CHs(2)
            Econ(r2,Er2i(pp2)) = Edis*pktsz*npktnod;
        end
        if Er2i(pp2) == CHs(2)
            Econ(r2,Er2i(pp2)) = (Edis*pktsz*npktnod*LEr2) + (Edis*pktsz*npktnod*(ceil(LEr2-rint2)));
            numpktsrxed(r2,2) = (LEr2-1)*npktnod;
            numpktstxed(r2,2) = ceil(LEr2-(LEr2/rint2))*npktnod;
        end
        else
        if Er2i(pp2) ~= CHs(2)
            Econ(r2,Er2i(pp2)) = Edis*pktsz*npktnod;
        end
        if Er2i(pp2) == CHs(2)
            Econ(r2,Er2i(pp2)) = (Edis*pktsz*npktnod*LEr2) + (Edis*pktsz*npktnod*(LEr2-1));
            numpktsrxed2(r2,2) = (LEr1-1)*npktnod;
            numpktstxed2(r2,2) = ceil(LEr1-(LEr1/rint2))*npktnod;
        end
        end
    end
    
    for pp3 = 1:LEr3
        if r2 > 50 && r2 <= 57
        if Er3i(pp3) ~= CHs(3)
            Econ(r2,Er3i(pp3)) = Edis*pktsz*npktnod;
        end
        if Er3i(pp3) == CHs(3)
            Econ(r2,Er3i(pp3)) = (Edis*pktsz*npktnod*LEr3) + (Edis*pktsz*npktnod*(ceil(LEr3-rint1)));
            numpktsrxed(r2,3) = (LEr3-1)*npktnod;
            numpktstxed(r2,3) = ceil(LEr3-(LEr3/rint1))*npktnod;
        end
        elseif r2 > 57 && r2 <= 60
        if Er3i(pp3) ~= CHs(3)
            Econ(r2,Er3i(pp3)) = Edis*pktsz*npktnod;
        end
        if Er3i(pp3) == CHs(3)
            Econ(r2,Er3i(pp3)) = (Edis*pktsz*npktnod*LEr3) + (Edis*pktsz*npktnod*(ceil(LEr3-rint2)));
            numpktsrxed(r2,3) = (LEr3-1)*npktnod;
            numpktstxed(r2,3) = ceil(LEr3-(LEr3/rint2))*npktnod;
        end
        else
        if Er3i(pp3) ~= CHs(3)
            Econ(r2,Er3i(pp3)) = Edis*pktsz*npktnod;
        end
        if Er3i(pp3) == CHs(3)
            Econ(r2,Er3i(pp3)) = (Edis*pktsz*npktnod*LEr3) + (Edis*pktsz*npktnod*(LEr3-1));
            numpktsrxed2(r2,3) = (LEr1-1)*npktnod;
            numpktstxed2(r2,3) = ceil(LEr1-(LEr1/rint2))*npktnod;
        end
        end
    end
    
    for pp4 = 1:LEr4
        if r2 > 50 && r2 <= 57
        if Er4i(pp4) ~= CHs(4)
            Econ(r2,Er4i(pp4)) = Edis*pktsz*npktnod;
        end
        if Er4i(pp4) == CHs(4)
            Econ(r2,Er4i(pp4)) = (Edis*pktsz*npktnod*LEr4) + (Edis*pktsz*npktnod*(ceil(LEr4-rint1)));
            numpktsrxed(r2,4) = (LEr4-1)*npktnod;
            numpktstxed(r2,4) = ceil(LEr4-(LEr4/rint1))*npktnod;
        end
        elseif r2 > 57 && r2 <= 60
        if Er4i(pp4) ~= CHs(4)
            Econ(r2,Er4i(pp4)) = Edis*pktsz*npktnod;
        end
        if Er4i(pp4) == CHs(4)
            Econ(r2,Er4i(pp4)) = (Edis*pktsz*npktnod*LEr4) + (Edis*pktsz*npktnod*(ceil(LEr4-rint2)));
            numpktsrxed(r2,4) = (LEr4-1)*npktnod;
            numpktstxed(r2,4) = ceil(LEr4-(LEr4/rint2))*npktnod;
        end
        else
        if Er4i(pp4) ~= CHs(4)
            Econ(r2,Er4i(pp4)) = Edis*pktsz*npktnod;
        end
        if Er4i(pp4) == CHs(4)
            Econ(r2,Er4i(pp4)) = (Edis*pktsz*npktnod*LEr4) + (Edis*pktsz*npktnod*(LEr4-1));
            numpktsrxed2(r2,4) = (LEr1-1)*npktnod;
            numpktstxed2(r2,4) = ceil(LEr1-(LEr1/rint2))*npktnod;
        end
    end
    end
    
    %Residual Energy Calculation
    for en1 = 1:LEr1
        Er(Er1i(en1)) = Er(Er1i(en1))-Econ(r2,Er1i(en1));
        Er1(en1) = Er(Er1i(en1));
        Err(r2,Er1i(en1)) = Er(Er1i(en1));
    end
    
    for en2 = 1:LEr2
        Er(Er2i(en2)) = Er(Er2i(en2))-Econ(r2,Er2i(en2));
        Er2(en2) = Er(Er2i(en2));
        Err(r2,Er2i(en2)) = Er(Er2i(en2));
    end
    
    for en3 = 1:LEr3
        Er(Er3i(en3)) = Er(Er3i(en3))-Econ(r2,Er3i(en3));
        Er3(en3) = Er(Er3i(en3));
        Err(r2,Er3i(en3)) = Er(Er3i(en3));
    end
    
    for en4 = 1:LEr4
        Er(Er4i(en4)) = Er(Er4i(en4))-Econ(r2,Er4i(en4));
        Er4(en4) = Er(Er4i(en4));
        Err(r2,Er4i(en4)) = Er(Er4i(en4));
    end
end

%Network Performance Evaluation Index
load('trstn.mat')
trstnd = n-length(malnod);
delt = 0:0.1:1;
figure;
plot(delt,sort(trstn),'k','linewidth',2);
xlabel('Delta');
ylabel('Number of Trusted Nodes');
title('Relationship between Delta and Trusted Nodes');

load('malnodn.mat')
malnodni = length(malnod);
psit = 0:0.1:1;
figure;
plot(psit,sort(malnodn),'k','linewidth',2);
xlabel('Psi');
ylabel('Number of Malicious Nodes');
title('Relationship between Psi and Malicious Nodes');

%Recognition and False Positive Percentage Calculations
load('RPpin.mat')
load('FPPi.mat')
RPpind = (maln/tot_mal_nodes)*100;
FPP = (length(nrmlnod)/(n-tot_mal_nodes))*100;

% load('RPi.mat')

% RPi(5,5) = RP;
% For a better RP range take nodes in Low Reputation Range as Malicious
% Nodes
% FPPi(5,5) = FPP;

figure;
plot(lts,sort(RPpin,'ascend'),'linewidth',2);
xlabel('Length of Time Series (n)');
ylabel('Recognition Percentage');
title('Relationship between eta and RP');
legend('η1=0.1','η2=0.2','η3=0.3','η4=0.4','η5=0.5','location','southeast');

figure;
plot(lts,sort(FPPi,'descend'),'linewidth',2);
xlabel('Length of Time Series (n)');
ylabel('False Positive Percentage');
title('Relationship between eta and FPP');
legend('η1=0.1','η2=0.2','η3=0.3','η4=0.4','η5=0.5','location','southeast');

RPpinn = RPpin./100;
figure;
plot3(lts,eta,RPpinn,'k->','linewidth',2);
xlabel('Length of Time Series (n)');
ylabel('η');
zlabel('Recognition Percentage (RP)');
title('Relationship between different η and RP');
legend('η1=0.1','η2=0.2','η3=0.3','η4=0.4','η5=0.5','location','southeast');

FPPin = FPPi./100;
figure;
plot3(lts,eta,FPPin,'k->','linewidth',2);
xlabel('Length of Time Series (n)');
ylabel('η');
zlabel('False Positive Percentage (FPP)');
title('Relationship between different η and FPP');
legend('η1=0.1','η2=0.2','η3=0.3','η4=0.4','η5=0.5','location','Northeast');

sRPpin = sort(RPpin);
sRPpinn = sRPpin(5,:);
nmalnod = 5:6:30;
figure;
plot(nmalnod,sort(sRPpinn,'ascend'),'k->','linewidth',2);
xlabel('Number of Mal Nodes');
ylabel('Recognition Percentage');
title('Number of Mal Nodes vs Recognition Percentage');

sFPPi = sort(FPPi);
sFPPin = sFPPi(5,:);
nnornod = 50:60:300;
figure;
plot(nnornod,sort(sFPPin,'descend'),'k->','linewidth',2);
xlabel('Number of Normal Nodes');
ylabel('False Positive Percentage');
title('Number of Normal Nodes vs False Positive Percentage');

figure;
for inde = 1:n
    for inde1 = 1:4
        if inde == CHs(inde1)
            plot(xloc(inde),yloc(inde),'go','linewidth',2);
            hold on;
        end
    end
end

for indr = 1:n
    for indr1 = 1:length(malnod)
        if indr == malnod(indr1)
            plot(xloc(indr),yloc(indr),'ro','linewidth',2);
            hold on;
        end
    end
end

for indb = 1:length(nornod)
    cindb = 0;cindb1 = 0;
    for indb1 = 1:4
        if indb ~= CHs(indb1) % || (indb == CHs(indb2))
            cindb = cindb+1;
        end
    end
    
    for indb2 = 1:length(malnod)
        if indb ~= malnod(indb2) % || (indb == CHs(indb2))
            cindb1 = cindb1+1;
        end
    end
    
    if cindb == 4 && cindb1 == (length(malnod))
        plot(xloc(indb),yloc(indb),'ko','linewidth',2);
        hold on;
    end
end

legend('green CHs','red Malnodes','black normnodes');
title('Wireless Sensor Network');

