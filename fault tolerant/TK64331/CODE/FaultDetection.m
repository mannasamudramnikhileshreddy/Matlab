function [Packets,Status,PKT2BS,PKTS,PKTS1,PKTS2,NPackets] = FaultDetection(Packets,CH,gridtype,n)
OPackets = Packets;

i1 = 1;
i2 = 1;
for i = 1:n
    if i ~= i1
        if gridtype(i) == i1
            CM(i2,:) = i;
            i2 = i2+1;
        end
        i1 = i1+1;
    end
end
z1 = 1;
for z = 1:length(status)
    if rmax >= 2
        if status(z) == 0
            sind(z1) = z;
            z1 = z1+1;
        end
    end
end
 
if rmax == 1
    for i = 7:12
        if i == 7
            Packets(1) = Packets(1)-1;
            Packets(13) = Packets(13)-1;
            Packets(7) = Packets(7)-1;
            PKTS1= (OPackets(7)-Packets(7))+((OPackets(1)-Packets(1))+(OPackets(13)-Packets(13)));
        elseif i == 8
            Packets(6) = Packets(6)-1;
            Packets(18) = Packets(18)-1;
            Packets(12) = Packets(12)-1;
            PKTS2 = (OPackets(12)-Packets(12))+(OPackets(6)-Packets(6))+(OPackets(18)-Packets(18));
        elseif i == 9
            Packets(5) = Packets(5)-1;
            Packets(11) = Packets(11)-1;
            Packets(17) = Packets(17)-1;
            PKTS2 = (OPackets(5)-Packets(5))+(OPackets(11)-Packets(11))+(OPackets(17)-Packets(17))+PKTS2;
        elseif i == 10
            Packets(4) = Packets(4)-1;
            Packets(10) = Packets(10)-1;
            Packets(16) = Packets(16)-1;
            PKTS2 = (OPackets(4)-Packets(4))+(OPackets(10)-Packets(10))+(OPackets(16)-Packets(16))+PKTS2;
        elseif i == 11
            Packets(3) = Packets(3)-1;
            Packets(9) = Packets(9)-1;
            Packets(15) = Packets(15)-1;
            PKTS2 = (OPackets(3)-Packets(3))+(OPackets(9)-Packets(9))+(OPackets(15)-Packets(15))+PKTS2;
        else
            Packets(2) = Packets(2)-1;
            Packets(14) = Packets(14)-1;
            Packets(8) = Packets(8)-1;
            PKTS = (OPackets(8)-Packets(8))+(OPackets(2)-Packets(2))+(OPackets(14)-Packets(14))+PKTS2+PKTS1;
        end
    end
else
    for i = 7:12
        if Status(7) ~= 0
            Packets(1) = Packets(1)-1;
            Packets(13) = Packets(13)-1;
            Packets(7) = Packets(7)-1;
            PKTS1= (OPackets(7)-Packets(7))+((OPackets(1)-Packets(1))+(OPackets(13)-Packets(13)));
        elseif Status(7) == 0
            if Status(2) ~= 0
                Packets(1) = Packets(1)-1;
                Packets(2) = Packets(2)-1;
                PKTS11 = (OPackets(1)-Packets(1))+(OPackets(2)-Packets(2));
            elseif Status(2) == 0
                Packets(1) = Packets(1)-1;
                PKTS11 = (OPackets(1)-Packets(1));
            end
            if Status(14) ~= 0
                Packets(13) = Packets(13)-1;
                Packets(14) = Packets(14)-1;
                PKTS12 = (OPackets(13)-Packets(13))+(OPackets(14)-Packets(14));
            elseif Status(14) == 0
                Packets(13) = Packets(13)-1;
                PKTS12 = (OPackets(13)-Packets(13));
            end
            PKTS1 = PKTS11+PKTS12;
        end
        if Status(12) ~= 0
            Packets(6) = Packets(6)-1;
            Packets(18) = Packets(18)-1;
            Packets(12) = Packets(12)-1;
            PKTS1= (OPackets(12)-Packets(12))+((OPackets(6)-Packets(6))+(OPackets(18)-Packets(18)));
        elseif Status(12) == 0
            if Status(5) ~= 0
                Packets(6) = Packets(6)-1;
                Packets(5) = Packets(5)-1;
                PKTS11 = (OPackets(6)-Packets(6))+(OPackets(5)-Packets(5));
            elseif Status(5) == 0
                Packets(6) = Packets(6)-1;
                PKTS11 = (OPackets(1)-Packets(1));
            end
            if Status(17) ~= 0
                Packets(18) = Packets(18)-1;
                Packets(17) = Packets(17)-1;
                PKTS12 = (OPackets(17)-Packets(17))+(OPackets(18)-Packets(18));
            elseif Status(17) == 0
                Packets(18) = Packets(18)-1;
                PKTS12 = (OPackets(18)-Packets(18));
            end
            PKTS1 = PKTS11+PKTS12;
        end
end
PKT2BS = PKTS;

for i1 = 1:length(CH)
    NPackets(i1) = OPackets(i1) - Packets(i1);
    if abs(NPackets(i1)) == 0
        Status(i1) = 0;
    else
        Status(i1) = 1;
    end
end
end