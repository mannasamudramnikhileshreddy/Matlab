for z = 1:length(ugt)
    if (z >= 3 || z <= 6) || (z >= 9 || z <= 12) || (z >= 15 || z <= 18)
        PKTS2 = 0;
        if Status(z) ~= 0
            if Status(z+6) ~= 0
                Packets(z+6) = Packets(z+6)-1;
                PKTS2 = OPackets(z+6)-NPackets(z+6);
            elseif Status(z+6) == 0
            end
            if Status(z-6) ~= 0
                Packets(z-6) = Packets(z-6)-1;
                PKTS2 = OPackets(z-6)-NPackets(z-6)+PKTS2;
            elseif Status(z-6) == 0
            end
            if Status(z+1) ~= 0
                Packets(z+1) = Packets(z+1)-1;
                PKTS2 = OPackets(z+1)-NPackets(z+1)+PKTS2;
            elseif Status(z+1) == 0
            end
            if Status(z-1) ~= 0
                Packets(z-1) = Packets(z-1)-1;
                PKTS2 = OPackets(z-1)-NPackets(z-1)+PKTS2;
            elseif Status(z-1) == 0
            end
            Packets(z) = Packets(z)-1;
            PKTS2 = OPackets(z)-NPackets(z)+PKTS2;
            PKTS = PKTS2;
    else
    end
    end
end