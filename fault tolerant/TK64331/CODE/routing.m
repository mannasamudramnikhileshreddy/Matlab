function packets = routing(CHs,gridtype,GI,packets)
for i = 1:length(CHs)
    for i1 = 1:length(GI)
        if i ~= 1
            if GI(i) == GI(8)
                packets = packets-1;
                break;
            end
            if (gridtype(CHs(i)) > 8) && (gridtype(CHs(i)) <= 12)
                for i2 = 1:7
                    if GI(i-i2) == GI(8)
                        packets = packets-1;
                        break;
                    end
                end
            end
            if (gridtype(CHs(i)) > 1) && (gridtype(CHs(i)) <= 6)
                for i3 = 1:7
                    if i ~= i3
                        if GI(abs(i-i3)) == GI(8)
                            packets = packets-1;
                            break;
                        end
                    end
                end
                if imag(GI(i)) == imag(GI(8))
                    packets = packets-1;
                end
            end
            if (gridtype(CHs(i)) > 13) && (gridtype(CHs(i)) <= 18)
                for i4 = 1:7
                    if GI(i-i4) == GI(8) 
                        packets = packets-1;
                        break;
                    end
                end
                if imag(GI(i)) == imag(GI(8))
                    packets = packets-1;
                end
            end
            if (gridtype(CHs(i)) == 1) || (gridtype(CHs(i)) == 7) || (gridtype(CHs(i)) == 13)
                i5 = 1;
                while(i5<=5)
                    if GI(i+i5) == GI(8)
                        packets = packets-1;
                        break;
                    end
                    i5 = i5+1;
                end
                if (gridtype(CHs(i)) == 1) || (gridtype(CHs(i)) == 13)
                    if imag(GI(i)) == imag(GI(8))
                        packets = packets-1;
                    end
                end
            end
        end
    end
end
end