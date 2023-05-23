function [Xmin,Ymin,XminI,YminI,XminID,YminID] = Paths(xlocs,ylocs,CHBS)
i7 = 1;
for i3 = 1:length(CHBS)
    i5 = 1;
    for i4 = 1:length(CHBS)
        if i3 ~= i4
            dxlCHBS(i5) = abs(xlocs(CHBS(i3)) - xlocs(CHBS(i4)));
            dylCHBS(i5) = abs(ylocs(CHBS(i3)) - ylocs(CHBS(i4)));
            i5 = i5+1;
            if i5 == length(CHBS)
                Xmin(i3) = min(dxlCHBS);
                Ymin(i3) = min(dylCHBS);
                for i6 = 1:length(dxlCHBS)
                    if dxlCHBS(i6) == Xmin(i3)
                        XminI(i3) = i6;
                        XminID(i3) = CHBS(i6);
                    end
                    if dylCHBS(i6) == Ymin(i3)
                        YminI(i3) = i6;
                        YminID(i3) = CHBS(i6);
                    end
                end
            end
        end
    end
    i7 = i7+1;
end
end