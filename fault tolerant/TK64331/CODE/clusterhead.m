function [CH,er] = clusterhead(gridtype,n,rnum,Eo,Er,gt)
i4 = 1;
for i3 = 1:n
    if gridtype(i3) == gt
        gi(i4) = i3;
        er(i4) = Er(i3);
        i4 = i4+1;
    end
end
Per = 1/(i4-1);
ER = length(er);
t = thresholdequation(Per,1,er,ER);
nrnum = rnum(gi(1:(i4-1)));
i6 = 1;
for i5 = 1:(i4-1)
    if er(i5) >= t(i5)
        PCH(i6) = gi(i5);
        if i6 >= 2
            if rnum(gi(i5)) == min(nrnum)
                CH = gi(i5);
            end
        end
        i6 = i6+1;
    end
end          
end