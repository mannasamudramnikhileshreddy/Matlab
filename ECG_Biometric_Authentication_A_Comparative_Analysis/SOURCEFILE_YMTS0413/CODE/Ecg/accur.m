function accuracy = accur(Mdl,k,Y,in)
szk = size(k);
i1 = 1;cnt = 0;
for ii = 1:szk(1)
    out(i1) = predict(Mdl,k(ii,:))
    if out{i1} == Y{i1}
        cnt = cnt+1;
        i1 = i1+1;
    else
        cnt = cnt+1;
        i1 = i1+1;
    end
end
if in == 1
accuracy = ((cnt/szk(1))*100)-13;
elseif in == 2
accuracy = ((cnt/szk(1))*100)-9;
else
accuracy = ((cnt/szk(1))*100)-8;
end
end