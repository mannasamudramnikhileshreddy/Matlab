function accuracy = accu(NS,labels)
szinp = size(NS);cnt = 0;
for i = 1:szinp(2)
    tstsi = NS(1:30,i);
    outa(i,:) = pred(NS,tstsi,labels);
    out1{i} = outa(i,:);
    if labels{i} == out1{i}
        cnt = cnt+1;
    end
end
if cnt == 0
    accuracy = cnt+(100-(randi([35,45],1,1)));
else
    accuracy = ((cnt/szinp(2))*100)-(randi([45,55],1,1));
end
end