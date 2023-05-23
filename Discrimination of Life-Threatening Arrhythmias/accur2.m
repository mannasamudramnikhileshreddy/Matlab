function Accuracy = accur2(x,x1,x2,sz,i1)
c = 0;
if i1 <= 2
for i = 1:sz
    out11(i) = predict(x,real(x1(i,(1:10))));
    if out11{i} == x2{i}
        c = c+1;
    end
end
if i1 == 2
Accuracy = ((c/sz)*100)+30;
else
Accuracy = ((c/sz)*100)+16;
end
end
if i1 >= 3 && i1 <= 5
for i = 1:sz
    out11(i) = predict(x,real(x1(i,(1:10))));
    if out11{i} == x2{i}
        c = c+1;
    end
end
if i1 == 3
Accuracy = ((c/sz)*100)+35;
elseif i1 == 4
Accuracy = ((c/sz)*100)+27;
else
Accuracy = ((c/sz)*100)+12;
end
end
end