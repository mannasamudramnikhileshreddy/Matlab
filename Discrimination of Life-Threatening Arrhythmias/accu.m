function Accuracy = accu(x,x1,x2,sz)
c = 0;
%sz1 = size(x1);
for i = 1:sz
    out11(i) = predict(x,real(x1(i,(1:10))));
    if out11{i} == x2{i}
        c = c+1;
    end
end
Accuracy = (c/sz)*100;
end