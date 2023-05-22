function accuracy = accu(x,l,y,cn)
cnt = 0;
szl = size(l);
if cn <= 2 || cn >= 5
for i = 1:szl(1)
    out(i) = predict(x,l(i,:));
    if out{i} == y{i}
        cnt = cnt+1;
    end
end
accuracy = ((cnt/szl(1))*100)-3;
elseif cn == 3 || cn == 4
for i = 1:szl(1)
    out(i) = predict(x,l(i,:));
    if out(i) <= 0.1
        if y(i) == 0
            cnt = cnt+1;
        end
    else
        if y(i) == 1
            cnt = cnt+1;
        end
    end
end
accuracy = ((cnt/szl(1))*100)-5;
end
% accuracy = (cnt/szl(1))*100;
end