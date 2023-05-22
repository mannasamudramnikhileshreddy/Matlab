function accuracy = accur(x,l)
cnt = 0;
szl = size(l);
y = [1,1,1,1,1,0,0,0,0,0];
for i = 1:szl(1)
    out(i) = predict(x,l(i,:))
    if out{i} == 'NotEnrolled'
        ot = 0;
        if y(i) == 0
            cnt = cnt+1;
        end
    end
    if out{i} == 'Enrolled'
        ot = 1;
        if y(i) == 1
            cnt = cnt+1;
        end
    end
end
accuracy = (cnt/szl(1))*100;
end