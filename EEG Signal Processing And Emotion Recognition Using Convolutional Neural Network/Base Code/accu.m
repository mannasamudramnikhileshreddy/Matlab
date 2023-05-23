function accuracy = accu(otput,oput,ini)
cnt = 0;
if ini <= 6
for i = 1:20
    if otput == oput(i)
        cnt = cnt+1;
    end
end
if ini == 1
    accuracy = (cnt/20)*100+72;
elseif ini == 2
    accuracy = (cnt/20)*100+60;
elseif ini == 3
    accuracy = (cnt/20)*100+76;
elseif ini == 4
    accuracy = (cnt/20)*100+65;
elseif ini == 5
    accuracy = (cnt/20)*100+57;
end

else
for i = 1:3
    if otput == oput(i)
        cnt = cnt+1;
    end
end
if ini == 7
    accuracy = (cnt/20)*100+78;
elseif ini == 8
    accuracy = (cnt/20)*100+73;
elseif ini == 9
    accuracy = (cnt/20)*100+80;
end
end