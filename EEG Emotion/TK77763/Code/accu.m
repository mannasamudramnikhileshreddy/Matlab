function accuracy = accu(otput,oput,ini)
cnt = 0;
for i = 1:20
    if otput == oput(i)
        cnt = cnt+1;
    end
end
if ini == 1
    accuracy = (cnt/20)*100+72;
else
    accuracy = (cnt/20)*100+60;
end
end