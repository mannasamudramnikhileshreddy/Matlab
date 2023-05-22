function Y = pred(inp,samp,out)
szinp = size(inp);Y = 0;
for i = 1:szinp(1)
    if samp(1:30) == inp(1:30,i)
        Y = out{i};
        break
    end
    if i == szinp(1)
        if Y == 0
            k = randi([1,szinp(1)],1,1);
            Y = out{k};
            break
        end
    end
end
end