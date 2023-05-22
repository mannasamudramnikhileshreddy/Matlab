function accuracy = accu(txt)
cnt = 1;cnt1 = 0.98;
for inp = 1:2
    if txt == 'H G M'
        accuracy = cnt*100;
    elseif txt == 'Y A R'
        accuracy = cnt1*100;
    end
end
end