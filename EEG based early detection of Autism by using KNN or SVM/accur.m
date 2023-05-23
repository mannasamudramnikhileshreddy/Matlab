function accuracy = accur(tsl,tsf,Mdl)
sztf = size(tsf);
cnt = 0;
for in = 1:sztf(1)
    ypredn = predict(Mdl,tsf(in,:));
    if ypredn == tsl{in}
        cnt = cnt+1;
    end
end
accuracy = (cnt/sztf(1))*100;
end