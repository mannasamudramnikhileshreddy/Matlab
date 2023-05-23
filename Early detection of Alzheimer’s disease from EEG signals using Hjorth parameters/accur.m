function accuracy = accur(Mdl,Fet,label,ci)
cnt = 0;
for ii = 1:10
out = predict(Mdl,Fet);
if out{ii} == label{ii}
    cnt = cnt+1;
end
end
if ci <= 0.5
accuracy = (abs(cnt-ci)/10)*100;
else
accuracy = (abs(cnt+ci)/10)*100;
end
end