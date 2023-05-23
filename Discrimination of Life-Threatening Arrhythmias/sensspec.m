function [sens,spec,TP,TN,FP,FN] = sensspec(y,y1,x2)
TP = 1;FP = 1;TN = 1;FN = 1;
for i = 1:length(x2)
    out11(i) = predict(y,real(y1(i,(1:10))));
    if out11{i} == x2{i}
        if out11{i} == x2{i}
           TP = TP+1;
        else
           FP = FP+1;
        end
    else
        if out11{i} == x2{i}
            FN = FN+1;
        else
            TN = TN+1;
        end
    end
end
sens = (TP/(TP+FN))*100;
spec = (TN/(TN+FP))*100;
end