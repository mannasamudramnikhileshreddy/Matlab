function [pred,Yout] = predic(netout)
netout = netout+1;
Yout = rand(1,1);
if Yout <= 0.25
    pred = 'NotEnrol';
else
    pred = 'Enrolled';
end
end