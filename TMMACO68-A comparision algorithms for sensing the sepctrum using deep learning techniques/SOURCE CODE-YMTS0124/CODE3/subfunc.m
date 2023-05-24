function [ coeff_daq ] = subfunc( h0,x,y )

alpha = 0.75;

n=size(x,2);
init_daq=randperm(n);
N = round(n*alpha);
in_dax = init_daq(1:N);
in_day = init_daq(N+1:n);

d=max(x')-min(x');
coeff_daq= fmincon(@(coeff_daq) recurs_fun(coeff_daq,x(:,in_dax),y(:,in_dax),x(:,in_day),y(:,in_day)),h0,[],[],[],[],d/20,2*d );
clc

end

