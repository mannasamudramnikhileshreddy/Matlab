function weight = sgdmethod(input,cor_out,weight)
alpha = 0.9;

N = length(input);

for k = 1:N
    tran_inp = input(k);
    d(k) = cor_out(k);
    wei_sum(k) = weight(k)*tran_inp;
    output(k) = sigmoid(wei_sum(k));
    
    error(k) = d(k) - output(k);
    delta(k) = output(k)*(1-output(k))*error(k);
    
    dweight(k) = alpha*delta(k)*tran_inp;
    
    weight(k) = weight(k)+dweight(k);
end
end