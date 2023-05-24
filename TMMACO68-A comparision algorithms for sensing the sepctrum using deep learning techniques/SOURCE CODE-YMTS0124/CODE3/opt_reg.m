
function ys=opt_reg(xs,x,y,h)

K=equ_dis(diag(1./h)*x,diag(1./h)*xs);
K=exp(-K/2);
ys = (y*K) ./ (sum(K));
