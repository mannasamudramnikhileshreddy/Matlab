function [ daq_ou ] = recurs_fun( h,x,y,xt,yt )

ys=opt_reg(xt,x,y,h);
e=(yt-ys);
daq_ou = e*e';

end

