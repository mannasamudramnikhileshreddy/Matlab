function pow_rat=carrier(pr)


nloop=1000;
P=12000;
for p=1:1:P;
    pow_db=2+0.001*p;                  
    Pr(p)=(sum(pr>pow_db))/nloop;       
end
pow_rat=Pr;