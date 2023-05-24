
function pr_db=pow_ana(y) 
 
 
mea_rat=mean(abs(y).^2);        
thrs=max(abs(y).^2);    
p=thrs./mea_rat;   
pr_db=10.*log10(p); 