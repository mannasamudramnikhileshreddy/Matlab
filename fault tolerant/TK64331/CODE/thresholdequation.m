function T = thresholdequation(Per,r,er,ER)
for re = 1:ER
    T(re) = (Per/(1-(Per*(r*(mod(1,Per))))))*er(re);
end
end