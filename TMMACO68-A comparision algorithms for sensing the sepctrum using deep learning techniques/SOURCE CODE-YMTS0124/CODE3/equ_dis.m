function dis_matrix=equ_dis(a,b)
array_manip = sum(a.*a,1); bb = sum(b.*b,1); ab = a'*b; 
dis_matrix = abs(repmat(array_manip',[1 size(bb,2)]) + repmat(bb,[size(array_manip,2) 1]) - 2*ab);

