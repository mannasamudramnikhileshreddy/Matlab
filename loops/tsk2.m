a=[1,1,2,2,3,3,3,3]
count1=0
count2=0
count3=0
for i=1:length(a)
        if a(i)==1
            count1=count1+1
            disp(a())   
        elseif(a(i)==2)
            count2=count2+1
        else
             count3=count3+1
        end
end