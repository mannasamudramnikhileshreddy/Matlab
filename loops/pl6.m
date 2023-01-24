%sum of all the even number between the given range
a=input('enter the the range of the number:')
count=0;
for i=1:a
    if mod(i,2)==0
        count=count+i
    end
end
        