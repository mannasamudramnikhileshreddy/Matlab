a=input('enter the array in the []:');
n=input('enter the number to search in the array[]:');
count=0;
for i=1:length(a)
    for j=1
        if (a(i)== n(j))
            count = count+1;
        end
    end
end 