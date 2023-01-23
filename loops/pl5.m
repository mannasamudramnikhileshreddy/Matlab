sum=0;
n=input('Tell me a number \n');
count=input('no. of digits\n')
x=n;
if (count>=2)
    for i=1:n
        if (n>0)
            r=round(rem(n,10)-0.1);
            sum=sum*10+r;
            n=round((n/10)-0.1);
            c=x-sum;
        else
            break
        end
    end
    if (c==0)
        disp('the no. is palindrome');
    else
        disp('the no. is not palindrome');
    end
else
    disp('the no. is not palindrome');
end