a=randi([-10,10],1,10)
for i=1:length(a)
    if a(i)>0
      formatspec='a is thenumber %d which is greater than o'
      fprintf('formatspec',a(i))
    elseif a(i)<0
      formatspec='a is thenumber %d which is less than than o'
      fprintf('formatspec',a(i))

    else
      formatspec='a is thenumber %d which is equal than o'
      fprintf('formatspec',a(i))
    end
end