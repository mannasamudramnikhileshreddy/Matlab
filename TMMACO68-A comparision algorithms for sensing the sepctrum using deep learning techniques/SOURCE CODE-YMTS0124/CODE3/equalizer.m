function rxout = equalizer(recv, ch_parameter)
a1=ch_parameter(1);
a2=ch_parameter(2);
d1=ch_parameter(3);
d2=ch_parameter(4);

copy1=zeros(size(recv));
for i=1+d1:length(recv)
    copy1(i)=a1*recv(i-d1);     
end
copy2=zeros(size(recv));
for i=1+d2:length(recv)
    copy2(i)=a2*recv(i-d2);
end
rxout=recv-copy1-copy2;