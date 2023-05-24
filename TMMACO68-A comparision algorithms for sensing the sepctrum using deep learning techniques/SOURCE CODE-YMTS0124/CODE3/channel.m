function [chout chout1] = channel(recv,ch_parameter,snr)

a1=ch_parameter(1);
a2=ch_parameter(2);
d1=ch_parameter(3);
d2=ch_parameter(4);

%%%%%%%%%%%%%%%% multipath spectrum sensing Funtionality
%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%
copy1=zeros(size(recv));
for i=1+d1:length(recv)
    copy1(i)=a1*recv(i-d1);
end
copy2=zeros(size(recv));
for i=1+d2:length(recv)
    copy2(i)=a2*recv(i-d2);
end
recv=recv+copy1+copy2;

%%%%%%%%%% Channel is modelled as Rayleigh Multipath  Fadding + AWGN Noise combined %%%%%%%%%%
% Gaussian noise of unit variance, 0 mean
  Gau_nt = 1/sqrt(2)*[randn(1,144) + j*randn(1,144)];
  channel_out = sqrt(144/128)*recv + 10^(-snr/20)*Gau_nt;
 channel_out1 = sqrt(144/128)*recv + 10^(-snr+10/20)*Gau_nt;

chout= channel_out;
chout1= channel_out1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


