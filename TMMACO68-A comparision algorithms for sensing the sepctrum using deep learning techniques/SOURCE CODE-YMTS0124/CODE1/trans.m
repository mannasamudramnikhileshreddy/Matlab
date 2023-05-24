%%%%%%%%%%%%%% TRansmitter script %%%%%%%%%%%%%%%%%%%%%%
data=t_data((d1-1)*dataframelength+1:d1*dataframelength);
constlen=7;
codegen = [171 133];                            %%% Polynomial
trellis = poly2trellis(constlen, codegen,171);
codedata = convenc(data, trellis);              %%% Convolutional encoding  for redundancy
matrix=reshape(codedata,length(codedata)/2,2);
intdat = matintrlv(matrix,1,dataframelength);   %%% Interleaving coded data
in=intdat(:);

                                    %%% M = 4 for QPSK

    M=4;
y = pskmod(in,M);                              %%% QPSK modulation with "NULL offset"
pilt=1+1j;

lendata=length(y);

k=1;
pvec=[1:6:96];
pilot112=zeros(1,112);
for i=1:length(pvec)
    pilot112(pvec(i))=pilt;                      %%% pilots
end
for j=1:112;
    if(pilot112(j)==0+0j)
        pilot112(j)=y(k);
        k=k+1;
    end
end
padding16=zeros(1,16);
 data128 = [ pilot112 padding16];
data128=reshape(data128,1,128);
recv=[];
ifft_sig=(fftlength/sqrt(112))*ifft(fftshift(data128.')).';
recv(1:16)=ifft_sig(113:fftlength);              %%% Adding Cyclic Extension
recv(17:144)=ifft_sig;
