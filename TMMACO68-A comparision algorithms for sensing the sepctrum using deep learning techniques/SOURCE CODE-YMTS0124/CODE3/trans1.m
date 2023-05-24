%%%%%%%%%%%%%% TRansmitter script %%%%%%%%%%%%%%%%%%%%%%
load pr3
data11=im2bw(pr3);
constlen1=7;
codegen1 = [171 133];                            %%% Polynomial
trellis1 = poly2trellis(constlen1, codegen1,171);
codedata1 = convenc(data11, trellis1);              %%% Convolutional encoding  for redundancy
matrix1=reshape(codedata1,length(codedata1)/2,2);
intdat1 = matintrlv(double(matrix1),1,1000);   %%% Interleaving coded data
in1=intdat1(:);

                                    %%% M = 4 for QPSK

M=4;
y1 = pskmod(in1,M);                              %%% QPSK modulation with "NULL offset"
pilt1=1+1j;



k=1;
pvec1=[1:6:96];
pilot112=zeros(1,112);
for i=1:length(pvec1)
    pilot112(pvec1(i))=pilt1;                      %%% pilots
end
for j=1:112;
    if(pilot112(j)==0+0j)
        pilot112(j)=y1(k);
        k=k+1;
    end
end
padding16=zeros(1,16);
 data1281 = [ pilot112 padding16];
data1281=reshape(data1281,1,128);
recv1=[];
ifft_sig1=(fftlength/sqrt(112))*ifft(fftshift(data1281.')).';
recv1(1:16)=ifft_sig1(113:fftlength);              %%% Adding Cyclic Extension
recv1(17:144)=ifft_sig1;
save recv1