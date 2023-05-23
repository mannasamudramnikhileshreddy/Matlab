load('bermaxwith')
load('ber4maxwith')
load('bermax')
load('ber4max')
load('snrwith')
load('snr4with')
load('snrwithout')
load('snr4without')
figure;
plot(sort(bermaxwith),'r');
title('Bit Errorrate vs Channel Bandwidth');
hold on;
plot(sort(ber4maxwith),'b');
figure;
plot(flip(sort(ber4max(1:4))),'r');
title('snr vs Channel Bandwidth');
hold on;
plot(flip(sort(bermax(1:4))),'b');
figure;
plot(sort(snrwith),'b');
title('Bit Errorrate vs Channel Bandwidth');
hold on;
plot(sort(snr4with),'r');
figure;
plot(flip(sort(snrwithout)),'b');
title('snr vs Channel Bandwidth');
hold on;
plot(flip(sort(snr4without)),'r');