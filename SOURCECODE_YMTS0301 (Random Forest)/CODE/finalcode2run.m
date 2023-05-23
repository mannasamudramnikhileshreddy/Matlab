clc;
clear all;
close all;
warning off;

[file,path] = uigetfile('*.*','Select Signal');
Name1 = [path,file];
a1=load(Name1);
a1 = struct2cell(a1);
a1 = cell2mat(a1);

interval1=1;
x1 =(1:size(a1)) * interval1;
figure;
plot(x1',a1');ylabel('Amplitude');
title('ECG1 signal');

% % % Applying IIR butterworth filter
data1=a1;
fc = 30;
fs = 1000;
[b1,a1] = butter(1,fc/(fs/2));
k1=imfilter(data1,b1);
figure;
freqz(b1,a1)
title('IIR butterworth filter');

% % % Daubechies wavelet
wname1 = 'db4';
X1 = dbwavf(wname1);
feat1=imfilter(k1,X1);
features1=(feat1(1:100))';

% fet={'Enrolled','Enrolled','Enrolled','Enrolled','Enrolled','Not enrolled','Not enrolled','Not enrolled','Not enrolled','Not enrolled'}';


load featur.mat
load fet.mat

YPred=networkrandom(features1);

k=char(YPred);
msgbox(k);
%%

a=featur;
a=imresize(a,[256 256]);
u_bw_filename = im2bw(a);

b=features1;
b=imresize(b,[256 256]);
u_GT_filename = im2bw(b);

u_GT = [((u_GT_filename)) > 0 ];
u_bw = [((u_bw_filename)) > 0 ];

temp_obj_eval = objective_evaluation_core(u_bw, u_GT);

disp('Accuracy--');
disp(temp_obj_eval.Accuracy);

disp('Sensitivity--');
disp(temp_obj_eval.Sensitivity*100);
