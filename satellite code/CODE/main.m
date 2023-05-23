clc;
clear all;
close all;warning off;
close all hidden;
%%
[file, path] = uigetfile('*','Select an image');  %Open file selection dialog box
input_img = imread([path,file]);  %Read image from graphics file
figure;imshow(input_img);title('Input Image');  %Display image
channel1 = input_img(:,:,1);
features1 = extractLBPFeatures(channel1);
channel2 = input_img(:,:,2);
features2 = extractLBPFeatures(channel2);
channel3 = input_img(:,:,3);
features3 = extractLBPFeatures(channel3);
features = [features1, features2, features3];
figure; bar(features); title('LBP Features of three channels');
tablefeatures = table(features);
% labels = {'Rugby Ball';'Rugby Ball';'Rugby Ball';'Rugby Ball';'Rugby Ball';'Football';'Football';'Football';'Football';'Football';'Football';'Human Face';'Human Face';'Human Face';'Human Face';'Human Face';'Satellite Image';'Satellite Image';'Satellite Image';'Satellite Image';'Satellite Image';'Satellite Image';'Satellite Image';'Satellite Image'};
% tablelabels = table(labels);
% load features;
% load labels;
% TExpanded = [tablefeatures, tablelabels];
load TExpanded;
%%
[trainedClassifier, validationAccuracy] = SVMClassifier(TExpanded);
YPred = trainedClassifier.predictFcn (tablefeatures); 
fprintf('Output of SVM Classifier is: %s\n',char(YPred));
fprintf('Accuracy of SVM Classifier is: %f\n\n',validationAccuracy*100);
