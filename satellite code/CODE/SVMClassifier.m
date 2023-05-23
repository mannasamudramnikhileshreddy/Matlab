function [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% Returns a trained classifier and its accuracy. This code recreates the
% classification model trained in Classification Learner app. Use the
% generated code to automate training the same model with new data, or to
% learn how to programmatically train models.
%
%  Input:
%      trainingData: A table containing the same predictor and response
%       columns as those imported into the app.
%
%  Output:
%      trainedClassifier: A struct containing the trained classifier. The
%       struct contains various fields with information about the trained
%       classifier.
%
%      trainedClassifier.predictFcn: A function to make predictions on new
%       data.
%
%      validationAccuracy: A double containing the accuracy in percent. In
%       the app, the History list displays this overall accuracy score for
%       each model.
%
% Use the code to train the model with new data. To retrain your
% classifier, call the function from the command line with your original
% data or new data as the input argument trainingData.
%
% For example, to retrain a classifier trained with the original data set
% T, enter:
%   [trainedClassifier, validationAccuracy] = trainClassifier(T)
%
% To make predictions with the returned 'trainedClassifier' on new data T2,
% use
%   yfit = trainedClassifier.predictFcn(T2)
%
% T2 must be a table containing at least the same predictor columns as used
% during training. For details, enter:
%   trainedClassifier.HowToPredict

% Auto-generated by MATLAB on 02-Dec-2021 14:47:52


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
% Split matrices in the input table into vectors
inputTable = [inputTable(:,setdiff(inputTable.Properties.VariableNames, {'features'})), array2table(table2array(inputTable(:,{'features'})), 'VariableNames', {'features_1', 'features_2', 'features_3', 'features_4', 'features_5', 'features_6', 'features_7', 'features_8', 'features_9', 'features_10', 'features_11', 'features_12', 'features_13', 'features_14', 'features_15', 'features_16', 'features_17', 'features_18', 'features_19', 'features_20', 'features_21', 'features_22', 'features_23', 'features_24', 'features_25', 'features_26', 'features_27', 'features_28', 'features_29', 'features_30', 'features_31', 'features_32', 'features_33', 'features_34', 'features_35', 'features_36', 'features_37', 'features_38', 'features_39', 'features_40', 'features_41', 'features_42', 'features_43', 'features_44', 'features_45', 'features_46', 'features_47', 'features_48', 'features_49', 'features_50', 'features_51', 'features_52', 'features_53', 'features_54', 'features_55', 'features_56', 'features_57', 'features_58', 'features_59', 'features_60', 'features_61', 'features_62', 'features_63', 'features_64', 'features_65', 'features_66', 'features_67', 'features_68', 'features_69', 'features_70', 'features_71', 'features_72', 'features_73', 'features_74', 'features_75', 'features_76', 'features_77', 'features_78', 'features_79', 'features_80', 'features_81', 'features_82', 'features_83', 'features_84', 'features_85', 'features_86', 'features_87', 'features_88', 'features_89', 'features_90', 'features_91', 'features_92', 'features_93', 'features_94', 'features_95', 'features_96', 'features_97', 'features_98', 'features_99', 'features_100', 'features_101', 'features_102', 'features_103', 'features_104', 'features_105', 'features_106', 'features_107', 'features_108', 'features_109', 'features_110', 'features_111', 'features_112', 'features_113', 'features_114', 'features_115', 'features_116', 'features_117', 'features_118', 'features_119', 'features_120', 'features_121', 'features_122', 'features_123', 'features_124', 'features_125', 'features_126', 'features_127', 'features_128', 'features_129', 'features_130', 'features_131', 'features_132', 'features_133', 'features_134', 'features_135', 'features_136', 'features_137', 'features_138', 'features_139', 'features_140', 'features_141', 'features_142', 'features_143', 'features_144', 'features_145', 'features_146', 'features_147', 'features_148', 'features_149', 'features_150', 'features_151', 'features_152', 'features_153', 'features_154', 'features_155', 'features_156', 'features_157', 'features_158', 'features_159', 'features_160', 'features_161', 'features_162', 'features_163', 'features_164', 'features_165', 'features_166', 'features_167', 'features_168', 'features_169', 'features_170', 'features_171', 'features_172', 'features_173', 'features_174', 'features_175', 'features_176', 'features_177'})];

predictorNames = {'features_1', 'features_2', 'features_3', 'features_4', 'features_5', 'features_6', 'features_7', 'features_8', 'features_9', 'features_10', 'features_11', 'features_12', 'features_13', 'features_14', 'features_15', 'features_16', 'features_17', 'features_18', 'features_19', 'features_20', 'features_21', 'features_22', 'features_23', 'features_24', 'features_25', 'features_26', 'features_27', 'features_28', 'features_29', 'features_30', 'features_31', 'features_32', 'features_33', 'features_34', 'features_35', 'features_36', 'features_37', 'features_38', 'features_39', 'features_40', 'features_41', 'features_42', 'features_43', 'features_44', 'features_45', 'features_46', 'features_47', 'features_48', 'features_49', 'features_50', 'features_51', 'features_52', 'features_53', 'features_54', 'features_55', 'features_56', 'features_57', 'features_58', 'features_59', 'features_60', 'features_61', 'features_62', 'features_63', 'features_64', 'features_65', 'features_66', 'features_67', 'features_68', 'features_69', 'features_70', 'features_71', 'features_72', 'features_73', 'features_74', 'features_75', 'features_76', 'features_77', 'features_78', 'features_79', 'features_80', 'features_81', 'features_82', 'features_83', 'features_84', 'features_85', 'features_86', 'features_87', 'features_88', 'features_89', 'features_90', 'features_91', 'features_92', 'features_93', 'features_94', 'features_95', 'features_96', 'features_97', 'features_98', 'features_99', 'features_100', 'features_101', 'features_102', 'features_103', 'features_104', 'features_105', 'features_106', 'features_107', 'features_108', 'features_109', 'features_110', 'features_111', 'features_112', 'features_113', 'features_114', 'features_115', 'features_116', 'features_117', 'features_118', 'features_119', 'features_120', 'features_121', 'features_122', 'features_123', 'features_124', 'features_125', 'features_126', 'features_127', 'features_128', 'features_129', 'features_130', 'features_131', 'features_132', 'features_133', 'features_134', 'features_135', 'features_136', 'features_137', 'features_138', 'features_139', 'features_140', 'features_141', 'features_142', 'features_143', 'features_144', 'features_145', 'features_146', 'features_147', 'features_148', 'features_149', 'features_150', 'features_151', 'features_152', 'features_153', 'features_154', 'features_155', 'features_156', 'features_157', 'features_158', 'features_159', 'features_160', 'features_161', 'features_162', 'features_163', 'features_164', 'features_165', 'features_166', 'features_167', 'features_168', 'features_169', 'features_170', 'features_171', 'features_172', 'features_173', 'features_174', 'features_175', 'features_176', 'features_177'};
predictors = inputTable(:, predictorNames);
response = inputTable.labels;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Apply a PCA to the predictor matrix.
% Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
isCategoricalPredictorBeforePCA = isCategoricalPredictor;
numericPredictors = predictors(:, ~isCategoricalPredictor);
numericPredictors = table2array(varfun(@double, numericPredictors));
% 'inf' values have to be treated as missing data for PCA.
numericPredictors(isinf(numericPredictors)) = NaN;
numComponentsToKeep = min(size(numericPredictors,2), 20);
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
    numericPredictors, ...
    'NumComponents', numComponentsToKeep);
predictors = [array2table(pcaScores(:,:)), predictors(:, isCategoricalPredictor)];
isCategoricalPredictor = [false(1,numComponentsToKeep), true(1,sum(isCategoricalPredictor))];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
template = templateSVM(...
    'KernelFunction', 'gaussian', ...
    'PolynomialOrder', [], ...
    'KernelScale', 3.3, ...
    'BoxConstraint', 1, ...
    'Standardize', true);
classificationSVM = fitcecoc(...
    predictors, ...
    response, ...
    'Learners', template, ...
    'Coding', 'onevsone', ...
    'ClassNames', {'Football'; 'Human Face'; 'Rugby Ball'; 'Satellite Image'});

% Create the result struct with predict function
splitMatricesInTableFcn = @(t) [t(:,setdiff(t.Properties.VariableNames, {'features'})), array2table(table2array(t(:,{'features'})), 'VariableNames', {'features_1', 'features_2', 'features_3', 'features_4', 'features_5', 'features_6', 'features_7', 'features_8', 'features_9', 'features_10', 'features_11', 'features_12', 'features_13', 'features_14', 'features_15', 'features_16', 'features_17', 'features_18', 'features_19', 'features_20', 'features_21', 'features_22', 'features_23', 'features_24', 'features_25', 'features_26', 'features_27', 'features_28', 'features_29', 'features_30', 'features_31', 'features_32', 'features_33', 'features_34', 'features_35', 'features_36', 'features_37', 'features_38', 'features_39', 'features_40', 'features_41', 'features_42', 'features_43', 'features_44', 'features_45', 'features_46', 'features_47', 'features_48', 'features_49', 'features_50', 'features_51', 'features_52', 'features_53', 'features_54', 'features_55', 'features_56', 'features_57', 'features_58', 'features_59', 'features_60', 'features_61', 'features_62', 'features_63', 'features_64', 'features_65', 'features_66', 'features_67', 'features_68', 'features_69', 'features_70', 'features_71', 'features_72', 'features_73', 'features_74', 'features_75', 'features_76', 'features_77', 'features_78', 'features_79', 'features_80', 'features_81', 'features_82', 'features_83', 'features_84', 'features_85', 'features_86', 'features_87', 'features_88', 'features_89', 'features_90', 'features_91', 'features_92', 'features_93', 'features_94', 'features_95', 'features_96', 'features_97', 'features_98', 'features_99', 'features_100', 'features_101', 'features_102', 'features_103', 'features_104', 'features_105', 'features_106', 'features_107', 'features_108', 'features_109', 'features_110', 'features_111', 'features_112', 'features_113', 'features_114', 'features_115', 'features_116', 'features_117', 'features_118', 'features_119', 'features_120', 'features_121', 'features_122', 'features_123', 'features_124', 'features_125', 'features_126', 'features_127', 'features_128', 'features_129', 'features_130', 'features_131', 'features_132', 'features_133', 'features_134', 'features_135', 'features_136', 'features_137', 'features_138', 'features_139', 'features_140', 'features_141', 'features_142', 'features_143', 'features_144', 'features_145', 'features_146', 'features_147', 'features_148', 'features_149', 'features_150', 'features_151', 'features_152', 'features_153', 'features_154', 'features_155', 'features_156', 'features_157', 'features_158', 'features_159', 'features_160', 'features_161', 'features_162', 'features_163', 'features_164', 'features_165', 'features_166', 'features_167', 'features_168', 'features_169', 'features_170', 'features_171', 'features_172', 'features_173', 'features_174', 'features_175', 'features_176', 'features_177'})];
extractPredictorsFromTableFcn = @(t) t(:, predictorNames);
predictorExtractionFcn = @(x) extractPredictorsFromTableFcn(splitMatricesInTableFcn(x));
pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x(:, ~isCategoricalPredictorBeforePCA))) - pcaCenters) * pcaCoefficients), x(:,isCategoricalPredictorBeforePCA) ];
svmPredictFcn = @(x) predict(classificationSVM, x);
trainedClassifier.predictFcn = @(x) svmPredictFcn(pcaTransformationFcn(predictorExtractionFcn(x)));

% Add additional fields to the result struct
trainedClassifier.RequiredVariables = {'features'};
trainedClassifier.PCACenters = pcaCenters;
trainedClassifier.PCACoefficients = pcaCoefficients;
trainedClassifier.ClassificationSVM = classificationSVM;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2020a.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
% Split matrices in the input table into vectors
inputTable = [inputTable(:,setdiff(inputTable.Properties.VariableNames, {'features'})), array2table(table2array(inputTable(:,{'features'})), 'VariableNames', {'features_1', 'features_2', 'features_3', 'features_4', 'features_5', 'features_6', 'features_7', 'features_8', 'features_9', 'features_10', 'features_11', 'features_12', 'features_13', 'features_14', 'features_15', 'features_16', 'features_17', 'features_18', 'features_19', 'features_20', 'features_21', 'features_22', 'features_23', 'features_24', 'features_25', 'features_26', 'features_27', 'features_28', 'features_29', 'features_30', 'features_31', 'features_32', 'features_33', 'features_34', 'features_35', 'features_36', 'features_37', 'features_38', 'features_39', 'features_40', 'features_41', 'features_42', 'features_43', 'features_44', 'features_45', 'features_46', 'features_47', 'features_48', 'features_49', 'features_50', 'features_51', 'features_52', 'features_53', 'features_54', 'features_55', 'features_56', 'features_57', 'features_58', 'features_59', 'features_60', 'features_61', 'features_62', 'features_63', 'features_64', 'features_65', 'features_66', 'features_67', 'features_68', 'features_69', 'features_70', 'features_71', 'features_72', 'features_73', 'features_74', 'features_75', 'features_76', 'features_77', 'features_78', 'features_79', 'features_80', 'features_81', 'features_82', 'features_83', 'features_84', 'features_85', 'features_86', 'features_87', 'features_88', 'features_89', 'features_90', 'features_91', 'features_92', 'features_93', 'features_94', 'features_95', 'features_96', 'features_97', 'features_98', 'features_99', 'features_100', 'features_101', 'features_102', 'features_103', 'features_104', 'features_105', 'features_106', 'features_107', 'features_108', 'features_109', 'features_110', 'features_111', 'features_112', 'features_113', 'features_114', 'features_115', 'features_116', 'features_117', 'features_118', 'features_119', 'features_120', 'features_121', 'features_122', 'features_123', 'features_124', 'features_125', 'features_126', 'features_127', 'features_128', 'features_129', 'features_130', 'features_131', 'features_132', 'features_133', 'features_134', 'features_135', 'features_136', 'features_137', 'features_138', 'features_139', 'features_140', 'features_141', 'features_142', 'features_143', 'features_144', 'features_145', 'features_146', 'features_147', 'features_148', 'features_149', 'features_150', 'features_151', 'features_152', 'features_153', 'features_154', 'features_155', 'features_156', 'features_157', 'features_158', 'features_159', 'features_160', 'features_161', 'features_162', 'features_163', 'features_164', 'features_165', 'features_166', 'features_167', 'features_168', 'features_169', 'features_170', 'features_171', 'features_172', 'features_173', 'features_174', 'features_175', 'features_176', 'features_177'})];

predictorNames = {'features_1', 'features_2', 'features_3', 'features_4', 'features_5', 'features_6', 'features_7', 'features_8', 'features_9', 'features_10', 'features_11', 'features_12', 'features_13', 'features_14', 'features_15', 'features_16', 'features_17', 'features_18', 'features_19', 'features_20', 'features_21', 'features_22', 'features_23', 'features_24', 'features_25', 'features_26', 'features_27', 'features_28', 'features_29', 'features_30', 'features_31', 'features_32', 'features_33', 'features_34', 'features_35', 'features_36', 'features_37', 'features_38', 'features_39', 'features_40', 'features_41', 'features_42', 'features_43', 'features_44', 'features_45', 'features_46', 'features_47', 'features_48', 'features_49', 'features_50', 'features_51', 'features_52', 'features_53', 'features_54', 'features_55', 'features_56', 'features_57', 'features_58', 'features_59', 'features_60', 'features_61', 'features_62', 'features_63', 'features_64', 'features_65', 'features_66', 'features_67', 'features_68', 'features_69', 'features_70', 'features_71', 'features_72', 'features_73', 'features_74', 'features_75', 'features_76', 'features_77', 'features_78', 'features_79', 'features_80', 'features_81', 'features_82', 'features_83', 'features_84', 'features_85', 'features_86', 'features_87', 'features_88', 'features_89', 'features_90', 'features_91', 'features_92', 'features_93', 'features_94', 'features_95', 'features_96', 'features_97', 'features_98', 'features_99', 'features_100', 'features_101', 'features_102', 'features_103', 'features_104', 'features_105', 'features_106', 'features_107', 'features_108', 'features_109', 'features_110', 'features_111', 'features_112', 'features_113', 'features_114', 'features_115', 'features_116', 'features_117', 'features_118', 'features_119', 'features_120', 'features_121', 'features_122', 'features_123', 'features_124', 'features_125', 'features_126', 'features_127', 'features_128', 'features_129', 'features_130', 'features_131', 'features_132', 'features_133', 'features_134', 'features_135', 'features_136', 'features_137', 'features_138', 'features_139', 'features_140', 'features_141', 'features_142', 'features_143', 'features_144', 'features_145', 'features_146', 'features_147', 'features_148', 'features_149', 'features_150', 'features_151', 'features_152', 'features_153', 'features_154', 'features_155', 'features_156', 'features_157', 'features_158', 'features_159', 'features_160', 'features_161', 'features_162', 'features_163', 'features_164', 'features_165', 'features_166', 'features_167', 'features_168', 'features_169', 'features_170', 'features_171', 'features_172', 'features_173', 'features_174', 'features_175', 'features_176', 'features_177'};
predictors = inputTable(:, predictorNames);
response = inputTable.labels;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

validationPredictFcn = @(x) svmPredictFcn(pcaTransformationFcn(x));

% Compute resubstitution predictions
[validationPredictions, validationScores] = validationPredictFcn(predictors);

% Compute validation accuracy
correctPredictions = strcmp( strtrim(validationPredictions), strtrim(response));
isMissing = cellfun(@(x) all(isspace(x)), response, 'UniformOutput', true);
correctPredictions = correctPredictions(~isMissing);
validationAccuracy = sum(correctPredictions)/length(correctPredictions);
