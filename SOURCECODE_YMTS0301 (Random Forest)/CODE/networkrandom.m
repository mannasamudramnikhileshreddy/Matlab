function YPred=networkrandom(as);
load featur.mat
load fet.mat
BaggedEnsemble = generic_random_forests(featur,fet,82,'classification');
YPred =predict(BaggedEnsemble,as);