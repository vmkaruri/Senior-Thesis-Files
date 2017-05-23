function [svm_struct, features_no] = train_svm(features)
     %load('training.mat');
     %[features, values] = svmextract('C:\Users\VonBass\Desktop\UCSB\isolated_digits_ti_test_endpt\MAN');
     %save('training.mat', 'features', 'values');
     %svm_struct = svmtrain(features, values, 'kernel_function', 'rbf');
     svm_struct = svmtrain(features(:, 1 : (size(features, 2) - 1)), features(:, size(features, 2)), 'kernel_function', 'linear'); % quadratic works best for one word, polynomial for sentences
     features_no = size(features, 2);  
end