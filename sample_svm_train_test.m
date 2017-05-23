function [confusion_matrix] = sample_svm_train_test(train_false_no, test_true_no, test_false_no, feature_select, test_method, varargin)

    % train
    load('train_compute_features_true.mat');
    training_true = result_matrix;
    load('train_compute_features_false.mat');
    training_false = result_matrix(1 : train_false_no, :);
    training_set = [training_true ; training_false];
    features = training_set;
    
    % feature selection
    if length(varargin) > 1
     choose = varargin{1};
    end
    
    features_matrix_new = [];
    if feature_select
        for i = 1:length(choose)
           features_matrix_new =  [features_matrix_new features(:, choose(i))];
        end 
    else
        features_matrix_new = features;
    end
    
    %[svm_struct, features_no] = train_svm(features_matrix_new);
    
    % test
    load('test_compute_features_true.mat');
    testing_true = result_matrix;
    load('test_compute_features_false.mat');
    testing_false = result_matrix(1 : test_false_no, :);
    testing_set = [testing_true ; testing_false];
    %test_svm(svm_struct, features_no, testing_set);
    
    if test_method == 0
       confusion_matrix = nearest_neighbour(training_set(:, 1 : (size(features, 2) - 1)), training_set(:, size(features, 2)), testing_set(:, 1 : (size(features, 2) - 1)), testing_set(:, size(features, 2)));
    elseif test_method == 1
       svm_struct = train_svm(training_set);
       confusion_matrix = test_svm(svm_struct, size(training_set,2), testing_set);     
    end
end