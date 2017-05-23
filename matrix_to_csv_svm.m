function matrix_to_csv_svm(train_true_no, train_false_no)

    %load('training.mat');  
    
    % train csv
    load('training_true.mat');
    training_true = result_matrix(1 : train_true_no, :);
    load('training_false.mat');
    training_false = result_matrix(1 : train_false_no, :);
    features_values = [training_true ; training_false];
  
    train_matrix = [[1 : size(features_values, 2)] ; features_values];
    csvwrite('train_svm.csv', train_matrix);
  
    % test csv
    load('testing_features_values1.mat');
    test_matrix = [[1 : size(features, 2)] ; features];
    csvwrite('test_svm.csv', [test_matrix [size(features, 2) + 1 ; values1]]);
end