function integrate_features()

    % train
    load('train_mfccs.mat');
    mfccs_features = features_matrix_new;
    mfccs_values = value_vector;
    
    load('train_mfccs_true.mat');
    mfcc_true = result_matrix;
    load('train_mfccs_false.mat');
    mfcc_false = result_matrix;
    
    load('train_compute_features.mat');
    compute_features = features_matrix_new;
    compute_values = value_vector;
    
    load('train_compute_features_true.mat');
    compute_true = result_matrix;
    load('train_compute_features_false.mat');
    compute_false = result_matrix;
    
    load('train_mir_others.mat');
    mir_features = features_matrix_new;
    mir_values = value_vector;
    
    load('train_mir_others_true.mat');
    mir_true = result_matrix;
    load('train_mir_others_false.mat');
    mir_false = result_matrix;
    
    load('train_mir_others_extra.mat');
    others_features = features_matrix_new;
    others_values = value_vector;
    
    load('train_mir_others_true_extra.mat');
    others_true = result_matrix;
    load('train_mir_others_false_extra.mat');
    others_false = result_matrix;
    
    len = min([size(mfccs_features, 1), size(compute_features, 1), size(mir_features, 1), size(others_features, 1)]);
    overall_train_features = [mfccs_features(1:len, :) compute_features(1:len, :) mir_features(1:len, :) others_features(1:len, :)];
    overall_train_values = mfccs_values(1:len, :);
    
    
    len = min([size(mfcc_true, 1), size(compute_true, 1), size(mir_true, 1), size(others_true, 1)]);
    overall_train_true = [mfcc_true(1:len, 1 : (size(mfcc_true, 2) - 1)) compute_true(1:len, 1 : (size(compute_true, 2) - 1)) mir_true(1:len, 1 : (size(mir_true, 2) - 1)) others_true(1:len, 1 : size(others_true, 2))];
    
    len = min([size(mfcc_false, 1), size(compute_false, 1), size(mir_false, 1), size(others_false, 1)]);
    overall_train_false = [mfcc_false(1:len, 1 : (size(mfcc_false, 2) - 1)) compute_false(1:len, 1 : (size(compute_false, 2) - 1)) mir_false(1:len, 1 : (size(mir_false, 2) - 1)) others_false(1:len, 1 : size(others_false, 2))];
    
    
   % test
    load('test_mfccs.mat');
    mfccs_features = features_matrix_new;
    mfccs_values = value_vector;
    
    load('test_mfccs_true.mat');
    mfcc_true = result_matrix;
    load('test_mfccs_false.mat');
    mfcc_false = result_matrix;
    
    load('test_compute_features.mat');
    compute_features = features_matrix_new;
    compute_values = value_vector;
    
    load('test_compute_features_true.mat');
    compute_true = result_matrix;
    load('test_compute_features_false.mat');
    compute_false = result_matrix;
    
    load('test_mir_others.mat');
    mir_features = features_matrix_new;
    mir_values = value_vector;
    
    load('test_mir_others_true.mat');
    mir_true = result_matrix;
    load('test_mir_others_false.mat');
    mir_false = result_matrix;
    
    load('test_mir_others_extra.mat');
    others_features = features_matrix_new;
    others_values = value_vector;
    
    load('test_mir_others_true_extra.mat');
    others_true = result_matrix;
    load('test_mir_others_false_extra.mat');
    others_false = result_matrix;
    
    len = min([size(mfccs_features, 1), size(compute_features, 1), size(mir_features, 1), size(others_features, 1)]);
    overall_test_features = [mfccs_features(1:len, :) compute_features(1:len, :) mir_features(1:len, :) others_features(1:len, :)];
    overall_test_values = mfccs_values(1:len, :);
    
    len = min([size(mfcc_true, 1), size(compute_true, 1), size(mir_true, 1), size(others_true, 1)]);
    overall_test_true = [mfcc_true(1:len, 1 : (size(mfcc_true, 2) - 1)) compute_true(1:len, 1 : (size(compute_true, 2) - 1)) mir_true(1:len, 1 : (size(mir_true, 2) - 1)) others_true(1:len, 1 : size(others_true, 2))];
    
    len = min([size(mfcc_false, 1), size(compute_false, 1), size(mir_false, 1), size(others_false, 1)]);
    overall_test_false = [mfcc_false(1:len, 1 : (size(mfcc_false, 2) - 1)) compute_false(1:len, 1 : (size(compute_false, 2) - 1)) mir_false(1:len, 1 : (size(mir_false, 2) - 1)) others_false(1:len, 1 : size(others_false, 2))];
    
    save('overall_feature_values.mat', 'overall_train_features', 'overall_train_true', 'overall_train_false', 'overall_train_values', 'overall_test_features','overall_test_true', 'overall_test_false', 'overall_test_values');
end