function train_test_mfccs() 
    load('train_mfccs.mat');
    train_data = features_matrix_new;
    train_values = value_vector;
    
    load('test_mfccs.mat');
    test_data = features_matrix_new;
    test_values = value_vector;
    
    nearest_neighbour(train_data, train_values, test_data, test_values)
end