function separate_true_false()
    load('training.mat');
    
    dataset = [features values];
    extract_target_values(dataset, 1, 'training_true.mat');
    extract_target_values(dataset, 0, 'training_false.mat');
    
    load('testing_features_values1.mat');
    
    dataset = [features values1];
    extract_target_values(dataset, 1, 'testing_true.mat');
    extract_target_values(dataset, 0, 'testing_false.mat');
end