function change_black_word()
    for black_listed = 1:1
        % train
        %svmextract('C:\Users\VonBass\Desktop\UCSB\isolated_digits_ti_train_endpt\MAN',int2str(black_listed), 'train_mfccs.mat');
        %%svmextract('C:\Users\VonBass\Desktop\UCSB\isolated_digits_ti_train_endpt\MAN',int2str(black_listed), 'train_mir_others_extra.mat');
        % test
        %svmextract('C:\Users\VonBass\Desktop\UCSB\isolated_digits_ti_test_endpt\MAN', int2str(black_listed),'test_mfccs.mat');
        %%svmextract('C:\Users\VonBass\Desktop\UCSB\isolated_digits_ti_test_endpt\MAN', int2str(black_listed),'test_mir_others_extra.mat');

        % separate true and false values for train and test
        %%extract_target_values('train_mir_others_extra.mat', 1, 'train_mir_others_true_extra.mat');
        %%extract_target_values('train_mir_others_extra.mat', 0, 'train_mir_others_false_extra.mat');
        %%extract_target_values('test_mir_others_extra.mat', 1, 'test_mir_others_true_extra.mat');
        %%extract_target_values('test_mir_others_extra.mat', 0, 'test_mir_others_false_extra.mat');
        
        feature_select_matrix = [1:13 14 17 18 29 31 39 61]; % from weka and mfccs
        test_300 = sample_svm_train_test_overall(300, 70, 600, 1, 1, 0, feature_select_matrix);
        test_350 = sample_svm_train_test_overall(350, 70, 600, 1, 1, 0, feature_select_matrix);
        test_400 = sample_svm_train_test_overall(400, 70, 600, 1, 1, 0, feature_select_matrix);
        
        overall_accuracy_300 = (test_300(1,1) + test_300(2,2)) / (test_300(1,2) + test_300(2,1) + test_300(1,1) + test_300(2,2));
        overall_accuracy_350 = (test_350(1,1) + test_350(2,2)) / (test_350(1,2) + test_350(2,1) + test_350(1,1) + test_350(2,2));
        overall_accuracy_400 = (test_400(1,1) + test_400(2,2)) / (test_400(1,2) + test_400(2,1) + test_400(1,1) + test_400(2,2));
        
        positive_accuracy_300 = test_300(1,1) / (test_300(1,1) + test_300(1,2));
        negative_accuracy_300 = test_300(2,2) / (test_300(2,1) + test_300(2,2));
        
        positive_accuracy_350 = test_350(1,1) / (test_350(1,1) + test_350(1,2));
        negative_accuracy_350 = test_350(2,2) / (test_350(2,1) + test_350(2,2));
        
        positive_accuracy_400 = test_400(1,1) / (test_400(1,1) + test_400(1,2));
        negative_accuracy_400 = test_400(2,2) / (test_400(2,1) + test_400(2,2));
        
        mat_file_name = strcat('blacklist_', int2str(black_listed), '_overall_features.mat');
        save(mat_file_name, 'test_300', 'test_350', 'test_400', 'overall_accuracy_300', 'overall_accuracy_350', 'overall_accuracy_400', 'positive_accuracy_300', 'negative_accuracy_300', 'positive_accuracy_350', 'negative_accuracy_350', 'positive_accuracy_400', 'negative_accuracy_400');
    end
end
