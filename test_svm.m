% function [confusion_matrix_svm] = test_svm(svm_struct, training_col_size, features
function [confusion_matrix_svm] = test_svm(svm_struct, features)
% performs SVM classification on provided svm_struct and on the provided
% features matrix and returns a confusion matrix confusion_matrix_svm

    %[features, values1] = svmextract('C:\Users\VonBass\Desktop\UCSB\isolated_digits_ti_train_endpt\MAN');
    %save('testing_features_values1.mat', 'features','values1');
    %load('testing.mat');
%     if size(features, 2) < training_col_size    
%         diff = training_col_size - size(features, 2);
%         features = [features zeros(size(features, 1), diff)];
%     elseif size(features, 2) > training_col_size
%         features = features(: , 1 : training_col_size);
%     end
        
    %values2 = svmclassify(svm_struct, features);
    values1 = features(:, size(features, 2));
    values2 = svmclassify(svm_struct, features(:, 1 : (size(features, 2) - 1)));
    % save('testing_values2.mat1','values2');
    
    FP = 0;
    FN = 0;
    TP = 0;
    TN = 0;
    
    for i =  1 : size(values1, 1)
        if values1(i, 1) == values2(i, 1)
            if values2(i, 1) == 1
                TP = TP + 1;
            else
                TN = TN + 1;
            end
        else
            if values2(i, 1) == 1
                FP = FP + 1;
            else
                FN = FN + 1;
            end
        end        
    end
    
    confusion_matrix_svm = [TP FN ; FP TN]
    save('testing_confusion_matrix.mat', 'confusion_matrix_svm');
end
