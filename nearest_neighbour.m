function [confusion_matrix_knn] = nearest_neighbour(train_data, train_values, test_data, test_values)
    model = fitcknn(train_data, train_values);
    predicted_values = predict(model, test_data);
    
    TP = 0;
    TN = 0;
    FP = 0;
    FN = 0;
    
    for i = 1 : size(test_values, 1)
        if test_values(i, 1) == predicted_values(i, 1)
            if predicted_values(i, 1) == 1
                TP = TP + 1;
            else
                TN = TN + 1;
            end
        else
            if predicted_values(i, 1) == 1
                FP = FP + 1;
            else
                FN = FN + 1;
            end
        end
    end
    
    confusion_matrix_knn = [TP FN ; FP TN]
    
end