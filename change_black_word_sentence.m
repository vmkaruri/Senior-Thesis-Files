function change_black_word_sentence(word_position, word_count)
     
    load('sentence_features_train.mat');
    values = [];
    count = 1;
    numzeros = 0;
    new_sentence_features = [];
    
    for i = 1 : size(sentence_features, 1) % no. of entries in the traning set
        if count == word_position
            values = [values ; 1];
            new_sentence_features = [new_sentence_features ; sentence_features(i, 1 : (size(sentence_features, 2) - 1))];
        elseif count ~= word_position && numzeros < 250
            values = [values ; 0];
            new_sentence_features = [new_sentence_features ; sentence_features(i, 1 : (size(sentence_features, 2) - 1))];
            numzeros = numzeros + 1;
        end
        
        count = mod(count + 1, word_count + 1);
        if count == 0
            count = 1;
        end
    end
    
    training_set = [new_sentence_features values];
    %training_set = new_sentence_features;
    save(strcat('sentence_features_', int2str(word_position), '.mat'), 'training_set');
    [svm_struct] = train_svm(training_set);
    load('sentence_features_test.mat');
    [confusion_matrix_svm] = test_svm(svm_struct, sentence_features);
    save('sentence_confusion_', int2str(word_position), '.mat'), 'confusion_matrix_svm');
        
end
