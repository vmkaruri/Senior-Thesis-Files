for i = 1 : 10
    confusion = train_test_sentences(i);
    save(strcat('sentence_confusion_rbf_', int2str(i), '.mat'), 'confusion');
end