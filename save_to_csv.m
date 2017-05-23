files = dir('3s and rest test/*.waV');
% fopen('testing_features_list_twos.csv','w');
for file = files'
    audio_file = file.name
    
    try
        features = extract_audio_features(audio_file);
    catch
        continue;
    end
    
    dlmwrite('testing_features_list_ones_twos.csv',features,'-append')
end
