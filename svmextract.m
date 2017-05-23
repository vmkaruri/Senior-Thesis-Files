function [features_matrix_new, value_vector] = svmextract(dir_path, black_listed, varargin)

    %// List all sub-directories under MainFolder recursively
    D = rdir(dir_path);             %// List of all sub-directories

    value_vector = [];
    features_matrix = [];
    features_matrix_new = [];
    raptured = 0;
    file_name = '';

    if length(varargin) > 0
        file_name = varargin{1};
    end
    
    for k = 1:length(D)

        %if k > 20
        %    break;
        %end

        currpath = D(k).name;                 %// Name of current directory

        %// Operate on all .waV files in currpath
        F = dir(fullfile(currpath, '*.waV')); %// Listing of all relevant files
        for m = 1:length(F)
            currfile = F(m).name;             %// Name of current file

            pattern = '._';
            if (strncmpi(pattern, currfile, 2))
                continue;
            %elseif ((~strncmpi('1', currfile, 1) && k > 10))
            %   continue;
            else            
                try
                    %features = extract_audio_features(strcat(currpath, '\', currfile));
                    features = extract_mfccs(strcat(currpath, '\', currfile));
                    % features = test_compute_feature(strcat(currpath, '\', currfile));
                    % features = extract_mir_others(strcat(currpath, '\', currfile));
                catch
                    continue;
                end
                            
               
                if ~raptured
                    %csvwrite('direct_train_svm.csv',features,'-append');  % write directly to file
                    
                    %if size(features, 2) < 13 % no. of features
                    %    len = 13 - size(features, 2);
                    %    padding = zeros(1, len);
                    %    features = [features padding];
                    %elseif size(features, 2) > 13
                    %    features = features(1 : 13);
                    %end
                    features_matrix = features;
                else
                    len = min(size(features, 2) , size(features_matrix, 2));
                    features = features(1 : len);
                    features_matrix = [features_matrix( : , 1 : len); features];
                    
                    %if size(features, 2) < 13 % no. of features
                    %    len = 13 - size(features, 2);
                    %    padding = zeros(1, len);
                    %    features = [features padding];                   
                    %elseif size(features, 2) > 13
                    %    features = features(1 : 13);
                    %end
                    
                    %features_matrix = [features_matrix ; features];
                end
                
                %if (strncmpi('1', currfile, 1) || strncmpi('2', currfile, 1) || strncmpi('3', currfile, 1))
                if (strncmpi(black_listed, currfile, 1))
                    if ~raptured
                        value_vector = 1;
                        raptured = 1;
                    else
                        value_vector = [value_vector; 1];
                    end
                    % display(strcat(currfile,' true'));
                else                   
                    if ~raptured 
                        value_vector = 0;
                        raptured = 1;
                    else
                        value_vector = [value_vector; 0];
                    end
                    % display(strcat(currfile, ' false'));
                end        
            end
             size(features_matrix, 1)
             size(value_vector, 1)
        end
    end
      % choose = [12,13,15,42,60,96,103,126];
      
      
      % feature selection
      %for i = 1:length(choose)
      % features_matrix_new =  [features_matrix_new features_matrix(:, choose(i))];
      %end
      
    features_matrix_new = features_matrix;
    
    if length(varargin) > 0
        save(file_name, 'features_matrix_new', 'value_vector');
    end
end
