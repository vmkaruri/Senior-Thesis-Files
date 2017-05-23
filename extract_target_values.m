function [result_matrix] = extract_target_values(mat_file, target_value, file_name)

    load(mat_file);
    result_matrix = [];
    
    for i = 1 : size(value_vector, 1)
        if value_vector(i, 1) == target_value
            result_matrix = [result_matrix ; [features_matrix_new(i, :) value_vector(i, 1)]];
        end
    end
    
    save(file_name, 'result_matrix');
end