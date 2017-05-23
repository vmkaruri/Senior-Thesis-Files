function determine_black_vector(audio_file, Fs, word_size, word_no) 
    % 0.3 secs works good for 'C:\Users\VonBass\Desktop\timit\test\dr1\felc0\sa1.wav'
    %[audio_signal, Fs] = audioread(audio_file); % read in the input signal
    audio_signal = audio_file;
    window_size = ceil(word_size * Fs); % in seconds
    hop_size = ceil(word_size * Fs); % in seconds
    num_blocks = ceil(length(audio_signal)/hop_size);
    
    test_scaling(audio_signal, Fs, signal_size)
    
    start_index = max(floor(window_size * (word_no - 1)), 1);
    end_index = min(ceil(window_size * word_no), length(audio_signal));
    audio_time_slice = audio_signal(start_index : end_index, :);
    
    
    % sound(audio_time_slice, Fs);
end