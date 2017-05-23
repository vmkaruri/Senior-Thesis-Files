function zcr_ = get_zero_crossings(audio_file)
    % assume the window size is 2 seconds
    % there are overlaps in windowing
    % assume the step of shif is 1 second

    [wav, fs] = audioread(audio_file);
    wav = wav / max(max(wav));
    window_length = 2 * fs;
    step = 1 * fs; 
    frame_num = floor((length(wav)-window_length)/step) + 1;

    zcr_ = zeros(frame_num, 1);
    wav_window2 = zeros(window_length, 1);
    pos = 1;

    for i=1:frame_num
       wav_window = wav(pos:pos + window_length-1);
       wav_window2(2:end) = wav_window(1:end-1);
       zcr_(i) = 1/2 * sum(abs(sign(wav_window) - ...
           sign(wav_window2))) * fs / window_length;
       pos = pos + step;
    end
end