function [mfccs] = extract_mfccs(audio_file, Fs)
    % [audio_signal, Fs] = audioread(audio_file); % read in the input signal
    audio_signal = audio_file;
    mfccs = melfcc(audio_signal, Fs, 'wintime', 0.3, 'hoptime', 0.3);
    mfccs = mfccs';
end

% 2.5, 1.25