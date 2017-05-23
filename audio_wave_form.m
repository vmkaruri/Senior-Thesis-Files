function audio_wave_form(audio_file)

    [signal,fs] = audioread(audio_file);

%If signal is Nx2 (two columns), extract one of them
    signal = signal(:,1);

    % If you have the Signal Processing Toolbox, enter
    plot(psd(spectrum.periodogram,signal,'Fs',fs,'NFFT',length(signal)));

end