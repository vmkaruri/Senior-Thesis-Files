function y = spectral_rolloff (x,fs)
    % Spectral Rolloff - Frequency below which percentage (usually 85%)  
    % of signal energy is contained 
    %   x = windowed input signal
    %   fs = sample frequency

    percentage = 0.85;      % Rolloff percentage
    N1 = 1024;           % FFT length
    N2 = N1/2;           % Half FFT length

    fx = fft(x,N1);       % Compute FFT of x
    f = abs(fx(1:N2));    % Find magnitude spectrum

    summed = sum(f);      % Sum all magnitudes (total energy)
    rolloff = percentage * summed;   % Percentage of total energy
    energy = 0;
    i = 1;     % Index to frequency array 
    while(energy <= rolloff)
        energy = energy + f(i);
        i = i+1;
    end

    y = i / N1;   % Final index relative to max index
    y = y*fs;   % Convert to Hz

end