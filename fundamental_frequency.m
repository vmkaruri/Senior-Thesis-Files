function [freq] = fundamental_frequency(Fs, y)

    if (mod(length(y), 2) == 0)
        ydft = fft(y);
        freq = 0:Fs/length(y):Fs/2;
        ydft = ydft(1:length(y)/2+1);
        plot(freq,abs(ydft))
        [maxval,idx] = max(abs(ydft));
        freq = freq(idx);
    else
        ydft = fft(y);
        freq = 0:Fs/length(y):Fs/2;
        ydft = ydft(1:floor(length(y)/2)+1);
        [maxval,idx] = max(abs(ydft));
        freq = freq(idx);  %this is frequency corresponding to max value
    end
end