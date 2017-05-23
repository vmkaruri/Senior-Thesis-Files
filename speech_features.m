function features = speech_features(sig, framesize, framestep, std_dev)
% Assumes 8kHz speech
%
% Features:
% 1 - formant frequency (Hz)
% 2 - confidence in formant frequency
% 3 - spectral entropy
%
% 4 - value of largest autocorrelation peak
% 5 - location of largest autocorrelation peak
% 6 - number of autocorelation peaks
%
% 7 - energy in frame
% 8 - time derivative of energy in frame
%
% 9:21 - Mel frequency cepstrum coefficients
%
% Spectral entropy can be used as a measure of how the energy in a voiced
% frame drops with frequency, because it tends to be roughly monotonic:
% lots of energy at low frequencies, less at higher.  An alternate would be
% to fit an exponential distribution, or even a gamma distribution.  I'll
% leave that as "future work."
sampling_rate = 16000;

num_frames = floor( (length(sig) - (framesize-framestep))/framestep);
features = zeros(21, num_frames);

%%%%%%%%%%  Spectrogram computation
% Use hanning window, to match Matlab's specgram() function.
h = hanning(framesize);
inv_length = 1/framesize;
nfft = framesize;

num_formant_fft = framesize*2;

% Minimum and maximum formant frequency, in Hz.
min_freq = 50;
max_freq = 500;

max_index = round(num_formant_fft/(min_freq / sampling_rate * framesize) + 1);
min_index = round(num_formant_fft/(max_freq / sampling_rate * framesize) + 1);


%%%%%%%%%%  Autocorrelation computation:
nacorr = framesize/2;
% compensate for frame effects: scale by 1/( 1/framesize -- framesize/framesize -- 1/framesize)
% compensating multiplier
comp = [framesize:-1:framesize-nacorr+1]/framesize;
comp = 1./comp';


%%%%%%%%%%  Energy computation
energy_deriv_filter = ((0:255) - 255/2)' / sum(abs((0:255) - 255/2));

for f = 1:num_frames
    % Grab frame and remove DC component
    frame = sig((f-1)*framestep+1:(f-1)*framestep+framesize);

    features(7, f) = sum(frame.^2) / framesize;

    frame = frame - sum(frame) / framesize;
    
    %%%%%%%%%%  Spectrogram
    % Compute normalized spectrogram
    tmp = fft(frame .* h,nfft); 
    spec = abs(tmp(1:(nfft/2)));
    normspec = spec / (sum(spec) + 1e-5);

    % Find the formant frequency
    big_fft = fft(sqrt(normspec), num_formant_fft);
    big_fft = real(big_fft(min_index:max_index));
    [features(2, f) maxloc] = max(big_fft);
    maxloc = maxloc + (min_index - 1);
    num_samples = num_formant_fft ./ (maxloc - 1);
    features(1, f) = num_samples .* (sampling_rate / framesize);

    % Compute the spectral entropy
    % Unscaled range is 0 to 5.5452.  Scale it to roughly 0..1.
    normspec(normspec < 1e-5) = 1e-5;
    features(3, f) = -sum(normspec .* log(normspec));
    if ~isfinite(features(3, f))
        % This becomes Inf on negotiation 109 somehow.
        % Actually, I now think that was memory corruption.
        error('Spectral entropy became infinite')
    end
    
    %%%%%%%%%%  Autocorrelation
    X = fft(frame,2*framesize);
    c = ifft(X.*conj(X));
    % Multiply by comp to compensate for frame effects
    acorr = real(c(1:nacorr)).*comp/(sum(frame.^2)+std_dev^2*framesize);
    
    % Compute the values of the autocorrelation peaks
    [peakvals peaklocs] = fast_find_acorr_peaks(acorr);
    [features(4, f) index] = max(peakvals);
    features(5, f) = peaklocs(index);
    features(6, f) = length(peakvals);
    
    
    %%%%%%%%%%  Energy
    
    features(8, f) = sum(frame.^2 .* energy_deriv_filter);

end

%%%%%%%%%%  MFCC
features(9:21, :) = mfcc(sig, sampling_rate, sampling_rate/framestep);
