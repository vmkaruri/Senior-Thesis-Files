function [onset,osr,D,tt,ff] = beat_onset(d,sr)
% [onset,osr,D,tt,ff] = beat_onset(d,sr)
%   Simple onset strength function for beat tracking.
%   d, sr define a waveform
%   or d is the name of a waveform file.
%   onset returns a sampled "onset strength" waveform
%   at a frame rate returned oin osr.
%   D returns the mel-spectrogram from which the onsets are
%   derived. tt labels the times of its columns, and ff the
%   frequencies of its rows.
% 2012-03-27 Dan Ellis dpwe@ee.columbia.edu

targetsr = 11025;

if ischar(d)
  % passed a file name - read it in
  fname = d;
  forcemono = 1;
  [d,sr] = audioread(fn,targetsr,forcemono);
else
  % passed actual waveform
  fname = '<noname>';
end

% Fix sampling rate if needed
if sr ~= targetsr
  d = resample(d,targetsr,sr);
end

% specgram: 256 bin @ 11kHz = 23 ms / 3 ms hop
swin = 256;
shop = 32;
% mel channels
nmel = 20;
% sample rate for specgram frames (granularity for rest of processing)
osr = targetsr/shop;

% Calculate spectrogram
D = specgram(d,swin,sr,swin,swin-shop);

% Calculate a time base for the spectrogram (to return)
tt = [0:size(D,2)-1]/osr;

% Construct db-magnitude-mel-spectrogram
[mlmx, ff] = fft2melmx(swin,sr,nmel);
D = 20*log10(max(1e-10,mlmx(:,1:(swin/2+1))*abs(D)));

% Only look at the top 60 dB
D = max(D, max(D(:))-60);

% The raw onset decision waveform - half-wave rectified first-order difference
mm = (mean(max(0,diff(D')')));

% dc-removed onset waveform
onset = filter([1 -1], [1 -.99],mm);

