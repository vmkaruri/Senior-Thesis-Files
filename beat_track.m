function [b,sgram,tt,ff] = beat_track(d,sr, VERBOSE)
% [b,sgram,tt,ff] = beat_track(d,sr, VERBOSE)
%    Wrapper to run beat tracking on a soundfile.
%    <d> is a waveform at <sr>, or is the name of a soundfile.
%    <b> returns the beat times.
% 2012-03-27 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3; VERBOSE = 0; end

if nargout == 0; VERBOSE = 1; end

% Maybe read in the waveform
if ischar(d); 
  targetsr = 11025;
  forcemono = 1;
  [d,sr] = audioread(d,targetsr,forcemono);
end

% Get the onset envelope
[onset, osr, sgram] =  beat_onset(d,sr);
% axes for sgram, in case we need them
tt = [0:size(sgram,2)-1]/osr;
ff = 1:size(sgram,1);

% Estimate the global tempo (and plot global a/c)
if VERBOSE
  subplot(211)
end
tempo = beat_tempo(onset, osr, VERBOSE); % plots autocorrelation if VERBOSE

% Do the beat time search
beats = beat_simple(onset, osr, tempo);

if VERBOSE
  disp(['Global tempo = ', num2str(tempo),' BPM']);

  % Plot the results
  subplot(212)
  beat_plot(beats, '', tt, ff, sgram);
end

if nargout < 1
  % Update the graphics
  drawnow();
  % Play the results
  %beat_play(beats, d, sr);
else
  b = beats;
end
