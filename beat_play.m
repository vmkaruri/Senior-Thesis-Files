function y = beat_play(b,d,sr)
% y = beat_play(b,d,sr)
%   Synthesize clicks for a beat track.
%   <b> is a list of beat times (in sec)
%   <d> is a maximum duration (in sec) or a waveform to add in.
%   <sr> is the sampling rate for the output (and for waveform <d>).
%   <y> is the synthesized waveform.
% 2012-03-27 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; d = 0; end
if nargin < 3; sr = 11025; end

% Maybe read in a soundfile for d?
if ischar(d); [d,sr] = audioread(d,sr,1); end

if length(d) == 1
  l = round(d*sr);
  d = zeros(l,1);
end

l = length(d);
% The waveform to use for each beat
tdur = 0.1;
tt = (0:round(tdur*sr))';
% 100ms pip @ 2khz
fblip = 2000;
blip = tt.*exp(-tt/((tdur*sr)/10)).*cos(2*pi*tt/sr*fblip)/200;
%% 100ms exponential decay white noise burst
%blip = exp(-tt/((tdur*sr)/20)).*randn(length(tt),1)/5;

lblip = length(blip);

bsamp = round(b*sr);

% remove beats that would run off end
if l > 0
  bsamp = bsamp(bsamp < (l-lblip));
else
  l = max(bsamp)+lblip;
end

% ensure d is long enough
d(l) = 0;

for bbx = 1:length(bsamp)
  bb = bsamp(bbx);
  d(bb+[1:lblip]) = d(bb+[1:lblip]) + blip;
end

if nargout == 0
  soundsc(d,sr);
else
  y = d;
end
