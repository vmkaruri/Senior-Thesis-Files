function [t,xcr] = beat_tempo(onsetenv,osr,disp)
% [t,xcr] = beat_tempo(onsetenv, osr, disp)
%   Estimate the tempo (in BPM) of the onset envelope onsetenv (at
%   frame rate osr).
%   xcr returns the cross-correlation from which the tempo peak was
%   picked.  We plot it if disp == 1.
% 2012-03-27 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3; disp = 0; end

% autoco out to 4 s
acmax = round(4*osr);

% Find rough global period
% Only use the 2nd 60 sec to estimate global pd (avoid glitches?)

maxd = 60;
maxt = 120; % sec
maxcol = min(round(maxt*osr),length(onsetenv));
mincol = max(1,maxcol-round(maxd*osr));

xcr = xcorr(onsetenv(mincol:maxcol),onsetenv(mincol:maxcol),acmax);

% find local max in the global ac
rawxcr = xcr(acmax+1+[0:acmax]);

% window it around default bpm
% BPM prior is centered on bpmmean, with a log_2-domain Gaussian
% halfwidth of bpmsd
bpmmean = 120;
bpmsd = 0.9;

% What BPM does each bin of the autocorrelation correspond to?
bpms = 60*osr./([0:acmax]+0.1);
% Calculate the log-normal windowing
xcrwin = exp(-.5*((log(bpms/bpmmean)/log(2)/bpmsd).^2));

% Apply the weighting
xcr = rawxcr.*xcrwin;

xpkix = min(find(xcr == max(xcr)));  

% Convert period to BPM
t = 60/(xpkix / osr);

if disp
  plot([1:length(xcr)]/osr,xcr);
  ax = axis;
  hold on;
  plot([1 1]*xpkix/osr, [ax(3) ax(4)], '-r');
  hold off;
  title(['Onset function autocorrelation, tempo = ', ...
         sprintf('%.1f',t),' BPM']);
  xlabel('lag / sec');
end
