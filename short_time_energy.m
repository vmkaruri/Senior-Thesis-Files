% Computation of ST-ZCR and STE of a speech signal.
%
% Functions required: zerocross, sgn, winconv.
%
% Author: Nabin Sharma
% Date: 2009/03/15
function [E] = short_time_energy(audio_file)
    [x,Fs] = audioread(audio_file); % word is: so
    x = x.';

    N = length(x); % signal length
    n = 0:N-1;
    ts = n*(1/Fs); % time for signal

    % define the window
    wintype = 'rectwin';
    winlen = 201;
    winamp = [0.5,1]*(1/winlen);

    
    % find the short time energy 
    E = energy(x,wintype,winamp(2),winlen);

    % time index for the ST-ZCR and STE after delay compensation
    out = (winlen-1)/2:(N+winlen-1)-(winlen-1)/2;
    t = (out-(winlen-1)/2)*(1/Fs);
end