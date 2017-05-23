% MatLab Code:
% Take the average and peak envelopes of an input signal, and plot all three
% together.
% Expects a .wav soundfile as input.
% usage: ampenv('myfile.wav')

% This beginning part just defines the function to be used in MatLab: takes a 
% "wav" sound file as an input, and spits out a graph.

function [maxenv] = amplitude_descriptor( file )
    % Reads in the sound file, into a big array called y.

    y = audioread( file );

    % Normalize y; that is, scale all values to its maximum. Note how simple it is 
    % to do this in MatLab.

    y = y/max(abs(y));

    % If you want to hear the sound after you read it in, uncomment this next line. 
    % Sound just plays a file to the speaker.

    %sound(y, 44100);

    % Go through the input signal, taking the maximum value of the previous 
    % peakwindowsize number of samples. We do this in the same way as we did the 
    % average, but now we use max instead of sum.

    for k = peakwindowsize:length(maxenv)							
     maxenv(k) = max(y(k-peakwindowsize:k));							
    end
end