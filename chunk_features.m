function [means, stds, others] = chunk_features(features, states_voiced, states_speaking, minutes_per_chunk)
% features = chunk_features(features, states_voiced, states_speaking, minutes_per_chunk)
%
% Returns various features over minutes_per_chunk minutes, default 5.
%
% If minutes_per_chunk is inf, then entire file is one chunk.
%
% means(f, t, s) is the mean of feature f over chunk t for speaker s.
% stds(f, t, s) is the standard deviation of feature f over chunk t for speaker s.
% others(f, t, s) is the value of feature f over chunk t for speaker s.
%
%  Eight Features (for means/stds):
%
%  1 - formant frequency (Hz)
%  2 - confidence in formant frequency
%  3 - spectral entropy
% 
%  4 - value of largest autocorrelation peak
%  5 - location of largest autocorrelation peak
%  6 - number of autocorelation peaks
% 
%  7 - energy in frame
%  8 - time derivative of energy in frame
%
% Other Features:
%  1 - Average length of voiced segment (seconds)
%  2 - Average length of speaking segment (seconds)
%  3 - Fraction of time speaking
%  4 - Voicing rate: number of voiced regions per second speaking
%  5 - Fraction speaking over: fraction of time that you & any other
%        speaker are speaking.
%  6 - Average number of short speaking segments (< 1 sec) per minute.
%        Only segments within 1 sec of similarly short segments of any
%        other speaker are included.

num_frames = size(states_speaking, 1);
num_speakers = size(states_speaking, 2);


if (minutes_per_chunk ==inf)
    frames_per_chunk = num_frames -1;
    minutes_per_chunk = round(num_frames/(8000/128)/60);
else
    frames_per_chunk = round(8000/128*60*minutes_per_chunk);
end
num_chunks = floor(num_frames / frames_per_chunk);

means = zeros(8, num_chunks, num_speakers);
stds = zeros(8, num_chunks, num_speakers);
others = zeros(6, num_chunks, num_speakers);

for s = 1:num_speakers
    regions_speaking{s} = states_to_regions(states_speaking(:, s));
end

for t = 1:num_chunks
    frames = frames_per_chunk*(t-1)+1:frames_per_chunk*t;
    
    for s = 1:num_speakers
        voiced = states_voiced(frames, s);
        speaking = states_speaking(frames, s);
        
        for feature = 1:8
            chunk = features(feature, frames, s);
            
            if feature == 5
                % "chunk" can be zero if there is no autocorrelation peak after
                % the first.  This shouldn't happen very often.
                chunk(chunk > 0) = 1 ./ chunk(chunk > 0);
            end
            
            % means(feature, :) = mean(x);
            chunk(voiced == 1) = [];
            if length(chunk) < 10
                warning(['length(chunk) < 10 chunk:' num2str(t) ' speaker: ' num2str(s) ' ' num2str(feature)]);
                chunk(voiced == 1) = 0;
            end
            if any(~isfinite(chunk))
                error('Feature isn''t finite');
            end
            means(feature, t, s) = mean(chunk);
            
            if feature == 1 || feature == 5 || feature == 7
                if (length(chunk(chunk > 0)) == 0)
                    stds(feature, t, s) = 0;
                else
                    stds(feature, t, s) = std(chunk) / mean(chunk(chunk > 0));
                end
            else
                stds(feature, t, s) = std(chunk);
            end
        end
        
        if (sum(voiced-1) >  1)
            % Average length of voicing segments (seconds)
            regions_voiced = states_to_regions(voiced);
            others(1, t, s) = mean(regions_voiced(2, :) - regions_voiced(1, :)) / (8000/128);

            % Average length (sec) and fraction of time speaking
            rs = regions_speaking{s};
            rs = rs(:, find(mean(rs) >= frames(1) & mean(rs) < frames(1) + frames_per_chunk));
            others(2, t, s) = mean(rs(2, :) - rs(1, :)) / (8000/128);
            others(3, t, s) = mean(speaking) - 1;

            % Voicing rate: number of voiced regions per second speaking
            others(4, t, s) = size(regions_voiced, 2) / length(find(speaking == 2));

            % Percent speaking over
            % speaking over = Any two speaking at once.
            other_speakers = [1:(s-1) (s+1):num_speakers];
            other_speaking = states_speaking(frames, [1:(s-1) (s+1):num_speakers]);
            others(5, t, s) = length(find(any(other_speaking == 2, 2) & speaking == 2)) / length(speaking);

            % Num short (< 1 sec) segments per minute, that are close (within
            % in second) to other speaker's short segments.

            if (sum(other_speaking-1) >  1)

                ind = find(rs(2, :) - rs(1, :) < 8000/128 & rs(2, :) - rs(1, :) > 8000/128/5);  % Find short segments
                delta_frames = 1 * 8000/128;  % Other segments must be within this much time of current segment.
                keepers = [];
                for j = ind
                    for other_s = other_speakers
                        ind_other = close_indices(regions_speaking{other_s}, rs(:, j), delta_frames);
                        if ~isempty(ind_other)
                            keepers = [keepers rs(:, j)];
                            break;
                        end
                    end
                end

                others(6, t, s) = size(keepers, 2) / minutes_per_chunk;
            else
                warning('other never speaks')
            end
        else
            warning('no voiced segment in chunk');
        end

    end
end

function indices = close_indices(other, this, delta_frames)

indices = find(other(2, :) - other(1, :) <= 62 & other(2, :) > this(1) - delta_frames & other(1, :) < this(2) + delta_frames);
