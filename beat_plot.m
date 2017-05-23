function beat_plot(b,s,tt,ff,D,yy)
% beat_plot(b,s,tt,ff,D,yy)
%   Plot vertical lines corresponding to the beat times in b on the
%   current plot.  Optional s is the plot style.
%   Optional tt,ff,D define underlying image for imagesc(tt,ff,D)
%   Optional yy is used as the y index for plotting the points (one
%   val) or line (two vals).
% 2012-03-27 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; s = ''; end
if nargin < 3; tt = []; end
if nargin < 4; ff = []; end
if nargin < 5; D = []; end
if nargin < 6; yy = []; end

if iscell(b)
  
  % passed a cell array of beats - plot each one
  nb = length(b);
  % If yy is passed in, it's the maximum value to spread the beat
  % tracks over
  ymax = max([1;yy]);
  
  for i = 1:nb
    % spread out the points between 0 and 1
    yy = ymax * ((i-0.5)/length(b))*[1;1];
    % recurse
    beat_plot(b{i}, s, tt, ff, D, yy);
  end  

else

  if length(s) == 0; s = '-r'; end

  if length(D) > 0
    imagesc(tt,ff,D);
    axis('xy');
    colormap(1-gray);
  end

  ax = axis;

  if length(yy) == 0; yy = [ax(3) ax(4)]'; end
  if length(yy) == 2; yy = [yy(1); yy(2)]; end % ensure a column
  
  hold on;
  plot([b;b], repmat(yy,1,length(b)), s);
  hold off;

end
