function  [Cc,d0] = dime(Z,emb,fileout)

% DIME(Z,EMB,FILEOUT)
%       Estimates correlation dimension of a scalar or vector
%     timeseries. For vectors  Z (N by n) where n is dimension of
%     each vector, N is a number of vectors (that is Z is a set of
%     N points in n-dimensional space.
%     This program first computes a "correlation sum" C  as a
%     function of distance  D.  C is equal to a number of pairs of
%     points of vector Z separated by a distance less than D in 
%     a phase space. Correlation dimension is equal to  d_corr  if
%     in a certain range of distances D the correlation sum behaves
%     like 
%               d_corr                        d(log(C))
%          C ~ D      ,    that is  d_corr =  --------
%                                             d(log(D)).
%     This program estimates this derivative (the slope of a curve
%     C(D) in log-log scale) for different ranges of distance D by
%     fitting a straight line over several points, corresponding to 
%     different distance bins.
%
%     If  Z  is a matrix (N by n), N > n  it is treated as a set of
%     N vectors in  n-dimensional phase space and distances  D  are
%     Eucledian distances in this phase space. In this case d_corr
%     approximates the dimension of a manifold containing the set Z.
%     If  Z  is a vector (N by 1) the program estimates its dimension
%     in the "embedding" k-dimensional phase space formed by each k
%     successive components of vector Z:
%     [Z(1) Z(2) ... Z(k);  Z(2) Z(3) ... Z(k+1);
%      Z(3) Z(4) ... Z(k+2); ... ;  Z(N-k+1) Z(N-k+2) ... Z(N)]
%     for different values of  k.
%
%   DIME(Z) performs correlation dimension computations for matrix
%     or vector Z  (for vectors the "embedding" procedure is applied).
%   DIME(FILEIN) takes matrix Z from a MAT-file  FILEIN.
%     If the needed dataset is stored in FILEIN under a name different
%     from Z it can be retrieved if the string  FILEIN is written as
%     'FILE(VAR)' or 'FILE VAR',  where  FILE is the name of a MAT-file
%     and  VAR  is the name of the variable substituting  Z.
%   DIME(Z,EMB) in case when Z is vector computes correlation dimension
%     in "embedding" phase space specified by  variable EMB: when EMB
%     is a vector (e.g. [3 4 6]) the computations are performed for 
%     "embedding" spaces of dimensions equal to each component of vector
%     EMB successively;  when EMB is a number (EMB>1) the program spans
%     all embedding dimensions from 2 to EMB. When EMB is not spesified
%     the default value (max. "embedding" dimension) is used.
%     The dependence of correlation dimension on the "embedding"
%     dimension  d_corr(d_emb)  is plotted. If d_corr tends to a limit
%     as  d_emb  goes to infinity this limiting value is assumed to be
%     the "true" dimension of the timeseries Z.
%   DIME(Z,FILEOUT) also saves main output variables ti a MAT-file
%     FILEOUT. The list of saved variables is specified by the string
%     "savevar" in the text of the program.

%  Connenction:  FITLINE.M  - routine for fitting a straight line
%   through the datapoints and estimating errors of the fit.

%  Kirill Pankratov, Feb 20, 1994.

% Computing parameters and defaults  ********************************
isnewfig = 1;   % New figure are created for plotting
embdflt = 8;    % Default for max embedding dimension
nd = 32;        % Number of distance bins
nd0 = nd/2-8;   % First of distance bins
naver = 6;      % Number of points for local corr. dimension
                % estimation (averaging)
scmult = 1.5;   % Multiplicator for scale (relative to variance of 
                %   the dataset)
af = .5;        % Smoothing coef. for range estimates
nfilt = 6;      % Number of smoothing iterations
lchunk = 25;    % Length of a chunk of a set for computing distances


 % Output variables if saved
savevar = 'Cc d0 p perr dest n0 pr prerr dm';

% Plotting parameters  **********************************************
figsz = [550 770];           % Size of figure (in pixels)

 % Axes positions (normalized):
posax1 = [.12 .55 .80 .38];  % First subplot (for the case of 2)
posax2 = [.12 .10 .80 .38];  % Second subplot (for the case of 2)
posax3 = [.21 .08 .62 .20];  % Third subplot

 % Axes titles and labels:
 % First subplot - C(d) in log-log scale
ttl1txt = 'Correlation  sum';
xlab1txt = 'Normalized  distance';
ylab1txt = 'C_corr';
 % Second subplot - local d_corr estimates with errorbars
ttl2txt = 'Local corr. dimension estimates for different ranges';
xlab2txt = 'Normalized  distance';
ylab2txt = 'dimension';
 % Third subplot - d_corr(d_emb)
ttl3txt = 'Correlation vs. embedding dimension';
xlab3txt = 'Embedding dimension';
ylab3txt = 'Corr. dimension';

 % Linestyles and colors
lnstcs = ['.'; 'o']; % Linestyle for C(d) (first subplot)
lnsta = '-';         % Linestyle for asymptotics lines
mrksz = 4;           % Marker size for first subplot


% Now begin input handling and computations ***********************
 % Handle input ..................................................
isfilein=0; isfileout=0; isembin=0; isemb=0; isrename=0;
if nargin==0   % Not enough input arguments
  disp([10 '  Error: there must be at least one input argument Z' ...
  ' (1-d timeseries or n-d dataset).' 10])
   return
end
if nargin >= 1  % Can be Z or input file name
  if isstr(Z)
    issep = (Z=='(')|(Z==')')|(Z==',')|(Z==' ');
    if any(issep), fnd = find(issep);
      s1 = Z(1:fnd(1)-1);
      s2 = fnd(1)+1:length(Z); s2 = s2(issep(s2)==0);
      s2 = Z(s2);
      if s2~=[], isrename=1; end
    else, s1 = Z;
    end
    if (exist(s1)==2)|(exist([s1 '.mat'])==2) isfilein=1;
    else
      disp([10 '  Error: file ' Z ' does not exist' 10])
      return
    end
  else
    if max(size(Z))<=1
    disp([10 ' Error: size of the dataset Z must be larger than 1' 10])
    return
    end
  end
end
if nargin >= 2  % Can be EMB or output file name
  if isstr(emb)
    fileout = emb; isfileout=1; emb=embdflt;
  else
    isembin=1;
    if any(emb<=0)
      disp([10 ' Embedding dimension must be positive!' 10])
      emb=embdflt
    end
  end
end
if nargin >=3
  if ~isstr(fileout)
    disp([10 '  Error: output filename must be a string.' ...
    10 '  Output is not saved.' 10])
  else, isfileout=1;
  end
end

if isfilein, eval(['load ' s1]);
  if isrename, eval(['Z=' s2 ';']); end
end
if ~isembin, emb = embdflt; end

 % Determine size of Z and embedding procedure ....................
sz0 = size(Z);                  % Initial size of a dataset Z
if sz0(2)>sz0(1), Z=Z'; end     % Make each point a row vector
if any(sz0==1),  isemb=1;  end  % 1-D timeseries, needs embedding
emb = emb(:);
if (size(emb,1)==1)&(emb>1),  emb = 2:emb;  end
sz0 = size(Z);         % Correct if necessary
if ~isemb, emb=1; end
lemb = length(emb);    % Number of the embedding dimensions
if sz0(1)<=max(emb)
  disp([10 '  Error: dataset Z is too short, can not compute.' 10])
  return
end

%size(Z), emb, [isfilein isemb isfileout isrename]
%break

 % Check number of subplots and adjust axes positions . . . .
if isemb     % If embedding is applied (1-D timeseries)
  np = 3;    % 3 subplots total
   % Now make room for the 3-rd axes: 
  part12 = 1-posax3(2)-posax3(4);
  posax1(2) = (1-part12)+posax1(2)*part12;
  posax1(4) = posax1(4)*part12;
  posax2(2) = (1-part12)+posax2(2)*part12;
  posax2(4) = posax2(4)*part12;
else         % No embedding, compute in actual phasespace
  np = 2;
end  %   np - numb. of subplots


 % Create necessary vectors and numbers:
d0 = 10.^((nd0+1:nd+nd0)/(nd+nd0))/10;  % Distance bins
C = zeros(size(d0));                    % Vector for corr. sum
nchunk = ceil(sz0(1)/lchunk);           % Numb. of "chunks"


% Begin computation of corr. dimension ******************************

for jemb = 1:lemb  % Begin embedding dimension cycle ````````````````0

   % Form matrix Z in the "embedding space" if needed . . . . . .
  if isemb
    cemb = emb(jemb);         % Current embedding dimension
    vz = (1:sz0(1)+1-cemb)';  % Auxillary (vert. size)
    eemb = ones(1,cemb);      % Auxillary (unit vector, hor. size)
    vemb = 0:cemb-1;          % Auxillary (hor. size)
    Z = Z(:,eemb);            % Make a matrix for embedding
    Z = Z(vz(:,eemb)+vemb(ones(size(vz)),:));  % Embedding itself
  end

   % Estimate the scales . . . . . . . . . . . . .
  sz = size(Z);    % Current size of Z
  zmean = mean(Z); % Mean of each component
  Z = Z-zmean(ones(sz(1),1),:);
  Z2 = sum(Z'.^2);    % Squared length of all vectors
  zvar = sum(Z.^2);   % Variances for each component

   % Take care of the distance bins vector
  dscale = sqrt(sum(zvar)/sz(1));
  dc = scmult*dscale*d0;  % Actual distance bins vector
  dc2 = dc.^2;  % Squared distances


  for jch1 = 1:nchunk % Begin 1-st chunk cycle ` ` ` ` ` ` ` ` ` ` 1
    num1 = (jch1-1)*lchunk+1:min(sz(1),jch1*lchunk);

    for jch2 = 1:jch1 % Begin 2-nd chunk cycle ` ` ` ` ` ` ` ` ` 2
      num2 = (jch2-1)*lchunk+1:min(sz(1),jch2*lchunk);

      % Create a distance matrix
      D2 = ones(size(num1'))*Z2(num2)+Z2(num1)'*ones(size(num2));
      D2 = D2-2*Z(num1,:)*Z(num2,:)';
           % Distance matrix for two chunks

      % Add to correlation sum
      D2 = D2(:);
      if jch1==jch2, cf=1; else, cf=2; end
      C = C+cf*sum(D2(:,ones(size(dc2)))<=dc2(ones(size(D2)),:));
                  % This is similar to a histogram
    end  %  End jch2 cycle ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 2

  end  %  End jch1 cycle ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 1

  C = C-sz(1);
  C = C+(C==0);             % Regularization (to exclude log(0))
  C = C/sz(1)/(sz(1)-1);   % Normalization - 1/N/(N-1)
  Cc(:,jemb) = C';
  if isemb, disp([' Done embedding dimension ' num2str(cemb)]), end
end   % End embedding cycle '''''''''''''''''''''''''''''''''''''''''0
%  End of main computation routine - correlation sum for distances d
%  Output here is Cc (normalized correlation sum) ****************EEE


 % Log. for plotting and fitting
Clog = log10(Cc);
d0log = log10(d0);


% Processing and fitting the results ********************************

 % Fitting straight lines to correlation sums as functions of d 
 % Corr. dimension estimates are   d(log(C))/d(log(d)), that is
 % the slope of the line fitted over numb. of  points  NAVER.
numest = 3:nd-naver-1;  % ## of points for which fit is made
lest = length(numest);  % Number of fits
p = zeros(lest,2*lemb);
perr = p;
dest = zeros(lest,1);
for jest = 1:lest  % Fitting cycle
  num = numest(jest):numest(jest)+naver-1;
  dest(jest) = mean(d0log(num));
  for je = 1:lemb
    dm = (je-1)*2+[1 2];
    [p(jest,dm),perr(jest,dm)] = fitline(d0log(num),Clog(num,je));
  end
end
dest = 10.^dest;
      % parameters of the curves:  y = p(1)*x+p(2)
      % perr - errors of the estimates of p(1) and p(2)

 % Find the "scaling ranges", where dimension estimates are stable
 % over a certain span od scales, look for minima of curvature
 % of function  log(C) vs. log(d) 
num = 2:lest-1;
curv = Clog(num-1,:)+Clog(num+1,:)-2*Clog(num,:);
curv = (curv.^2)./(1+(Clog(num-1,:)-Clog(num+1,:)).^2/4);
 % Simple filtering (smoothing) of curvature estimates
num = 2:lest-3;
for jiter = 1:nfilt
 curv(num,:) = (1-af)*curv(num,:)+af/2*(curv(num-1,:)+curv(num+1,:));
end
 % Find minima of smoothed curvatures:
cond = (curv(num,:)<curv(num-1,:))&(curv(num,:)<curv(num+1,:));
[fnd,dm] = find(cond);  % dm - ## of embedding dimensions
fnd = fnd+2;        % ## of minima points
lfnd = length(fnd); % Numb. of minima


% Plotting the results **********************************************
if isnewfig, figure, end   % Create a new figure if needed

 % Plot correlation sums vs. distance in log-log scale ............
axvec = [min(d0) 1 min(min(Cc)) 1];  % Axes limits
lns = ['. '; '. '];                     % "Default" linestyle
lns(1:size(lnstcs,1),1:size(lnstcs,2)) = lnstcs;

 % Begin 1-st subplot . . . . .
subplot(np,1,1)
posf = get(gcf,'pos');
posf(3:4) = figsz;
set(gcf,'pos',posf)
 % Adjust color order
colo = get(gca,'colororder');
colo = reshape([colo'; colo'],3,2*size(colo,1))';
set(gca,'colororder',colo)

loglog(d0,Cc,lns(1,:))  % Plot lines
hold on
ln1=loglog(d0,Cc,lns(2,:));  % Plot markers
set(ln1,'markersize',mrksz)
axis(axvec)
set(gca,'units','norm','pos',posax1)
title(ttl1txt)
xlab1 = text('string',xlab1txt);
set(xlab1,'pos',[.4 axvec(3)],'vertical','base')
ylab1 = get(gca,'ylabel');
set(ylab1,'string',ylab1txt)

 % Calculate and plot the d_corr estimate over different ranges .......
 % Compute dimension estimates over "scaling regions"
 % and put asymptitics lines in log-log plot
for jr = 1:lfnd  % For all scaling region  
  n0 = numest(fnd(jr));
  num = n0-naver:n0+naver;
  num = num(num>0&num<nd);
   % The fit itself:
  [pr(jr,:),prerr(jr,:)] = fitline(d0log(num),Clog(num,dm(jr)));
  nump = num(1)-2:max(num)+2;   % Plotting ##
  nump = nump(nump>0&nump<nd);
  ya = 10.^(pr(jr,1)*d0log(nump)+pr(jr,2));  % Asymptotics
  plot(d0(nump),ya,lnsta)   % Plot asymptotics . . . . . . .
   % Put numbers (dimension estimates) near asympt. lines:
  numt = min(jr+3,length(nump));
  ta(jr) = text(d0(nump(numt)),ya(numt-1),sprintf('%5.2f',pr(jr,1)));
  set(ta(jr),'fontsize',10)
end

 % Find maximum dimension estimates for each embedding space
pa = zeros(lfnd,max(dm));  % Auxillary matrix
pa((1:lfnd)'+(dm-1)*lfnd) = pr(:,1);
pa = max(pa);
num = (1:size(pr,1))'*ones(1,size(pa,2)); % ## of elements
[fndm,num] = find(pr(num)==pa(ones(lfnd,1),:));
%[fndm,num] = find(pr(:,ones(1,max(dm)))==pa(ones(lfnd,1),:))
prm = pr(fndm,:);       % max. slope
prmerr = prerr(fndm,:); % err. of max slope
numd = dm(fndm);


 % Plot dimension estimates with errorbars ............................
subplot(np,1,2)
errorbar(dest(:,ones(1,lemb)),p(:,1:2:lemb*2),perr(:,2:2:lemb*2))
ylim = floor(min(min(p(:,1:2:lemb*2))));
ylim = [ylim ceil(max(max(p(:,1:2:lemb*2))))];
set(gca,'Xlim',axvec(1:2),'Ylim',ylim,'Xscale','log','pos',posax2)
grid
 % Labeling . . . . . . . . .
for jemb = 1:lemb  % Labeling the dimension estimates curves
  num = min(1+3*jemb,size(p,1));
  te = text(dest(num),p(num,1+2*(jemb-1))+.1,num2str(emb(jemb)));
end
ttl2 = get(gca,'title');  % Title of the 2-nd plot
set(ttl2,'units','norm','pos',[.5 1.0],'string',ttl2txt)
xlab2 = text('string',xlab2txt);
set(xlab2,'pos',[.4 ylim(1)],'vertical','base')
ylab2 = get(gca,'ylabel');
set(ylab2,'string',ylab2txt)


 % Plot correlation vs. embedding dimension .........................
if isemb
  subplot(np,1,3)
  set(gca,'units','norm','pos',posax3,'box','on')
  hold on
  errorbar(emb(numd),prm(:,1),prmerr(:,2)) % Plot itself
  plot(0:max(emb),0:max(emb),lnsta)   % For comparison:
                                      %  d_corr=d_emb  line
  axis([0 max(emb) 0 ceil(max(prm(:,1)))]) % Unequal axes limits
  grid
   % Labeling . . . . . . . . .
  ttl3 = get(gca,'title');
  set(ttl3,'string',ttl3txt)
  xlab3 = get(gca,'xlabel');
  set(xlab3,'string',xlab3txt)
  ylab3 = get(gca,'ylabel');
  set(ylab3,'string',ylab3txt)
  txt31 = text(emb(1)+.3,prm(1,1),'d_corr(d_emb)');
  txt32 = text(.6,.4,'d_corr=d_emb');
  hold off
end

orient tall  % For printing

% Save output
if isfileout, eval(['save ' fileout ' ' savevar]), end




