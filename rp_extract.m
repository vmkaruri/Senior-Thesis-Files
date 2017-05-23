function varargout = rp_extract(varargin)

% RP_EXTRACT - Rhythm Patterns Audio Feature Extraction
% version: 0.6411
% (pls. edit also global version variable below help text)
%
%                             (c) 2005 - 2008 by Thomas Lidy
%   Institute of Software Technology and Interactive Systems
%                            Vienna University of Technology
%             http://www.ifs.tuwien.ac.at/mir/downloads.html
%
% based on:
%  Musik Analysis Toolbox for Matlab Version 0.4
%  (c) 2002 Elias Pampalk
%  Austrian Research Institute for Artificial Intelligence
%  http://www.oefai.at/~elias/
%
% and further development on ma_gui_v04 until version 0.48: 
%  2003 - 2005 by          Thomas Lidy, Stefan Leitich
%
%
% RP_EXTRACT: 
%   Music Analysis, Audio Indexing, Similarity
%   Extraction of the following features:
%   - Rhythm Patterns
%   - Statistical Spectrum Descriptor
%   - Rhythm Histogram
%   - Approximate bpm measure
%
% VERSION HISTORY:
%
% 0.51  18/05/05 - first stable rp_extract (tml)
% 0.51b 06/07/05 - corrected a bug, appearing in matlab < 7 only (') (tml)
%                  improved options path processing issues (tml)
%                  process now input filelists with multiple TAB-sep. columns, take only 1st column (tml)
% 0.52  30/08/05 - option to include DC componont in modulation spectrum
%                  option to extract phase as separate feature vector
% 0.53  21/11/05 - changed space processing in read_options_from_file
%                  read_options_from_file now initializes with default options (to set options that are not given in file)
% 0.54  13/02/06 - on the fly decoding of .ogg files through oggdec
% 0.55  14/02/06 - optional resampling of all files through sox
% 0.551 14/02/06 - improved error handling with too short audio files
% 0.56  09/03/06 - writing $DATA_TYPE header line with content type information
% 0.57  21/04/06 - solved bug in phon calculation when n_bark_bands < 24
% 0.58  17/07/06 - added flac support (rn)
% 0.59	18/07/06 - added support for Monkey's Audio Codec (.ape), mac binary is required (rn)
% 0.60	21/07/06 - modified BPM output in that it produces valid SOMLib output files (rn)
% 0.61  08/02/07 - two minor bugs solved (tml)
% 0.62  29/06/07 - enabled extraction of SSD and RH features only (without RP) (tml)
% 0.621 22/08/07 - corrected bug in 3 argument processing (tml)
% 0.63  14/08/08 - introduced temporal statistics over SSD features: TSSD (tml)
%                  changed header: added $DATA_DIM (instead of including it in $DATA_TYPE)
% 0.64  15/08/08 - introduced temporal statistics over RH features: TRH (tml)
%                  introduced Modulation Frequency Variance Descriptor: MVD
%                  added version info in vector header: $EXTRACTOR
%
% USAGE:
%
% 1) Batch Audio Feature Extraction
%    for multiple files, results written to vector files
%
%    RP_EXTRACT(filelist, dirs [, opts])
%
%    RP_EXTRACT(filelist, optionsfile)
%
%    RP_EXTRACT(optionsfile)
%
%    filelist:  list of audio files, either:
%               filelistarray | filelistfile
%
%    dirs:      directory options, either:
%               dirs_struct | dirs_optionsfile
%
%    opts:      processing options, either:
%               opts_struct | opts_optionsfile
%               (if omitted, default options will be used)
%
%
%   filelistarray: Matlab cell array containing all files to be processed
%   filelistfile:  Ascii file listing all files to be processed
%   
%   dirs_struct:   Matlab struct, containing input/output directory etc.
%   dirs_optionsfile: Ascii file containing directory options
%
%   opts_struct:   Matlab struct with feature extraction options
%   opts_optionsfile: Ascii file containing feature extraction options
%
%   optionsfile:   Ascii file containing both directory and processing options
%
%   If the only argument is optionsfile, it has to include reference to 
%   a filelist Ascii file as well.
% 
%      For more information on directory and processing options
%      see rp_options.txt example file
%
%
% 2) Audio Feature Extraction from single audio file
%    results are returnd as Matlab vectors
%
%    [feat] = RP_EXTRACT('func','process_audiofile',input_filename [, opts])
%
%    input_filename:  complete path to audio file
%
%    opts:     Matlab struct containing feature extraction options 
%              (if omitted, default options will be used)
%
%   return value [feat] is a Matlab struct containing the following fields:
%    feat.rp:  vector with Rhythm Patterns features
%    feat.ssd: vector with Statistical Spectrum Descriptor
%    feat.rh:  vector with Rhythm Histogram features
%    feat.bpm: single value containing approximate bpm peak value

global extr_version;
extr_version = 'Matlab rp_extract v 0.6411 by tml';

% ARGUMENT PROCESSING
% -------------------

if (nargin == 0)
	help rp_extract;
else


% INITIALIZATION OF GLOBAL CONSTANTS
% ----------------------------------

global CONST_bark;
global CONST_spread;
global CONST_phon;
global CONST_loudn_freq;
global CONST_eq_loudness;
global CONST_loudn_bark;
%global CONST_filt;

% CRITICAL BARK BANDS
CONST_bark = [100 200 300 400 510 630 770 920 1080 1270 1480 1720 2000 2320 2700 3150 3700 4400 5300 6400 7700 9500 12000 15500];

% SPREADING FUNCTION FOR SPECTRAL MASKING
% CONST_spread contains matrix of spectral frequency masking factors
n_bark_bands = 24;
for i = 1:n_bark_bands,
	CONST_spread(i,:) = 10.^((15.81+7.5*((i-(1:n_bark_bands))+0.474)-17.5*(1+((i-(1:n_bark_bands))+0.474).^2).^0.5)/10);
end


% PHON CALCULATION
CONST_phon = [3 20 40 60 80 100 101];

% LOUDNESS LEVELS
% 22 frequencies corresponding to following loudness sensation values 
CONST_loudn_freq = [31.62, 50, 70.7, 100, 141.4, 200, 316.2, 500, 707.1, 1000, 1414, 1682, 2000, 2515, 3162, 3976, 5000, 7071, 10000, 11890, 14140, 15500]; 

CONST_eq_loudness = [
    55   40   32   24   19   14   10    6    4    3    2    2    0   -2   -5   -4    0    5   10   14   25   35; 
    66   52   43   37   32   27   23   21   20   20   20   20   19   16   13   13   18   22   25   30   40   50;
    76   64   57   51   47   43   41   41   40   40   40   39.5 38   35   33   33   35   41   46   50   60   70;
    89   79   74   70   66   63   61   60   60   60   60   59   56   53   52   53   56   61   65   70   80   90;
    103   96   92   88   85   83   81   80   80   80   80   79   76   72   70   70   75   79   83   87   95  105;
    118  110  107  105  103  102  101  100  100  100  100   99   97   94   90   90   95  100  103  105  108  115];

% We have the loudness values for the frequencies in CONST_loudn_freq
% now we calculate in CONST_loudn_bark a matrix of loudness sensation values for the bark bands margins

i=1; j=1;
for bsi = CONST_bark,
    while j<length(CONST_loudn_freq) & bsi>CONST_loudn_freq(j),
        j=j+1;  
    end
    j=j-1;
    if ~isempty(find(CONST_loudn_freq==bsi)),  % loudness value for this frequency already exists
        CONST_loudn_bark(:,i) = CONST_eq_loudness(:,find(CONST_loudn_freq==bsi));
    else
        w1 = 1 / abs(CONST_loudn_freq(j) - bsi);
        w2 = 1 / abs(CONST_loudn_freq(j + 1) - bsi);
        CONST_loudn_bark(:,i) = (CONST_eq_loudness(:,j).*w1 + CONST_eq_loudness(:,j+1).*w2) / (w1 + w2);
    end
    i=i+1;
end 




% CONTINUATION OF ARGUMENT PROCESSING
% -----------------------------------

   if (strcmp(varargin(1),'func')) % ENABLES CALL OF SUB-FUNCTIONS
   	
	varargin = varargin(2:end);
	
	if ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
	
	try
		if (nargout)
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
		feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end
	
	end    
   else	
   	% PROCESS ARGUMENTS FOR FEATURE EXTRACTION	
   	
	if (nargin > 3)
		message('ERROR: inapropriate number of arguments. see ''help rp_extract''');
		return
	end
	
	if (nargin==1)	% all options, including filelistfile, are read from optionsfile
		optfile=varargin{1};
		if (exist(optfile) ~= 2)
			message('ERROR: Argument has to be existing and valid options file!');
			return 
		end
		
		[dirs, opts] = read_options_from_file(optfile);
		
		if isempty(dirs)
			message('ERROR: directory options [dirs] not provided in options file!');
			return
		end
		if isempty(opts)
			opts = get_default_options;
			message('Processing options [opts] not provided in options file - setting default options.');			
		end	
		
		if ~isfield(dirs,'filelistfile')
			message('ERROR: filelistfile in optionsfile not set! Provide filename for filelist!');
			return
		end
		
		filelist=dirs.filelistfile;		
	else
		
		filelist=varargin{1};
		
		if (iscellstr(filelist))
			% filelist is given directly as MATLAB cell string array
		elseif (isdir(filelist))
			% checkdir()
			message('Processing of directory contents currently not supported!');
			return
		elseif (exist(filelist) ~= 2)
			message('ERROR: Argument 1: Filelist invalid or not found.');
			return
		end 
		
		if (nargin==2)	% in this case 2nd argument ought to be optionsfile
			optfile=varargin{2};
			if (isstruct(optfile))
				message('ERROR: When providing 2 arguments only, 2nd argument must be filename of options file!');
				return 
			elseif (exist(optfile) ~= 2)
				message('ERROR: Argument 2: options file does not exist!');
				return 
			end
			
			[dirs, opts] = read_options_from_file(optfile);
			
			if isempty(dirs)
				message('ERROR: directory options [dirs] not provided in options file!');
				return
			end
			if isempty(opts)
				opts = get_default_options;
				message('Processing options [opts] not provided in options file - setting default options.');			
			end		
		else	% nargin == 3
			dirs=varargin{2};
			opts=varargin{3};
			
			if ~isstruct(dirs)
				if (exist(dirs) ~= 2)
					message('ERROR: Argument 2: options file for directory options does not exist!');
					return 
				end
				
				[dirs, dummy] = read_options_from_file(dirs);
			
				if isempty(dirs)
					message('ERROR: loading of directory options [dirs] failed!');
					return
				end
			end
	
			if ~isstruct(opts)
				if (exist(opts) ~= 2)
					message('ERROR: Argument 3: options file with processing options [opts] does not exist!');
					return 
				end
				
				[dummy, opts] = read_options_from_file(opts);
			
				if isempty(dirs)
					message('ERROR: loading of processing options [opts] failed!');
					return
				end
			end	
		
		end	% if (nargin==2)
        
        % if projectname is not set, derive it from filelist filename
        if ~isfield(dirs,'projectname') || isempty(dirs.projectname)
            % this however is not possible, when filelist is a cell arry
            if (iscellstr(filelist))
                message('ERROR: projectname in [dirs] not set!');
                return
            else
                b = max(findstr(filelist,filesep)) + 1;
                if isempty(b), b=1; end
                e = max(findstr(filelist,'.')) - 1;
                dirs.projectname = filelist(b:e);
            end
        end
	end
	
	filelist
	dirs
	opts
	
	% CALL FEATURE EXTRACTION
	process_filelist(filelist, dirs, opts);
	
	% when called from M2K, use exit to return
	if (~iscellstr(filelist))
		if ~isempty(findstr(filelist,'mirex.opt')), exit; end
	end
   end
end  % if nargin




function opts = get_default_options
    % what to extract/write
    opts.extract_rp = 1;
    opts.extract_pha = 0;   % Phase
    opts.return_full_complex_rp = 0;
    opts.extract_ssd = 1;
    opts.extract_tssd = 1;
    opts.extract_rh = 1;
    opts.extract_trh = 1;
    opts.extract_mvd = 1;  % modulation frequency variance descriptor (vertical variance on Rhythm Pattern)
    opts.extract_bpm = 1;
    opts.write_log = 0;

    % pre-processing options
    opts.resample = 0;              % do resampling before processing: 0 = no, > 0 = resampling frequency in Hz
    
    % processing options
    opts.skip_leadin_fadeout = 1;   % how many sample windows to skip at the beginning and the end
    opts.step_width = 3;            % >=1! each step_width'th sample window is analyzed
    
    opts.n_bark_bands = 24;         % 15 or 20 or 24 (for 11, 22 and 44 kHz audio respectively)
    opts.mod_ampl_limit = 60;
    
    opts.include_DC = 0;

    % enable/disable parts of feature extraction 
    opts.spectral_masking = 0;      % [S3]
    opts.transform_db = 1;          % [S4] advisable only to turn off when [S5] and [S6] are turned off too
    opts.transform_phon = 1;        % [S5] if disabled, sone_transform will be disabled too
    opts.transform_sone = 1;        % [S6] only applies if opts.transform_phon = 1
    opts.fluctuation_strength_weighting = 1;  % [R2] Fluctuation Strength weighting curve
    opts.blurring = 1;              % [R3] Gradient+Gauss filter

    

function [dirs, opts] = read_options_from_file(filename)
   
    dirs=[];
    opts=[];
    opts=get_default_options; % to set [new] fields that are not set in options file
    sec=0; 	% section in options file
    l=0;	% line number
    
    fid=fopen(filename);
    while 1    	
        line = fgetl(fid);
        if ~ischar(line), break, end;
	l = l + 1;
	
	if (~isempty(line))
	    % cut off comments
		p=findstr(line,';');
		if ~isempty(p), line=line(1:p(1)-1); end
		p=findstr(line,'%');
		if ~isempty(p), line=line(1:p(1)-1); end
	end
	
        
	if ((~isempty(line)) && (line(1) ~= ';')) % omit blank lines and lines that start with ;
                   
	    if (line(1) == '[')
            if (strcmp(line(2:5),'dirs')), sec=1;
            elseif (strcmp(line(2:5),'opts')), sec=2; 
            else
                message(['ERROR: Line ',num2str(l),': Invalid section (line starting with ''['') in options file!']);
                return
            end
	    else
		% seperate options key and value
	    p=findstr(line,'=');
        if (length(p) ~= 1)
			message(['ERROR: Line ',num2str(l),': Either 0 or more than one ''='' sign found in options file!']);
			return
        end
		
		opt_id=line(1:p(1)-1);
        opt_id=strrep(opt_id,' ',''); % remove blanks        
		opt_val=line(p(1)+1:end);
        
		if (sec==1)
			dirs=setfield(dirs,opt_id,opt_val);
		elseif 	(sec==2)
            opt_val=strrep(opt_val,' ',''); % remove blanks
			opt_val=str2num(opt_val);
			opts=setfield(opts,opt_id,opt_val);
		end
	    end
	    
        end
    end
    fclose(fid);

    
function [bool,opts] = check_options(opts)
	
	bool = 0;
	
	% Initialize options that have not been provided
	%if ~isfield(opts,'do_resampling'), opts.do_resampling=0; end

    if ~isfield(opts,'include_DC'), opts.include_DC=0; end
    
	if (opts.step_width < 1)
		message('ERROR: option step_width must be 1 or higher!');
		return
	end
	
	bool = 1;   % return true
	
	
function [bool,dirs] = check_dirs(dirs)	
	
	bool = 0;

	if ~isfield(dirs,'binpath'), dirs.binpath=''; end
    if ~isfield(dirs,'input_basedir'), dirs.input_basedir=''; end % in that case let's hope that full path filenames are provided in filelist
    if ~isfield(dirs,'output') || isempty(dirs.output), dirs.output='.'; end
    
	if ~isempty(dirs.binpath)
		if ~strcmp(dirs.binpath(end),filesep), dirs.binpath(end+1) = filesep; end 	
	end
	if ~isempty(dirs.input_basedir)
        if ~strcmp(dirs.input_basedir(end),filesep), dirs.input_basedir(end+1) = filesep; end 
        
        if ~(isdir(dirs.input_basedir))
		    message('ERROR: input_basedir directory does not exist!');
		    return
        end
	end
    
	if ~strcmp(dirs.output(end),filesep), dirs.output(end+1) = filesep; end 
		
	if ~(isdir(dirs.output))
		message('ERROR: output directory does not exist!');
		return
	end

	bool = 1;   % return true
	
	
	    
% DO FEATURE EXTRACTION for multiple files

function process_filelist(filelist, dirs, opts)
    global extr_version;  % to include program version info from above in output files

    errors = 0;     % count no. of errors during processing

	if (nargin < 3) opts = get_default_options; end	
	
	if ~iscell(filelist)	% READ FILELIST from file
		filelistfile = filelist;
		clear filelist;
		n = 0;
		fid=fopen(filelistfile);
		if (~isnumeric(fid) || (fid == -1)),
			message('ERROR: process_filelist: Parameter filelist: filelist-file could not be opened!');
			return
		else
			while 1
				line = fgetl(fid);
				if ~ischar(line), break, end;
				
				n = n + 1; 
                
                % if line contains TAB, suppose that those are additional columns with class information etc.
                t=findstr(line,char(9));
                if ~isempty(t), line = line(1:t(1)-1); end
                
				filelist{n} = line;          
			end
		end
		fclose(fid);	
	end	
	
	[bool,opts] = check_options(opts);
	if (bool == 0) return; end
	
	[bool,dirs] = check_dirs(dirs);
	if (bool == 0) return; end

	n_files = size(filelist,2);
	
	out.log = [dirs.output,dirs.projectname,'.log'];
	%out.matlab_workspace = [dirs.output,dirs.projectname,'.mat'];     
        %out.svml = [dirs.output,dirs.projectname,'.svm'];
        %out.svml_titles = [dirs.output,dirs.projectname,'.svm.txt'];           
    
	if (opts.write_log)
		diary (out.log); diary on;
	end	
	
	time_start = clock;
	message(['== BEGIN OF FEATURE EXTRACTION at ',datestr(time_start,31)]);
	message('-----------------------------------------------------');
	message(['INPUT  DIR: ',dirs.input_basedir]);
	message(['OUTPUT DIR: ',dirs.output]);
	message(['Checking input files: ',num2str(n_files),' files']);
	
	bytes_total = 0;
	bytes_done = 0;
	
	for cur_file = 1:n_files,
		input_file = [dirs.input_basedir,filelist{cur_file}];
		if exist(input_file) == 2
			d = dir(input_file);
			bytes_total = bytes_total + d.bytes;   % NOTE: FE duration estimation will be inaccurate when different file formats (WAV/MP3/AU) are mixed
		else
			message(['ERROR: File does not exist!:']);
			message(input_file);
			% continue anyway, File will be skipped due to error later
		end
	end
	
	message(['Total bytes: ',num2str(bytes_total),' bytes in all files.']);
	message(['   ']);
	
	message(' * Writing output file headers...');
    
    % INITIALIZE SETTINGS for different FILE FORMATS
    
    n_file_formats = 7;

    % RP:   Rhythm Pattern
    % SSD:  Statistical Spectrum Descriptor
    % RH:   Rhythm Histogram
    % BPM:  Beats per Minute
    % TSSD: temporal evolution of SSD
    % TRH:  temporal evolution of RH
    % MVD:  Modulation Frequency Variance Descriptor
    
	tid{1} = 'rp';  writeyn(1)=opts.extract_rp;     vecdimx{1} = opts.mod_ampl_limit;    vecdimy{1} = opts.n_bark_bands; 
    tid{2} = 'ssd'; writeyn(2)=opts.extract_ssd;    vecdimx{2} = get_n_statistical_feat; vecdimy{2} = opts.n_bark_bands; 
    tid{3} = 'rh';  writeyn(3)=opts.extract_rh;     vecdimx{3} = opts.mod_ampl_limit;    vecdimy{3} = 1;
    tid{4} = 'bpm'; writeyn(4)=opts.extract_bpm;    vecdimx{4} = 1;    					 vecdimy{4} = 1;
    tid{5} = 'tssd';writeyn(5)=opts.extract_tssd;   vecdimx{5} = get_n_statistical_feat; vecdimy{5} = opts.n_bark_bands*get_n_statistical_feat; 
    tid{6} = 'trh'; writeyn(6)=opts.extract_trh;    vecdimx{6} = get_n_statistical_feat; vecdimy{6} = opts.mod_ampl_limit; 
    tid{7} = 'mvd'; writeyn(7)=opts.extract_mvd;    vecdimx{7} = get_n_statistical_feat; vecdimy{7} = opts.mod_ampl_limit; 
	
	for i=1:n_file_formats
	  if (writeyn(i))
		outfile = [dirs.output,dirs.projectname,'.',lower(tid{i})];
        message(['   - ',tid{i},' ',outfile]);		
		
		fid = fopen(outfile,'w');
		fprintf(fid,'$TYPE vec\n');
        fprintf(fid,['$DATA_TYPE audio-',tid{i},'\n']);    %content type header
        fprintf(fid,['$DATA_DIM ',num2str(vecdimy{i}),'x',num2str(vecdimx{i}),'\n']);    %content type header
        fprintf(fid,['$EXTRACTOR ',extr_version,'\n']);
		fprintf(fid,'$XDIM %d\n',n_files);
		fprintf(fid,'$YDIM 1\n');
		fprintf(fid,'$VEC_DIM %d\n',vecdimx{i}*vecdimy{i});
		fclose(fid);			
	  end
	end
	
	for cur_file = 1:n_files,
	   try   % in case an error occurs, skip file and proceed to the next one (see catch block below)                                             
		input_file = [dirs.input_basedir,filelist{cur_file}];	

		time_elapsed = etime(clock,time_start);
		
		if (bytes_done==0),
			message([' - Loading file ',num2str(cur_file),' of ', num2str(n_files)]);
		else
			t_expect = round(time_elapsed/bytes_done*(bytes_total - bytes_done));
			message([' - Loading file ',num2str(cur_file),' of ', num2str(n_files),' - Estimated time remaining: ', hminsecstr(t_expect)]);
		end
		
		message(['   ',filelist{cur_file}]);
		
		% CALL FEATURE EXTRACTION
		feat = process_audiofile(input_file, opts);
	
		% WRITE RESULT VECTORS TO VECTOR FILES	
		for i=1:n_file_formats
	  	  if (writeyn(i))
            feat_type = lower(tid{i});
            outfile = [dirs.output,dirs.projectname,'.',feat_type];
            
            if (~isfield(feat,feat_type))
                message([' ERROR: Feature vector ',upper(tid{i}),' was not returned from feature extraction!']);
                % as file does not need to be skipped, just this vector format, don't use error function, but count error
                errors=errors + 1;
            else
                vector = getfield(feat,feat_type);
			
       			if (vector==0)
				    message([' ERROR: ',upper(tid{i}),' vector returned 0-vector from feature extraction!']);
                    % as file does not need to be skipped, just this vector format, don't use error function, but count error
                    errors=errors + 1;
                else
			        message([' * Writing to ',upper(tid{i}),' feature vector file...']);
			
			        fid = fopen(outfile,'a');
		            for j=1:size(vector,2)    
            			fprintf(fid,'%f ',vector(j)); 
        		    end
        		    fprintf(fid,'%s\n',filelist{cur_file});         
        		    fclose(fid);   
                end
		    end % if (~isfield)
          end % if (writeyn(i))
    	end % for
		
		direntry = dir(input_file);
    	bytes_done = bytes_done + direntry.bytes;     

	   catch
		message(lasterr);
		message('Skipping file and proceeding with next file.');
        errors=errors + 1;
	   end  % end try/catch
	end  % end for
		
	time_end = clock;
	time_proc = hminsecstr(etime(time_end, time_start));
	
	message(['== END OF FEATURE EXTRACTION at ',datestr(time_end,31),' - Duration: ',time_proc,' (hour:min:sec)']);
    
    if (errors > 0)
        message(['There were ',num2str(errors),' errors during feature extraction. You might examine the log and output files.']);
    end
		
	if (opts.write_log)
		diary off;
	end
   
	
	
function output_file = decode_mp3(input_file);
    
    % Detecting if EXTERNAL programs exist:
    % As shell return values differ among the programs, and even vary between Windows and Linux versions,
    % this is the somewhat unconventional solution: The amount of output text is measured, and if it is
    % high enough (> 150), the program is said to be existent. 'Not found' error strings contain 42 
    % characters in Linux and 85 characters in a German WindowsXP

    %output_file = [dirs.output,dirs.projectname,'.tmp.wav'];    
    output_file = [tempname,'.wav'];
    dirs.binpath='';
    
    [rval,outtext]=system([dirs.binpath,'mpg123']);   
         
    if size(outtext,2) > 150
        message(['   Decoding mp3 to wav with mpg123 ...']);                 
        system([dirs.binpath,'mpg123 -q -w "',output_file,'" "',input_file,'"']);	
    else
    
    	[rval,outtext]=system([dirs.binpath,'lame --version']);  
    
    	if size(outtext,2) > 150
            message(['   Decoding mp3 to wav with lame ...']); 
            system([dirs.binpath,'lame --quiet --decode "',input_file,'" "',output_file,'"']);            
        else    
	        message('ERROR: Neither LAME nor MPG123 for MP3 decoding found! Cannot decode mp3!');
            output_file=[];
            return
        end
    end

function output_file = decode_ape(input_file);
     
    % Detecting if EXTERNAL programs exist:
    % - see above: decode_mp3
  
    output_file = [tempname,'.wav'];
    dirs.binpath='';
    
    [rval,outtext]=system([dirs.binpath,'flac']);   
         
    if size(outtext,2) > 150
        message(['   Decoding ape to wav with mac -d ...']);                 
        system([dirs.binpath,'mac "',input_file,'" "',output_file,'" -d']);	
    else     
	    message('ERROR: mac binary for decoding .ape files not found! Cannot decode ape!');
		output_file=[];
		return
    end

     
function output_file = decode_flac(input_file);
    
    % Detecting if EXTERNAL programs exist:
    % - see above: decode_mp3
  
    output_file = [tempname,'.wav'];
    dirs.binpath='';
    
    [rval,outtext]=system([dirs.binpath,'flac']);   
         
    if size(outtext,2) > 150
        message(['   Decoding flac to wav with flac -d ...']);                 
        system([dirs.binpath,'flac -d -f -o "',output_file,'" "',input_file,'"']);	
    else     
	    message('ERROR: flac binary for decoding flac files not found! Cannot decode flac!');
		output_file=[];
		return
    end

   
    
function output_file = decode_ogg(input_file);
    
    % Detecting if EXTERNAL programs exist:
    % - see above: decode_mp3
  
    output_file = [tempname,'.wav'];
    dirs.binpath='';
    
    [rval,outtext]=system([dirs.binpath,'oggdec']);   
         
    if size(outtext,2) > 150
        message(['   Decoding ogg to wav with oggdec ...']);                 
        system([dirs.binpath,'oggdec -Q -o "',output_file,'" "',input_file,'"']);	
    else     
	    message('ERROR: program oggdec for decoding ogg files not found! Cannot decode ogg!');
		output_file=[];
		return
    end


function output_file = resample(input_file, input_extension, hz);

    if  (strcmp(lower(input_extension),'au'))
        output_file = [tempname,'.au'];
    else
        output_file = [tempname,'.wav'];
    end
    
    dirs.binpath='';
    
    [rval,outtext]=system([dirs.binpath,'sox -h']);
    
    if size(outtext,2) > 150
        message(['   Resampling to ',hz,' Hz/mono with sox ...']);
        system([dirs.binpath,'sox "',input_file,'" -r ',hz,' -c 1 "',output_file,'"']);
    else
        message('ERROR: Cannot do resampling: program sox was not found!');
        output_file=[];
        return;
    end

    
    
% DO FEATURE EXTRACTION for 1 file

function feat = process_audiofile(input_filename, opts)

	%global CONST_filt;
	
	if (nargin < 2) opts = get_default_options; end

	% return values in case of an error	
	if (opts.extract_rp), feat.rp = 0; end
	if (opts.extract_ssd), feat.ssd = 0; end
	if (opts.extract_rh), feat.rh = 0; end
	if (opts.extract_bpm), feat.bpm = 0; end

	% PRE-CALCULATIONS FROM WAVE HEADER DATA
	input_extension = lower(input_filename(max(findstr(input_filename,'.'))+1:end));

    %message(input_filename);   
    
    temp_file = 0;
    
	% call external decoder to decode mp3 files to temporary wav file
	if (strcmp(input_extension,'mp3')) 
		input_filename = decode_mp3(input_filename);
		if isempty(input_filename), return; end
		input_extension = 'wav';
		temp_file = 1;
	end

  	% call external decoder to decode ape files to temporary wav file
	if (strcmp(input_extension,'ape')) 
		input_filename = decode_ape(input_filename);
		if isempty(input_filename), return; end
		input_extension = 'wav';
		temp_file = 1;
	end
 
  	% call external decoder to decode flac files to temporary wav file
	if (strcmp(input_extension,'flac')) 
		input_filename = decode_flac(input_filename);
		if isempty(input_filename), return; end
		input_extension = 'wav';
		temp_file = 1;
	end
    
 	% call external decoder to decode ogg files to temporary wav file
	if (strcmp(input_extension,'ogg')) 
		input_filename = decode_ogg(input_filename);
		if isempty(input_filename), return; end
		input_extension = 'wav';
		temp_file = 1;
	end
    
    % do resampling of audio if option is set
    if (opts.resample > 0)
        input_filename = resample(input_filename, input_extension, num2str(opts.resample));
        if isempty(input_filename), return; end
		temp_file = 1;
    end    
    
    filesize = 0;
    sample_freq = 0;
	% Read WAVE header information
	if (strcmp(input_extension,'au')) 
		[y Fs] = audioread(input_filename);
		sample_freq = Fs;
	else	% assume wav
		[y Fs] = audioread(input_filename);
		sample_freq = Fs;
    end		
	
    filesize = size(y);
	duration =  filesize(1)/sample_freq;	

	message(['   ','Processing Wave Signal: Size: ',num2str(filesize(1)),' samples, ',minsecstr(duration),' min:sec']); 
    
	if filesize(2) > 1
        	message(['   ','WARNING: File has ',num2str(filesize(2)),' channels. Averaging channels.']);
    	end		
		
	% segment_size should always be ~6 sec, fft_window_size should always be ~ 23ms
	if (sample_freq <= 11025),	segment_size = 2^16;   fft_window_size = 256;
	elseif (sample_freq <= 22050)	segment_size = 2^17;   fft_window_size = 512;
	else 				segment_size = 2^18;   fft_window_size = 1024;  % assume 44100 Hz	
	end
	
	pre.fft_window_size = fft_window_size;
	
	% calculate frequency values on y-axis (for bark scale calculation)
	pre.freq_axis=sample_freq/fft_window_size*(0:fft_window_size/2);	

	% modulation frequency x-axis (after 2nd fft)
	mod_freq_res = 1 / (segment_size/sample_freq);  % resolution of modulation frequency axis (0.17 Hz)             
	mod_freq_axis = mod_freq_res*(0:256);		% modulation frequencies along x-axis from index 1 to 257)
	% "the modulation frequencies [...] are in the range from 0 to 43 Hz with an accuracy of 0.17 Hz"
	warning off;	% ignore div/0 for mod_freq_axis(1)
	pre.fluct_curve = 1./(mod_freq_axis/4 + 4./mod_freq_axis);	% fluctuation strength weighting curve along modulation freq. axis
	warning on;
	
	% filter pre-calculation (dep. on n_bark_bands and mod_ampl_limit)
    CONST_filt = [.05 .1 0.25 0.5 1 0.5 .25 0.1 0.05];   % BLUR FILTER
	pre.blur1 = filter2(CONST_filt,eye(opts.n_bark_bands));
	pre.blur1 = pre.blur1./repmat(sum(pre.blur1,2),1,opts.n_bark_bands);
	pre.blur2 = filter2(CONST_filt,eye(opts.mod_ampl_limit)); 
	pre.blur2 = pre.blur2./repmat(sum(pre.blur2,2),1,opts.mod_ampl_limit);
	
	% WAVE PROCESSING ROUTINE - INITIALIZATION
	
	seg_pos = [1 segment_size]; % segment position in wave file (samples) (init with start position)
    	step_width = opts.step_width;	% TODO IDEA: automatic step width
	skip_seg   = opts.skip_leadin_fadeout;
    	%end_of_wave = filesize(1) - (segment_size * opts.skip_leadin_fadeout);
    
    if ((opts.skip_leadin_fadeout > 0) | (step_width > 1))        		
		if (duration < 45)
			% if file is too small, don't skip leadin/fadeout and set step_width to 1
			step_width = 1;
			skip_seg = 0;
			%end_of_wave = filesize(1);
			message(['   ','WARNING: File is < 45 sec. Ignoring step width > 1 and skip leadin/fadeout options.']);
		else
			seg_pos = seg_pos + segment_size*skip_seg;	% advance by # of skip_seg segments (i.e. skip lead_in)
		end
	end       

	% calculate number of segments
	n_segments =  floor( (floor( (filesize(1)-(skip_seg*2*segment_size)) / segment_size ) - 1 ) / step_width ) + 1;
	
	message_progress(0);

	% WAVE PROCESSING ROUTINE - SEGMENT LOOP	
	
	%while (seg_pos(2) + win_num*meta_window_samples) < end_of_wave, 
	for seg_id = 1:n_segments		% seg_id = segment index
		
		% Read 6 sec. WAVE segments
				
		if  (strcmp(input_extension,'au')) 
			wav = auread(input_filename,seg_pos); 
		else
			wav = wavread(input_filename,seg_pos); 
		end	% auread/wavread return values in the range [-1,+1]		
		
		if (filesize(2) > 1)   % when there are more channels (e.g. stereo file), take the mean of all channels
			wav = mean(wav,2);
		end    
		
		% Extract features from segment
		
		%[feat_rp_seg, feat_ssd_seg, feat_rh_seg] = extract_features(wav, opts, pre);
		feat_seg = extract_features(wav, opts, pre);
		
		if (seg_id == 1)	% init feature matrix for all segments of an audio file
			if (opts.extract_rp), feat_rp_list = zeros(n_segments, size(feat_seg.rp(:),1)); end
			if (opts.extract_rh | opts.extract_trh), feat_rh_list = zeros(n_segments, size(feat_seg.rh(:),1)); end
			if (opts.extract_ssd | opts.extract_tssd), feat_ssd_list= zeros(n_segments, size(feat_seg.ssd(:),1)); end
            if (opts.extract_mvd), feat_mvd_list = zeros(n_segments, size(feat_seg.mvd(:),1)); end
		end
		
		if (opts.extract_rp), feat_rp_list(seg_id,:) = feat_seg.rp(:)'; end
		if (opts.extract_rh | opts.extract_trh), feat_rh_list(seg_id,:) = feat_seg.rh(:)'; end
		if (opts.extract_ssd | opts.extract_tssd), feat_ssd_list(seg_id,:) = feat_seg.ssd(:)'; end
		if (opts.extract_mvd), feat_mvd_list(seg_id,:) = feat_seg.mvd(:)'; end
        
		message_progress(seg_id/n_segments*100);
		seg_id = seg_id + 1;
		seg_pos = seg_pos + segment_size*step_width;
				
		%if (debug.playwav), wavplay(wav/(0.0875*2^15),sample_freq); end
	end % for
		
	if (temp_file)
		delete(input_filename);
	end
	
	if (n_segments < 1)
		message(['ERROR: number of segments calculation for audio file was ',num2str(n_segments),'.']);
		return
	else
		% build SUMMARIZATION of all segments of an audio file
		% (note: do not omit dimension (1), in case that list had only 1 entry, median would be calculated from row!)
		if (opts.extract_rp), feat.rp = median(feat_rp_list,1); end
		if (opts.extract_ssd), feat.ssd = mean(feat_ssd_list,1); end
        if (opts.extract_rh), feat.rh = median(feat_rh_list,1); end
		if (opts.extract_mvd), feat.mvd = mean(feat_mvd_list,1); end
        
        if (opts.extract_tssd)  % temporal statistical features over SSD features
            tssd = statistical_features(feat_ssd_list');          
            feat.tssd = tssd(:)'; 
        end
        
        if (opts.extract_trh)  % temporal statistical features over SSD features
            trh = statistical_features(feat_rh_list');          
            feat.trh = trh(:)'; 
        end
        		
		if (opts.extract_bpm)	 % [BPM] derive beats per minute
			max_ampl = find(feat.rh==max(feat.rh));  % find index (indices) of modulation frequency with maximum energy    
			feat.bpm = round(max_ampl(1) * mod_freq_res * 60);
		end
		
		message_progress(100);
		message([' complete. (',num2str(n_segments),' segments)']);  
	end





% EXTRACT features of a 6 second segment

function feat = extract_features(wavsegment, opts, precalc)

	% adjust hearing threshold	
	wavsegment = 0.0875*wavsegment*2^15;  % TODO: investigate this
	
	% [S1]
	%disp('S1')
	spectrogr = fft_real_hann(wavsegment, precalc.fft_window_size);
	
	% [S2]
	%disp('S2')
	matrix = transform_bark(spectrogr, precalc.freq_axis, opts.n_bark_bands);
	
	if (opts.spectral_masking)	% [S3]
		%disp('S3');
		matrix = spectral_masking(matrix);
	end
	
	if (opts.transform_db)		% [S4]
		%disp('S4');
		matrix = transform_db(matrix);
	end
	
	if (opts.transform_phon)	% [S5]
		%disp('S5');
		matrix = transform_phon(matrix);
	end
	
	if (opts.transform_sone)	% [S6]
		%disp('S6')
		matrix = transform_sone(matrix);
	end
	
	if (opts.extract_ssd | opts.extract_tssd)	% [SSD] Statistical Spectrum Descriptor
		%disp('SSD')
        % compute statistics/variances on the 24 critical bands
        % Note that the full matrix with 511 values per band is used.
		feat.ssd = statistical_features(matrix);
    end
	
    
    if (opts.extract_rp | opts.extract_pha | opts.extract_rh | opts.extract_trh | opts.extract_mvd)
        %remainder has to be processed only for these features
        
        % [R1]
        %disp('R1')

        feature_part_xaxis1 = 1:opts.mod_ampl_limit;    % take first (opts.mod_ampl_limit) values of fft result including DC component
        feature_part_xaxis2 = 2:opts.mod_ampl_limit+1;  % leave DC component and take next (opts.mod_ampl_limit) values of fft result
        if (opts.include_DC)  feature_part_xaxis_rp = feature_part_xaxis1; else feature_part_xaxis_rp = feature_part_xaxis2; end

        % GET complex rhythm_patterns

        rhythm_patterns = modulation_ampl(matrix) ./ 256;  % TODO: investigate why 256 here?

        if (opts.return_full_complex_rp)
            feat.rp_complex = rhythm_patterns;
        end

        if (opts.extract_rp | opts.extract_mvd)
            % take only a part of the rhythm patterns FFT result as actual feature set
            feat.rp = abs(rhythm_patterns(:,feature_part_xaxis_rp));	
        end
        
        if (opts.extract_mvd)
            % compute statistics/variances between the critical bands alongside the 60 different modulations frequencies
            % (matrix has to be transposed, as it is transposed back in the statistical_features function)
            feat.mvd = statistical_features(feat.rp');
        end

        if (opts.extract_pha)   % Phase information
            feat.pha = angle(rhythm_patterns(:,feature_part_xaxis1));
        end

        if (opts.extract_rh | opts.extract_trh)	% [RH] Rhythm Histogram Features
            %disp('RH');
            % ALWAYS skip DC component 
            feat.rh = sum(abs(rhythm_patterns(:,feature_part_xaxis2)));
        end

        if (opts.extract_rp & opts.fluctuation_strength_weighting)	% [R2]
            %disp('R2');
            % feature_part_xaxis has to be considered in the indices of fluct_curve
            feat.rp =  fluctuation_strength_weighting(feat.rp, precalc.fluct_curve(feature_part_xaxis_rp));
        end

        if (opts.extract_rp & opts.blurring)	% [R3]
            %disp('R3');
            feat.rp = blur_patterns(feat.rp, precalc.blur1, precalc.blur2);
        end
	
    end % if (opts.extract_rp | opts.extract_pha | opts.extract_rh)
    
	%message('extraction fin.');
	
    

% [S1] spectrogram: real FFT with hanning window and 50 % overlap
function spectrogr = fft_real_hann(wavsegment, fft_window_size)

	segment_size = size(wavsegment,1);
		
	n_iter = segment_size/fft_window_size*2 - 1;   % # of iterations with 50 % overlap
	
	w = hann(fft_window_size);	
	spectrogr = zeros(fft_window_size/2 + 1, n_iter);
	
	idx = 1:fft_window_size;
	for i=1:n_iter,	% stepping through the wave segment, building spectrum for each window		
		spectrogr(:,i) = periodogram(wavsegment(idx),w);
		idx = idx + fft_window_size/2;
	end

	
% [S2] transform to Bark scale (sum frequency bands to critical bands)
function matrix = transform_bark(spectrogr, freq_axis, n_bands)
	global CONST_bark;
	
	matrix = zeros(n_bands, size(spectrogr,2));
	matrix(1,:) = sum(spectrogr(find(freq_axis < CONST_bark(1)),:),1);
	for i=2:n_bands
		matrix(i,:) = sum(spectrogr(find(freq_axis >= CONST_bark(i-1) & freq_axis < CONST_bark(i)),:));
	end
	

% [S3] Spreading function for Spectral Masking
function matrix = spectral_masking(matrix)
	global CONST_spread;
	
	spread = CONST_spread(1:size(matrix,1),:); % TODO: investigate - this does not work when n_bark_bands <> 24
	
	matrix = spread * matrix;		
	
	
% [S4] transform to log. deciBel scale
function matrix = transform_db(matrix)

	matrix(matrix<1)=1;
	matrix = 10*log10(matrix);
	
	
% [S5] Phon transformation: equal loudness curves
function matrix = transform_phon(matrix)	% TODO: investigate
		
	global CONST_phon;
	global CONST_loudn_freq;
	global CONST_eq_loudness;
	global CONST_loudn_bark;

            % number of bark bands, matrix length in time dim
            [n_bands, t] = size(matrix);
    
            % DB-TO-PHON BARK-SCALE-LIMIT TABLE
            % introducing 1 level more with level(1) being infinite
            % to avoid (levels - 1) producing errors like division by 0
            
            %%%table_dim = size(CONST_loudn_bark,2);
            table_dim = n_bands;
            cbv = [repmat(Inf,table_dim,1) CONST_loudn_bark(:,1:n_bands)'];
            phons = [0 CONST_phon];                      
           
           % init lowest level = 2
           levels = repmat(2,n_bands,t); 
           
           for lev = 2:7,
               db_thislev = repmat(cbv(:,lev),1,t);      
               levels(find(matrix(:,1:t) > db_thislev)) = lev + 1;                
           end
           
           % the matrix 'levels' stores the correct Phon level for each datapoint
           
           cbv_ind_hi = sub2ind([table_dim 7],repmat([1:table_dim]',1,t),levels);
           cbv_ind_lo = sub2ind([table_dim 7],repmat([1:table_dim]',1,t),levels-1); % TODO: optimize either hi or lo away
           % for all data values we have now the indices in cbv for looking up the equal loudness values for the correct level
           
           ifac = (matrix(:,1:t) - cbv(cbv_ind_lo)) ./ (cbv(cbv_ind_hi) - cbv(cbv_ind_lo)); % interpolation factor % OPT: pre-calc diff
           ifac(levels==2) = 1;    % keeps the upper phon value;
           ifac(levels==8) = 1;    % keeps the upper phon value;
           if (sum(levels==8) > 0)
               %message('  ERROR: dB limit in Phon calculation!');
	       error('dberror','  ERROR: dB limit in Phon calculation!');
           end
           matrix(:,1:t) = phons(levels - 1) + ifac .* ((phons(levels) - phons(levels - 1))); % OPT: pre-calc diff
           
        


% [S6] Sone transformation: specific loudness sensation
function matrix = transform_sone(matrix)
	
	idx = matrix>=40;
	matrix(idx)  = 2.^((matrix(idx)-40)/10);
	matrix(~idx) = (matrix(~idx)/40).^2.642;
	
	
% [R1] FFT to get Rhythm Patterns
function rhythm_patterns = modulation_ampl(matrix)  % returns complete rhythm_patterns in complex form

	% find an appropriate (i.e. calculation efficient) fft size (power of 2)
	fft_size = 2^(nextpow2(size(matrix,2)));
	
	rhythm_patterns = zeros(size(matrix,1), fft_size);
	
	for b = 1:size(matrix,1)
		%fftcomplex = fft(matrix(b,:), fft_size);
		%rhythm_patterns(b,:) = abs(fftcomplex/256);   % TODO: investigate why 256 here?
        rhythm_patterns(b,:) = fft(matrix(b,:), fft_size);
	end

	
% [R2] Fluctuation Strength weighting curve
function rhythm_patterns = fluctuation_strength_weighting(rhythm_patterns, fluct_curve)
	% fluct_curve was pre-calculated in initialization after reading Wave header 
	% (only a part is considered here; see extract_features())	

	for b = 1:size(rhythm_patterns,1)
		rhythm_patterns(b,:) = rhythm_patterns(b,:) .* fluct_curve;
	end
	

% [R3] Gradient+Gauss filter	
function rhythm_patterns = blur_patterns(rhythm_patterns, blur1, blur2)

	for i=2:size(rhythm_patterns,2),
                rhythm_patterns(:,i-1) = abs(rhythm_patterns(:,i) - rhythm_patterns(:,i-1));
        end
        rhythm_patterns = blur1 * rhythm_patterns * blur2;
	

% [SSD] Statistical Spectrum Descriptor (7 statistical features for each ROW of the matrix)
function statmat = statistical_features(matrix)  
    
    %Note: according to the Master's thesis of Bernhard Pflugfelder, variance and skewness are little useful!

    mat=matrix';
    statmat=zeros(size(mat,2),7);
    statmat(:,1)=mean(mat)';
    statmat(:,2)=var(mat)';
    warning off;
    statmat(:,3)=skewness(mat)';    % Warning: div/0 Warnings appear and result is NaN when an entire ROW is 0
    statmat(:,4)=kurtosis(mat)';    % Warning: div/0 Warnings appear and result is NaN when an entire ROW is 0    
    warning on;
    statmat(:,5)=median(mat)';
    statmat(:,6)=min(mat)';         
    statmat(:,7)=max(mat)';
 
    % set NaN values to 0
    statmat(find(isnan(statmat)))=0;
    
    % normalize columns
    %for i = 1:7
    %    statmat(:,i) = statmat(:,i)/max(statmat(:,i));
    %end 	

% IMPORTANT: edit this function to return always the correct number of statistical features in function above
function n = get_n_statistical_feat      
    n = 7;
        

	
	
% --------------------------------------------------------------------
%                       HELPER FUNCTIONS
% --------------------------------------------------------------------

	
function message(msg)
    disp(msg);
    %global MSG_commandline MSG_gui handles;
    %if (MSG_commandline), disp(msg); end
    %if (MSG_gui), buffertext(msg,handles); end
    
    
function message_progress(percent)
    %global MSG_commandline MSG_gui handles;
    %if (MSG_commandline)
        if (percent > 0), fprintf('\b\b\b\b\b\b\b\b\b'); end
        fprintf('  %5.1f %%',percent);          
    %end   
    %if (MSG_gui) set(handles.txt7,'string',['   Progress: ',num2str(percent, '%3.1f'),' %']); end

function s=minsecstr(seconds)       % (by tml) converts seconds to a min:sec string
	min = floor(seconds/60);
	sec = round(seconds - min*60);
	s = sprintf('%2.2d:%2.2d',min,sec);

function s=hminsecstr(seconds)       % (by tml) converts seconds to a hour:min:sec string
	min = floor(seconds/60);
	sec = round(seconds - min*60);
	hour = floor(min/60);
	min = min - hour*60;
	s = sprintf('%d:%2.2d:%2.2d',hour,min,sec);		
	
