function varargout = reasspecgram(varargin)
%REASSPECGRAM spectrogram reassignment
%  RS = REASSPECGRAM(SIG,WIN,OVLAP,NFFT,FS,'PROPERTYNAME',PROPERTYVALUE)
%  RS = REASSPECGRAM(SIG,WIN,OVLAP,NFFT,FS,PROPERTYSTRUCTURE)
% [RS,S] = REASSPECGRAM(...)
% [RS,FNEW,TNEW] = REASSPECGRAM(...)
% [RS,FNEW,TNEW,S,FORIG,TORIG] = REASSPECGRAM(...)
%
% INPUT:
% SIG   - analysed signal of length N points.
% WIN   - smoothing window. If WIN is a vector, the signal SIG is divided into
%         segments of length equal to the length of the window, and each
%         segment is filtered with WIN. If WIN is an integer, then SIG is
%         divided into segments of length equal to the specified integer and
%         filtered with Hamming window. If no window is provided, then Hamming
%         window of length N/10 is used by default.
% OVLAP - number of overlapping points between two adjacent windows. Is used to
%         compensate energy loss at the ends of smoothing windows. Should be
%         specified as an integer smaller than the length of the window. If no
%         OVLAP is given, the default value is used to obtain a 50% overlap.
% NFFT  - number of Fourier transform points. If no input is provided, the
%         default value is chosen as next power of 2 greater then the length of
%         the window.
% FS    - sampling frequency in Hz.
%
% PROPERTIES:
% By default all properties are unset. Property can be passed as name-value pair
% or as a structure. If field requires logical true/false input, every non-empty,
% finite and not 'NaN' input is treated as logical true, e.g. true, 'yes', 1
% etc. will be validated as logical true.
% 'psd'    : true | false
%            If true, power spectral density (PSD) is reassigned. If not,
%            squared absolute value of short-time Fourier transform is used.
% 'pad'    : 'zeros' | 'const' | 'periodic' |'symmetric'
%            If specified before processing a signal will be padded according to
%            the chosen method. Padding can be done with zeros ('zeros'); with
%            constant  equal to the first and last values of the signal
%            ('const'); 'periodic' continues signal in both directions
%            periodically; 'symmetric' reflects signal symmetrically on both
%            sides. If padding is chosen it is added on both slides, the length
%            of the padding is half of the window length.
% 'crop'   : Sometimes there might be extreme values in the final time-frequency
%            representation. You can "crop" such values, by specifying
%            percentile, above which the values should not be considered
%            (analogy with outliers), value of 99.5-99.9% is reasonable.
% 'size'   : You can change the size of the output TFR matrix.
%            Specify a new size as two-element vector: number of rows
%            (frequency) and number of columns (time). Not recommended, because
%            the number of elements in the TFR matrix is fixed.
% 'step'   : Similar to 'size' parameter, but instead of specifying new size you
%            provide new frequency and time sampling.
% 'interp' : Interpolate neighbouring points. Can be useful when using new
%            size/sampling, because it increases the number of points.
%            Otherwise, it is not recommended to use it, because the idea of
%            reassignment is lost like that.
%
% OUTPUT:
% RS    : A matrix with reassigned spectrogram. Here, rows are frequencies and
%         columns are time points.
% S     : A matrix with original (not reassigned) spectrogram. Here, rows are
%         frequencies and columns are time points.
% FNEW  : Frequency vector corresponding to rows in reassigned spectrogram
%         matrix.
% TNEW  : Time vector corresponding to columns in reassigned spectrogram matrix.
% FORIG : Frequency vector corresponding to rows in spectrogram matrix.
% TORIG : Time vector corresponding to columns in spectrogram matrix.
%
%
% This function computes reassigned version of the spectrogram. The algorithm is
% based on [1], some parts are based on [2], like interpolation part. The idea
% is to first compute conventional spectrogram, then find optimal (in a sense of
% energy) time and frequency positions and reassigns values in the spectrogram
% to this new positions.
%
% You have a choice between reassigning spectrogram either as short-time Fourier
% transform (STFT) or power spectrum density (STFT normalized by window energy
% and sampling rate). Controlled by 'pad' property.
%
% You can pad signal before any processing, it can be done by providing a type
% of padding in 'pad' property.
%
% Sometimes, extreme values (outliers) appear on a spectrogram (reassigned and
% original), for instance, when using too many overlapping points between the
% segments. You can 'crop' them by providing to a desired value of percentile of
% the data, above which all the values will be set to the highest values. It
% works similarly to thresholding. To use this feature specify percentile as
% 'crop' property.
%
% It is possible to obtain spectrogram on finer grid with this method. For that,
% you have to specify in the properties three additional parameters: new steps
% or new dimensions of time and frequency, specify if it is new steps or
% dimensions. Usage of this feature is not recommended, because even though it
% is possible to increase number of points in the output matrix the original
% number of points is not changing. So, you might obtain very sparse matrix. To
% partially solve it you can use interpolation (set 'interp' parameter), which
% interpolates additional points with values of its neighbouring points. A use
% of this feature is also not advised, because it decreases the sharpness of
% reassigned spectrogram and contradicts the whole idea of the reassignment.
% 
% References
% 1. Auger, F. et al. Time-Frequency Reassignment and Synchrosqueezing: An
% Overview. IEEE Signal Processing Magazine 30, 32?41 (2013). 
% 2.Fulop, S. A. & Fitz, K. Algorithms for computing the time-corrected
% instantaneous frequency (reassigned) spectrogram, with applications. The
% Journal of the Acoustical Society of America 119, 360 (2006).
%
% Copyright Mariia Fedotenkova, 2015, INRIA Nancy.
% Licensed for use under GNU General Public License, Version 2.  See LICENSE for
% details.


% check inputs and outputs
narginchk(1,15);
nargoutchk(1,6);
% distribute inputs
[sig,win,Nw,ovlap,nfft,fs,opts]=parse_inpts(varargin{:});
opts = check_opts(opts);


% make sure window is a column vector
win = win(:);
% construct additional windows for reassignment
[Twin,Dwin] = reassignment_get_windows(win,fs);


% pad signal with half of the window on both sides in order to avoid edge
% effects and to have time starting from zero, not half of the window
if opts.pad
    sig = reassignment_pad_signal(sig,Nw,opts.pad);
end


% compute three STFTs
Sw = reassignment_get_stft(sig,win,ovlap,nfft);
Stw = reassignment_get_stft(sig,Twin,ovlap,nfft);
Sdw = reassignment_get_stft(sig,Dwin,ovlap,nfft);


% nr. of time points and frequency bins
[frow,tcol] = size(Sw);


% squared short-time Fourier transform
S = abs(Sw).^2;
% if the output is psd
if opts.psd
    % normalized spectrogram by window's energy
    S = S/(win'*win);
    % divide over sampling frequency to get PSD (Power/freq)
    S = S/fs;
end
% if the signal is real, we take only half of the spectrogram, so we have
% to multiply this half by two, except DC and Nyquist
if frow ~= nfft
    S = [S(1,:); 2*S(2:end-1,:); S(end,:)];
end


% original time and frequency vectors (indices)
% MATLAB indexation drives me crazy, I start from zero here
Forig = (0:frow-1)';
Torig = (0:(tcol-1))*(Nw-ovlap);
% original time and frequency vectors (seconds and Hz)
% !don't take half of the window for now, it will be returned in the end
forig = Forig*fs/nfft;
torig = Torig/fs;


% get new reassigned vectors of frequencies (in Hz) and times (in s)
[fhat,that] = reassignment_get_displacements(Sw,Stw,Sdw,torig,forig);


% create vectors of time and frequency with higher spacing or just by
% rounding fhat,that to the closest bin in original vectors
[Snew,Fhat,That,fnew,tnew] = reassignment_get_new_vectors(fhat,that,...
        forig,torig,S,opts);


% reassign spectrogram values to new locations
% turn all matrices into column vectors
That = That(:);
Fhat = Fhat(:);
Snew = Snew(:);
sz = [length(fnew) length(tnew)];


if opts.interp
    RS = zeros(sz);
    alpha = That - floor(That);
    beta  = Fhat - floor(Fhat);
    for k = 1:length(Snew)
        RS(floor(Fhat(k)),floor(That(k))) = RS(floor(Fhat(k)),floor(That(k))) + ((1-alpha(k)) * (1-beta(k)) * Snew(k));
        RS(ceil(Fhat(k)),floor(That(k)))  = RS(ceil(Fhat(k)),floor(That(k))) + ((1-alpha(k)) * (beta(k)) * Snew(k));
        RS(floor(Fhat(k)),ceil(That(k)))  = RS(floor(Fhat(k)),ceil(That(k))) + ((alpha(k)) * (1-beta(k)) * Snew(k));
        RS(ceil(Fhat(k)),ceil(That(k)))   = RS(ceil(Fhat(k)),ceil(That(k))) + ((alpha(k)) * (beta(k)) * Snew(k));
    end
else
    RS = accumarray([round(Fhat) round(That)],Snew,sz);
end


% crop reassigned and conventional sectrograms to have only values below
% specified percentile
if opts.crop
    RS = crop_matrix(RS,opts.crop);
    S  = crop_matrix(S,opts.crop);
end

% add half of the window to time vector if no padding was used
if ~opts.pad
    torig = torig + floor(Nw/2)/fs;
    tnew  = tnew  + floor(Nw/2)/fs;
end


% ------------------------ distribute the outputs -------------------------
switch nargout
    case 1
        varargout = {RS};
    case 2
        varargout = {RS,S};
    case 3
        varargout = {RS,fnew,tnew};
    case 6
        varargout = {RS,fnew,tnew,S,forig,torig};
    otherwise
        error('Wrong number of outputs. See help for more information.')
end



%--------------------------------------------------------------------------
function [sig,win,Nw,ovlap,nfft,fs,opts]=parse_inpts(varargin)
sig = varargin{1};
N  = length(sig);
% if window is not defined, use hamming window of length=N/10
if nargin < 2 || isempty(varargin{2})
    Nw = floor(N/10);
    win = hamming(Nw);
elseif isscalar(varargin{2})
    Nw = varargin{2};
    win = hamming(Nw);
else
    win = varargin{2}(:);
    Nw = length(win);
end
% if no overlap provided, use 50% by default
if nargin < 3 || isempty(varargin{3})
    ovlap = floor(0.5*Nw);
else
    ovlap = varargin{3};
    if ovlap >= Nw
        error('Overlap should be less or equal than window length');
    end

end
% if no nfft provided, use next power of 2 greater than the length of the
% window. If Nfft greater than the length of the signal, use the latter
% instead.
if nargin < 4 || isempty(varargin{4})
    nfft = 2^nextpow2(Nw);
else
    nfft = varargin{4};
    if nfft > 2^nextpow2(N)
        nfft = 2^nextpow2(N);
        warning('NFFT is more than signal length. Using signal length instead.');
    end
end
% if no sampling rate provided use length of the signal by default -> the
% tend = 1s.
if nargin < 5 || isempty(varargin{5})
    fs = N;
    opts = struct;
elseif nargin >= 5
    fs = varargin{5};
    if fs < 0 || ~isscalar(fs)
        error('Sampling rate should be positive scalar');
    end
    % the rest (if any) are options for new spacing, turn it into structure
    if nargin > 5
        try
            opts = struct(varargin{6:end});
        catch
            error('Specify properties as one or more name-value pairs.')
        end
    else
        opts = struct;
    end
end


%--------------------------------------------------------------------------
function opts = check_opts(opts)
% first check that all the fields are valid
val_flds = {'new_sampling','step','size','interp','psd','crop','pad'};
opt_flds = fieldnames(opts);
not_valid = opt_flds(~ismember(opt_flds,val_flds));
if ~isempty(not_valid)
    error('"%s" not a valid property.\n', not_valid{:});
end

% check step/size property (if any)
if isfield(opts,'size') && isfield(opts,'step')
    % we should chose only one way to specify output matrix size
    error('"size" and "step" properties are mutually exclusive, you should chose only one.');
elseif isfield(opts,'size')
   % check the 'size' property
   if any(size(opts.size)~=[1 2]) && any(size(opts.size)~=[2 1])
       % 'size' should be two-element vector
       error('New size should be a two-element vector.')
   end
   if any(~isfinite(opts.size)) || ~isreal(opts.size) || ischar(opts.size)
       % it also should be finite and real
       error('Size should contain only real finite values');
   end
   opts.new_sampling = true;
elseif isfield(opts,'step')
   % the same goes for 'step'
   if any(size(opts.step)~=[1 2]) && any(size(opts.step)~=[2 1])
       error('New sampling step should be a two-element vector.')
   end
   if any(~isfinite(opts.step)) || ~isreal(opts.step) || ischar(opts.step)
       error('Sampling step should contain only real finite values');
   end
   opts.new_sampling = true;
elseif ~isfield(opts,'size') && ~isfield(opts,'step')
   % if neither step or size are specified, then we don't want to change
   % size of the output matrix
   opts.new_sampling = false;
end

% check interpolation property (if any)
if isfield(opts,'interp')
    % if 'interp' property was set check that it meets all requirements
    if isempty(opts.interp)||any(isnan(opts.interp))||all(~opts.interp)
        % we don't want interpolation if it is false, empty or NaN
        opts.interp = false;
    else
        % otherwise turn it to logical true
        opts.interp = true;
    end
else
    % if it wasn't passed, set it to false
    opts.interp = false;
end

% check 'return PSD' property (if any)
if isfield(opts,'psd')
    if isempty(opts.psd)||any(isnan(opts.psd))||all(~opts.psd)
        % we don't want PSD if it is false, empty or NaN
        opts.psd = false;
    else
        % otherwise turn it to logical true
        opts.psd = true;
    end
else
    opts.psd = false;
end

% check crop property (if any)
if isfield(opts,'crop')
    if isempty(opts.crop)||~isscalar(opts.crop)
        % we don't want to crop the matrix
        opts.crop = false;
    elseif opts.crop<0 || opts.crop>100
        % percentile should be between 0 and 100%
        error('Percentile should be specified as a number between 0 and 100%%.');
    end
else
    opts.crop = false;
end

% check padding property (if any)
if isfield(opts,'pad')
    if isempty(opts.pad)
        % we don't want any padding
        opts.crop = false;
    elseif ~ischar(opts.pad)
        % should be a string
        error('Specify padding type as a string');
    end
else
    opts.pad = false;
end
