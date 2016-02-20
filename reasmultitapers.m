function varargout = reasmultitapers(varargin)
%REASMULTITAPERS multitaper spectrogram reassignment
%  RS = REASMULTITAPERS(SIG,NSEQ,TAPERS,OVLAP,NFFT,FS,'PROPERTYNAME',PROPERTYVALUE)
%  RS = REASMULTITAPERS(SIG,NSEQ,TAPERS,OVLAP,NFFT,FS,PROPERTYSTRUCTURE)
% [RS,S] = REASMULTITAPERS(...)
% [RS,FNEW,TNEW] = REASMULTITAPERS(...)
% [RS,FNEW,TNEW,S,FORIG,TORIG] = REASMULTITAPERS(...)
%
% INPUT:
% SIG    - analysed signal of length N points.
% NSEQ   - length of tapers in samples. If no NSEQ is provided, then default
%          value is NSEQ = N/10.
% TAPERS - taper functions, can be in one of the two forms:
%          (1) As a vector [NW K], where NW is a time half bandwidth product and
%          K is a number of tapers. NW must be strictly less than NSEQ/2 and K
%          is usually chosen to be 2*NW-1 (to provide good energy
%          concentration).
%          (2) Slepian sequences (e.g. obtained from MATLAB dpss function).
%          Should be a matrix with NSEQ rows and K columns, where K is a number
%          of tapers.
%          If nothing is specified default value is TAPERS=[3 5].
% OVLAP  - number of overlapping points between two adjacent moving windows. Is 
%          used to compensate energy loss at the ends of smoothing windows.
%          Should be specified as an integer smaller than the length of the
%          window. If no OVLAP is given, the default value is used to obtain a
%          50%.
% NFFT   - number of Fourier transform points. If no input is provided, the
%          default value is chosen as next power of 2 greater then the length of
%          the window.
% FS     - sampling frequency in Hz. Default is FS = 1.
%
% PROPERTIES:
% By default all properties are unset (except 'mean' property). Property can be
% passed as name-value pair or as a structure. If field requires logical
% true/false input, every non-empty, finite and not 'NaN' input is treated as
% logical true, e.g. true, 'yes', 1 etc. will be validated as logical true.
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
%'mean'    : Method of averaging spectrograms obtained with multitapers into one
%            time-frequency representation. The possible values are:
%            'mean' (for arithmetic mean),'geom' (geometric mean),'min' (takes
%            the min value among all spectrograms),'median' (uses median). For
%            more information on different averaging methods see [3].
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
% This function computes reassigned version of the multitaper spectrogram. The
% algorithm is based on [3], some parts, like interpolation part are based on
% [2].
% 
% The main idea is the same as in reassigned spectrogram (see help reasspecgram
% for more details), the only difference is that we compute it for each taper
% and then average over tapers. It is more computationally expensive (the method
% requires computations of 3 STFT for each taper), but it has advantage of lower
% varianse and betters localisation than conventional spectrogram.
% 
% References
% 1. Auger, F. et al. Time-Frequency Reassignment and Synchrosqueezing: An
% Overview. IEEE Signal Processing Magazine 30, 32?41 (2013). 
% 2. Fulop, S. A. & Fitz, K. Algorithms for computing the time-corrected
% instantaneous frequency (reassigned) spectrogram, with applications. The
% Journal of the Acoustical Society of America 119, 360 (2006).
% 3. Xiao, J. & Flandrin, P. Multitaper Time-Frequency Reassignment for
% Nonstationary Spectrum Estimation and Chirp Enhancement. IEEE Transactions on
% Signal Processing 55, 2851?2860 (2007).
%
% See also: reasspecgram.m
%
% Copyright Mariia Fedotenkova, 2016, INRIA Nancy.
% Licensed for use under GNU General Public License, Version 2.  See LICENSE for
% details.

% check inputs and outputs
narginchk(1,15);
nargoutchk(1,6);
% distribute inputs
[sig,Nseq,taps,ovlap,nfft,fs,opts]=parse_inpts(varargin{:});
opts = reassignment_check_opts('tapers',opts);

% generate taper or check provided
tapers = reassignment_get_tapers(taps,Nseq);
% number of tapers
K = size(tapers,2);

% pad signal with half of the window on both sides in order to avoid edge
% effects and to have time starting from zero, not half of the window
if opts.pad
    sig = reassignment_pad_signal(sig,Nseq,opts.pad);
end

% get the number of rows and columns in the output matrix
tcol = fix((length(sig)-ovlap)/(Nseq-ovlap));
if ~isreal(sig)     % analytical signal
    frow = nfft;
else                % real signal
    if rem(nfft,2)
        frow = (nfft+1)/2;
    else
        frow = nfft/2+1;
    end
end

% initialize empty matrices for spectrogram and reassigned spectrogram
Stapers  = zeros(frow,tcol,K);
RStapers = zeros(frow,tcol,K);

% now we'll have to loop over all tapers
%(unless one day I come up with something better)
for k = 1:K
    % tapers are nothing else but window functions
    win = tapers(:,k);
    
    % do the reassignment
    [RStapers(:,:,k),fnew,tnew,Stapers(:,:,k),forig,torig] = ...
        reassignment_core(sig,win,ovlap,nfft,fs,opts);
end

% average spectrograms
S  = reassignment_get_mean(Stapers,opts.mean);
RS = reassignment_get_mean(RStapers,opts.mean);


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

end


%--------------------------------------------------------------------------
function [sig,Nseq,tapers,ovlap,nfft,fs,opts]=parse_inpts(varargin)
sig = varargin{1};
N  = length(sig);
% if length of tapers is not provided, set Nw=N/10
if nargin < 2 || isempty(varargin{2})
    Nseq = floor(N/10);
else
    Nseq = varargin{2};
end
% if no tapers provided, set [NW K] = [3 5]
if nargin < 3 || isempty(varargin{3})
    tapers = [3 5];
else
    tapers = varargin{3};
end
% if no overlap provided, use 50% by default
if nargin < 4 || isempty(varargin{4})
    ovlap = floor(0.5*Nseq);
else
    ovlap = varargin{4};
    if ovlap >= Nseq
        error('Overlap should be less or equal than window length');
    end

end
% if no nfft provided, use next power of 2 greater than the length of the
% window. If Nfft greater than the length of the signal, use the latter
% instead.
if nargin < 5 || isempty(varargin{5})
    nfft = 2^nextpow2(Nseq);
else
    nfft = varargin{5};
    if nfft > 2^nextpow2(N)
        nfft = 2^nextpow2(N);
        warning('NFFT is more than signal length. Using signal length instead.');
    end
end
% if no sampling rate provided use length of the signal by default -> the
% tend = 1s.
if nargin < 6 || isempty(varargin{6})
    fs = 1;
    opts = struct;
elseif nargin >= 6
    fs = varargin{6};
    if ~isscalar(fs) || fs < 0
        error('Sampling rate should be positive scalar');
    end
    % the rest (if any) are options for new spacing, turn it into structure
    if nargin > 6
        try
            opts = struct(varargin{7:end});
        catch
            error('Specify properties as one or more name-value pairs.')
        end
    else
        opts = struct;
    end
end
end
