function varargout = reasmultitapers(varargin)
%REASMULTITAPERS multitaper spectrogram reassignment
%  RS = REASMULTITAPERS(SIG,WIN,OVLAP,NFFT,FS,'PROPERTYNAME',PROPERTYVALUE)
%  RS = REASMULTITAPERS(SIG,WIN,OVLAP,NFFT,FS,PROPERTYSTRUCTURE)
% [RS,S] = REASMULTITAPERS(...)
% [RS,FNEW,TNEW] = REASMULTITAPERS(...)
% [RS,FNEW,TNEW,S,FORIG,TORIG] = REASMULTITAPERS(...)
%
% INPUT:
% SIG    - analysed signal of length N points.
% NSEQ   - length of tapers in samples. If no NSEQ is provided, then default
%          value is NSEQ = N/10.
% TAPERS - taper functions, can be in one of the two forms:
%          (1) Slepian sequences (e.g. obtained from MATLAB dpss function).
%          Should be a matrix with NSEQ rows and K columns, where K is a number
%          of tapers.
%          (2) As a vector [NW K], where NW is a time half bandwidth product and
%          K is a number of tapers. NW must be strictly less than NSEQ/2 and K
%          is usually chosen to be 2*NW-1 (to provide good energy
%          concentration).
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

% check inputs and outputs
narginchk(1,15);
nargoutchk(1,6);
% distribute inputs
[sig,Nseq,taps,ovlap,nfft,fs,opts]=parse_inpts(varargin{:});
opts = reassignment_check_opts(opts);

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
S  = reassignment_get_mean(Stapers,opt.mean);
RS = reassignment_get_mean(Stapers,opt.mean);


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
