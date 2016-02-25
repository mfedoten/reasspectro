function S = reassignment_get_stft(sig,win,ovlap,nfft)
% This function computes short-time Fourier transform (STFT) to be used in
% spectrogram reassinment method.
%
% INPUT:  signal  - signal to analyze;
%         window  - smoothing sliding window used in STFT;
%         overlap - sliding window overlap, in time points;
%         nfft    - number of Fourier transform points.
% OUTPUT: S - STFT matrix.
%
% Only signal is required parameter, others are optional. If only the 
% signal is provided window is chosen to be hamming window of N/10 length 
% with 50% overlap and nfft is next power of 2 greater than the length of 
% the window. The returned STFT matrix is not normalized.
%
%
% (C) Mariia Fedotenkoava 2016.



% ---------------------------- arguments check ----------------------------
% display help if no input provided
if nargin < 1
    help reassignment_get_stft
    return
else
    % nr. of points in the signal
    sig = sig(:);
    N  = length(sig);
end
% if window is not defined, use hamming window of length=N/10
if nargin < 2 || isempty(win)
    Nw = floor(N/10);
    win = hamming(Nw);
else
    win = win(:);
    Nw = length(win);
end
% if no nfft provided, use next power of 2 greater than the length of the
% window. If Nfft greater than the length of the signal, use the latter
% instead.
if nargin < 3 || isempty(nfft)
    nfft = 2^nextpow2(Nw);
else
    if nfft>N
        nfft=N;
        warning('NFFT is more than signal length. Using signal length instead.');
    end
end
% if no overlap provided, use 50% by default
if nargin < 4 || isempty(ovlap)
    ovlap = floor(0.5*Nw);
else
    if ovlap >= Nw
        error('Overlap should be less or equal than window length');
    end
end
% ------------------------- end of argument check -------------------------


% number of columns in the output matrix (time points)
tcol = fix((N-ovlap)/(Nw-ovlap));
% number of rows in the output matrix (freq. bins)
if ~isreal(sig)
    frow = nfft;
else
    if rem(nfft,2)
        frow = (nfft+1)/2;
    else
        frow = nfft/2+1;
    end
end

% cut signal into segments of window length, taking overlap into account, 
% and stack them into matrix of size [Nw,tcol]
rowindex = (1:Nw)';
colindex = 1 + (0:(tcol-1))*(Nw-ovlap);
X = zeros(Nw,tcol);
X(:) = sig(rowindex(:,ones(1,tcol))+colindex(ones(Nw,1),:)-1);

% multiply each segment with smoothing window and take Fourier transform
S = win(:,ones(1,tcol)).*X;
S = fft(S,nfft);
S = S(1:frow,1:tcol);

