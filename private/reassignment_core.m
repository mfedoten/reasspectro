function [RS,fnew,tnew,S,forig,torig] = reassignment_core(sig,win,ovlap,nfft,fs,opts)
% This functions does actually the main job in spectrogram reassignment:
% - creates additional windows;
% - computes three STFTs;
% - does the reassignment;
% - computes new vectors.
%
% Inputs and outputs are described in reasspecgram.m (see for more info).
%
% (C) Mariia Fedotenkoava 2016.

% get window size
Nw = length(win);

% construct additional windows for reassignment
[Twin,Dwin] = reassignment_get_windows(win,fs);

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