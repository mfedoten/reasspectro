function [fhat,that] = reassignment_get_displacements(Sw,Stw,Sdw,torig,forig)
% This function computes frequency and time displacements.
%
% INPUT:  Sw     - STFT with the original window;
%         Stw    - STFT with the time window;
%         Sdw    - STFT with time derivative of the window;
%         S      - Spectrogram normalized;
%         Forig  - original frequency vector;
%         Forig  - original time vector;
%         thres  - values lower than threshold will not be reassigned.
% OUTPUT: Fhat - matrix, of the same size as Sw, where each points defines
%                "reassigned" frequency for each point in Sw;
%         That - the same as Fhat, but for time;
%
%
% (C) Mariia Fedotenkova 2016.


% find where STFT is non-zero
idx = find(Sw);

% compute time and frequency displacements
Stw(idx) = real(Stw(idx)./Sw(idx));             % for time
Sdw(idx) = imag(Sdw(idx)./Sw(idx));             % for frequency

% nr. of time points and frequency bins
[frow,tcol] = size(Sw);

% turn freq. (Hz) and time (s) vectors into grid
that = ones(frow,1)*torig;
forig = forig(:);
fhat = forig*ones(1,tcol);

% construct reassigned frequency and time vectors
that = that + Stw;
fhat = fhat - Sdw;
