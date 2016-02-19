function tfrAvg = reassignment_get_mean(tfr,method)
% Function to average multiple TFRs obtained with different tapers into a single
% TFR, averaging is done according to a provided method. For more information on
% averaging methods see [1].
%
% INPUT:
% tfr    : TFRs to be averaged. Should have size of Nf x Nt x K, where Nf is the
%          number of frequency bins in the TFR, Nt is number of time points in
%          the TFR and K is the number of tapers. Averaging along K.
% method : Averaging method. The possible values are:
%          'mean','geom','min','median'.
% 
% OUTPUT:
% tfrAvg : final (averaged) TFR.
%
% References:
% 1. Xiao, J. & Flandrin, P. Multitaper Time-Frequency Reassignment for
% Nonstationary Spectrum Estimation and Chirp Enhancement. IEEE Transactions on
% Signal Processing 55, 2851?2860 (2007).
%
% (C) Mariia Fedotenkoava 2016.

switch(method)
    case 'mean'
        tfrAvg = mean(tfr,3);
    case 'geom'
        tfrAvg = exp(mean(log(max(tfr,1e-12)),3));
    case 'min'
        tfrAvg = min(tfr,[ ],3);
    case 'median'
        tfrAvg = median(tfr,3);
end
end