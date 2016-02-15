function [RStm,Stm] = reastapers(sig,t,nfft,ovlap,Nh,M,tm)

[h,dh,tt] = hermf(Nh,M,tm);

K = fix((length(sig)-ovlap)/(length(h)-ovlap));
St = zeros(nfft/2+1,K,M) ;
RSt = zeros(nfft/2+1,K,M) ;

for k = 1:M
    fs = 1/(t(2) - t(1));
    
    win  = h(k,:).';
%     Dwin = dh(k,:).';
    Nw   = length(win);

%     if rem(Nw,2)
%         tw = (-(Nw-1)/2:(Nw-1)/2)';  % window is odd
%     else
%         tw = (-Nw/2:Nw/2-1)';        % window is even
%     end
    %     Twin = win.*tt.';
%     Twin = win.*tw;

    [Twin,Dwin] = reassignment_get_windows(win,fs);
    
    % compute three STFTs
    Sw = reassignment_get_stft(sig,win,ovlap,nfft);
    Stw = reassignment_get_stft(sig,Twin,ovlap,nfft);
    Sdw = reassignment_get_stft(sig,Dwin,ovlap,nfft);
    
    
    % nr. of time points and frequency bins
    [frow,tcol] = size(Sw);
    
    
    % squared short-time Fourier transform
    S = abs(Sw).^2;
    % normalized spectrogram by window's energy
    S = S/(win'*win);
    % divide over sampling frequency to get PSD (Power/freq)
    S = S/fs;
    
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
    [Snew,Fhat,That,fnew,tnew] = reassignment_get_new_vectors_tmp(fhat,that,...
        forig,torig,S);
    
    
    % reassign spectrogram values to new locations
    % turn all matrices into column vectors
    That = That(:);
    Fhat = Fhat(:);
    Snew = Snew(:);
    sz = [length(fnew) length(tnew)];
    
    
    RS = accumarray([round(Fhat) round(That)],Snew,sz);
    
    % distribute outputs
    RSt(:,:,k) = RS;
    St(:,:,k) = S;
end

Stm = mean(St,3) ;
RStm = mean(RSt,3) ;