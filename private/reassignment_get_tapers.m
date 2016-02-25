function tapers=reassignment_get_tapers(taps,N)
% Function to generate tapers from given parameters or check provided tapers.
%
% INPUT:
% taps : (1) two element vector [NW K], where NW is time half-bandwidth and K is
%        number of tapers;
%        OR
%        (2) matrix with Slepian sequences, where nr. of rows is number of time
%        points and nr. of columns is number of tapers.
% N    : length of each taper in samples.
% 
% OUTPUT:
% tapers : final Slepian sequences.
%
% (C) Mariia Fedotenkoava 2016.

% get size
sz = size(taps);
% it should be two-element vector
if length(sz) > 2
    error(['There''s something wrong with the tapers. ',...
        'Why are they %d-dimensional?'],length(sz));
end

if all(sz == [1 2]) || all(sz==[2 1])
    % if 2-element vector is provided
    if taps(2) > 2*taps(1) - 1
        % print warning if K>2*NW-1
        warning(...
            ['Desired number of tapers (%d) is larger than recommended (%d)',...
            '. Tapers will not be sufficiently concentrated in frequency.'], ...
            taps(2),2*taps(1) - 1);
    end
    % generate tapers
    tapers = dpss(N,taps(1),taps(2));
else    % if dpss are provided
    % tapers' length should match NSEQ
    if max(sz) ~= N
        error('Length of provided tapers do not match sliding window length (NSEQ).');
    end
    % tapers should be NSEQxK matrix
    if sz(1) < sz(2)
        tapers = taps.';
    else
        tapers = taps;
    end        
end

end