function [varargout] = reassignment_get_new_vectors(fhat,that,forig,torig,S,opt)
% This function creates new frequency and time vectors. It transforms time
% and  frequency vectors returned by reassignment method according to the
% new spacing.
% If last input is not specified, the function rounds the output from
% reassignment procedure to the closest value in original vectors.
% If new spacing is specified, the function creates vectors of new time and
% frequency indices, but with higher sampling.
%
% Reassignment procedure returns vectors of time (in seconds) and
% frequencies (in Hz), which specify "real" position of each point in
% spectrogram. This function first creates the output vectors of time and
% frequencies with desired sampling It also takes care that frequencies
% higher than Fmax and lower than 0 will be wrapped and time larger than
% Tmax and smaller than 0 will be truncated. Then it creates vectors of
% corresponding integer indices. The last steps is to map output vectors
% from reassignment to newly created vectors.
% NOTE: if you input higher sampling, time vector limits will change to
% [0 tend]
%
% INPUT:  fhat, that   - output times and frequencies from reassignment;
%         forig, torig - original times and frequencies;
%         fs           - sampling rate;
%         Inputs required only if new spacing is desired:
%         tend/tshift  - a scalar; if new spasing is not required it's a
%                        time shift of time vector, equal to half of the
%                        window; if new spacing is required it is a
%                        finishing time in seconds;
%         new_spacing  - a two-element vector defining new spacing, first
%                        element for frequency, second - for time.
%         spacing_flag - string, defines how to interpret previous vector.
%         'step' - for new spacing, 'size' - for the size of the output
%         matrix.
% OUTPUT: Fhat, That - new indices for reassignment;
%         fnew,tnew  - new times and frequencies (optional);
%
%
% (C) Mariia Fedotenkova 2016.

if ~opt.new_sampling
    % construct new time and frequency vectors
    fnew = forig;
    tnew = torig;
    % calculate time and frequency steps;
    df = forig(2) - forig(1);
    dt = torig(2) - torig(1);
else
    % get time and frequency steps
    if isfield(opt,'step')
        df = opt.step(1);
        dt = opt.step(2);
    elseif isfield(opt,'size')
        df = (forig(end) - forig(1))/(opt.size(1)-1);
        dt = (torig(end) - torig(1))/(opt.size(2)-1);
    end
    % construct new higher spaced vectors
    fnew = (forig(1):df:forig(end))';
    tnew = torig(1):dt:torig(end);
end

% make sure that spectrogram values (S) at frequencies outside [0 Fend] and
% times outside [t1 tend] are not reassigned, they are set to zero
idx = find(fhat<0 | fhat>fnew(end) | that<tnew(1) | that>tnew(end));
fhat(idx) = 0;
that(idx) = 0;
S(idx) = 0;


% get new reassignment vectros
Fhat = 1+fhat/df;
That = 1+that/dt;


% ----------------------------- return output -----------------------------
switch nargout
    case 3
        varargout = {S,Fhat,That};
    case 5
        varargout = {S,Fhat,That,fnew,tnew};
    otherwise
        error('Wrong number of outputs. Number of output elements can be either 2 or 4.')
end



