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