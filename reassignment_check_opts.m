function opts = reassignment_check_opts(opts)
% Helper function to check options for reassignment

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