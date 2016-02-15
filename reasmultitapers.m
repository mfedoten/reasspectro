function varargout = reasmultitapers(varargin)

% check inputs and outputs
narginchk(1,15);
nargoutchk(1,6);
% distribute inputs
[sig,win,Nw,ovlap,nfft,fs,opts]=parse_inpts(varargin{:});
opts = check_opts(opts);

% get spectrograms

% reassign spectrograms

% average spectrograms

% return outputs