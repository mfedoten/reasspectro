function plot_spectro(t,f,S,varargin)
% Function to plot spectrograms
%
% INPUTS
% t : time of spectrogram;
% f : frequency;
% S : spectrogram matrix;
%
% OPTIONS (as structure or name-value pair)
% Nw     : length of the window; used to plot COI (if not provided COI is not 
%          plotted);
% tReal  : original time; used to plot COI (if not provided COI is not plotted);
% type   : type of plots: 'image', 'contour', 'pcolor', default is 'image';
% ncont  : number of countour levels (if contour is chosen as display method),
%          if not provided, nr. of contours is chosen automatically;
% flim   : frequency limits, given as two-element vctor: [fmin fmax] OR as a
%          scalar value: fmax. If you chose to plot frequency axis in dBs you
%          should provide limits in dBs as well;
% font   : ticks font size, axes font is computed automatically as font+2;
% hax    : handle to axes, where you want to plot it, otherwise creates new
%          figure;
% dbFreq : plot frequency in dB (skips f=0Hz);
% dbPow  : plot Power in dB (assign eps to zero values);
% colmap : colormap.
% 
% Note: if you want to plot frequency in dBs, you shouldn't use 'image' plot
% type, use 'contour' or 'pcolor' instead.
%
% (C) Mariia Fedotenkova 2016.


% get options into structure
opts = struct(varargin{:});

% check frequency for inf values
if ~all(isfinite(f))
    error('Frequency vector contains Inf values');
end

% transform frequencies to dB
if isfield(opts,'dbFreq') && ~isempty(opts.dbFreq) && opts.dbFreq
    S(f==0,:) = [];
    f(f==0) = [];
    f = 10*log10(f);
elseif abs(sum(diff(diff(f)))) > eps
    opts.dbFreq = true;    
else
    opts.dbFreq = false;
end

% transform power to dB
if isfield(opts,'dbPow') && ~isempty(opts.dbPow) && opts.dbPow
    S(S==0) = eps;
    S = 10*log10(S);
end

% chose plotting method
if ~isfield(opts,'type') || isempty(opts.type) || ~any(strcmpi(opts.type,...
        {'image','contour','pcolor'}))
    % set default plotting method
    if opts.dbFreq
        opts.type = 'pcolor';
    else
        opts.type = 'image';
    end
    % if non-linear frequencise + image -> switch to pcolor
elseif opts.dbFreq && all(strcmpi(opts.type,'image'))
    warning(['Image will not give good results with non-linear frequencies.',...
        ' Switching to pcolor']);
    opts.type = 'pcolor';
end
opts.type

% use provided figure, if asked
if isfield(opts,'hax') && ~isempty(opts.hax)
    ha = opts.hax;
else
    figure('Units','Centimeters');
    fpos = get(gcf,'Position');
    fpos = [0.6*fpos(1:2) 1.2*fpos(3) 1.2*fpos(3)];
    set(gcf,'Position',fpos);
    ha = axes('Units','Centimeters','Position',[0.1*fpos(3:4) 0.8*fpos(3:4)]);
end

switch opts.type
    case 'image'
        set(gcf,'Render','painters');
        imagesc(t,f,S,'Parent',ha);
        axis xy;
    case 'contour'
        if isfield(opts,'ncont') && ~isempty(opts.ncount)
            contourf(ha,t,f,S,opts.ncont,'EdgeColor','None');
        else
            contourf(ha,t,f,S,'EdgeColor','None');
        end
    case 'pcolor'
        pcolor(ha,t,f,S); 
        shading flat;
end

% set Y-limits; plot until max freq., if specified
if isfield(opts,'flim') && ~isempty(opts.flim)
    if isscalar(opts.flim)
        ylim([f(1) opts.flim]);
    elseif length(opts.flim)==2
        ylim([opts.flim(1) opts.flim(2)]);
    end
else
    ylim([f(1) f(end)]);
end

% set desired colormap
if isfield(opts,'colmap') && ~isempty(opts.colmap)
    colormap(ha,opts.colmap);
end

% change font sizes
if isfield(opts,'font') && ~isempty(opts.font)
    % smallest font is for axes ticks
    fs_ticks = opts.font;
    set(ha,'FontSize',fs_ticks);
else
    fs_ticks = get(ha,'FontSize');
end
fs_labels = fs_ticks + 2;

% anotate the plots
ylabel('Frequency (Hz)', 'FontSize', fs_labels);
xlabel('Time (s)', 'FontSize', fs_labels);

% Colorbar
pos = get(gca,'Position');
cb = colorbar; pause(0.1);
if verLessThan('matlab','8.4')
    set(cb, 'TickLength', [0 0], 'FontSize', fs_ticks);
    poscb = get(cb, 'Position');
    set(cb, 'Position', [poscb(1)+poscb(3) poscb(2) poscb(3) poscb(4)]);
else
    set(cb, 'TickLength', 0, 'FontSize', fs_ticks, 'Box', 'off');
    cb.Ruler.Axle.Visible = 'off';
    cb.Ruler.SecondaryLabel.HorizontalAlignment = 'left';
end
set(gca,'Position',pos);

% cone of influence (plot only if exists)
if isfield(opts,'Nw') && isfield(opts,'tReal')
    fs = 1/(opts.tReal(2)-opts.tReal(1));
    % color to plot COI
    ccoi = 'w';
    if opts.Nw/2/fs > t(1)
        disp_coi(gca,opts.tReal,opts.Nw,fs,[f(1) f(end)],ccoi,min(S(:)));
    end
end

end


function disp_coi(hAx,tReal,Nw,fs,flim,cc,minLvl)
axes(hAx);
hold on;
% add ones due to imagesc properties
hPatch = patch(tReal(1) + [0 Nw/2/fs Nw/2/fs 0],...
    [flim(1) flim(1) flim(2) flim(2)], minLvl);
hatchfill(hPatch, 'cross', 45, 10, cc);
plot([Nw/2/fs Nw/2/fs],flim,cc,'linewidth',1.5);
hPatch = patch(tReal(end) - [0 Nw/2/fs Nw/2/fs 0],...
    [flim(1) flim(1) flim(2) flim(2)], minLvl);
hatchfill(hPatch, 'cross', 45, 10, cc);
plot(tReal(end)-[Nw/2/fs Nw/2/fs],flim,cc,'linewidth',1.5);
end


