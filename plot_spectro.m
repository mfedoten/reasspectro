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
% flim   : max. frequency to be displayed;
% font   : ticks font size, axes font is computed automatically as font+2;
% hax    : handle to axes, where you want to plot it, otherwise creates new
%          figure;
% colmap : colormap.

% get options into structure
opts = struct(varargin{:});

% use provided figure, if asked
if isfield(opts,'hax') && ~isempty(opt.hax)
    ha = opts.hax;
else
    figure('Units','Centimeters');
    fpos = get(gcf,'Position');
    fpos = [0.6*fpos(1:2) 1.2*fpos(3) 1.2*fpos(3)];
    set(gcf,'Position',fpos);
    
    ha = axes('Units','Centimeters','Position',[0.1*fpos(3:4) 0.8*fpos(3:4)]);
end
set(gcf,'Render','painters');

% set font sizes
if isfield(opts,'font') && ~isempty(opts.font)
    % smallest font is for axes ticks
    fs_ticks = opts.font;
    set(ha,'FontSize',fs_ticks);
else
    fs_ticks = get(ha,'FontSize');
end
fs_labels = fs_ticks + 2;

% chose plotting method
if ~isfield(opts,'type') || ~any(strcmpi(opts.type,{'image','contour','pcolor'}))
    opts.type = 'image';
end

axes(ha);
switch opts.type
    case 'image'
        imagesc(t,f,S);
        axis xy;
    case 'contour'
        if isfield(opts,'ncont') && ~isempty(opts.ncount)
            contourf(t,f,S,opts.ncont,'EdgeColor','None');
        else
            contourf(t,f,S,'EdgeColor','None');
        end
    case 'pcolor'
        pcolor(t,f,S); 
        shading flat;
end

% plot until max freq. if specified
if isfield(opts,'flim') && ~isempty(opts.flim)
    ylim([f(1) opts.flim]);
end

% set desired colormap
if isfield(opts,'colmap') && ~isempty(opts.colmap)
    colormap(opts.colmap);
end

% anotate the plots
ylabel('Frequency (Hz)', 'FontSize', fs_labels);
xlabel('Time (s)', 'FontSize', fs_labels);

% Colorbar
pos = get(gca,'Position');
cb = colorbar; pause(0.5);
if verLessThan('matlab','8.4')
    set(cb, 'TickLength', [0 0], 'FontSize', 11);
    poscb = get(cb, 'Position');
    set(cb, 'Position', [poscb(1)+poscb(3) poscb(2) poscb(3) poscb(4)]);
else
    set(cb, 'TickLength', 0, 'FontSize', 12, 'Box', 'off');
    cb.Ruler.Axle.Visible = 'off';
    cb.Ruler.SecondaryLabel.HorizontalAlignment = 'left';
    cb.FontSize = fs_ticks;
end
set(gca,'Position',pos);

% cone of influence (plot only if exists)
if isfield(opts,'Nw') && isfield(opts,'tReal')
    fs = 1/(opts.tReal(2)-opts.tReal(1));
    % color to plot COI
    ccoi = 'w';
    if opts.Nw/2/fs > t(1)
        disp_coi(gca,opts.tReal,opts.Nw,fs,fs/2,ccoi,min(S(:)));
    end
end

end


function disp_coi(hAx,tReal,Nw,fs,fmax,cc,minLvl)
axes(hAx);
hold on;
% add ones due to imagesc properties
hPatch = patch(tReal(1) + [-1 Nw/2/fs Nw/2/fs -1],...
    [-1 -1 fmax+1 fmax+1], minLvl);
hatchfill(hPatch, 'cross', 45, 10, cc);
plot([Nw/2/fs Nw/2/fs],[-1 fmax+1],cc,'linewidth',1.5);
hPatch = patch(tReal(end) - [-1 Nw/2/fs Nw/2/fs -1],...
    [-1 -1 fmax+1 fmax+1], minLvl);
hatchfill(hPatch, 'cross', 45, 10, cc);
plot(tReal(end)-[Nw/2/fs Nw/2/fs],[-1 fmax+1],cc,'linewidth',1.5);
end


