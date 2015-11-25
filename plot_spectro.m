function plot_spectro(t,f,S,varargin)
% Function to plot spectrograms
%
% INPUTS
% t : time of spectrogram;
% f : frequency;
% S : spectrogram matrix;
%
% OPTIONS (as structure)
% Nw    : length of the window (used to plot COI);
% tReal : original time (used to plot COI);
% type  : type of plots: 'image', 'contour', 'pcolor'
% ncont : number of countour levels, if contour is chosen as display method
% flim  : max. frequency to be displayed;
% font  : size of the ticks font;
% hax   : handle to axes, where you want to plot it, otherwise creates new
%         figure

% get options into structure
opts = struct(varargin{:});

% get default font
fs_default = get(0,'DefaultAxesFontSize');
% set font sizes
if isfield(opts,'font')
    % smallest font is for axes
    if nargin > 6 && font_size
        set(0,'DefaultAxesFontSize',font_size);
    end
    fs_ticks  = get(0,'DefaultAxesFontSize');
    fs_labels = fs_ticks + 4;
end

% use provided figure, if asked
if isfield(opts,'hax')
    ha = opts.hfig;
else
    figure('Units','Centimeters');
    fpos = get(gcf,'Position');
    fpos = [0.6*fpos(1:2) 1.2*fpos(3) 1.2*fpos(3)];
    set(gcf,'Position',fpos);
    
    ha = axes('Units','Centimeters','Position',[0.1*fpos(3:4) 0.8*fpos(3:4)]);
end
set(ha,'Render','painters');

% chose plotting method
if ~isfield(opts,'type') || ~any(strcmpi(opts.type,{'image','contour','pcolor'}))
    opts.type = 'image';
end

switch type
    case 'image'
        
% plot full TFR
if nargin > 5 && flim 
    imagesc(t,f(f<=flim),S(f<=flim,:));
else
    imagesc(t,f,S);
end

% normal y-direction
axis xy;
axis square;
% decent tick lengths
set(gca,'ticklength',.5*get(gca,'ticklength'));

% anotate the plots
ylabel('Frequency', 'FontSize', fs_labels);
xlabel('Time', 'FontSize', fs_labels);
if nargin > 7 && ischar(titleString)
    fs_title  = fs_labels + 2;
    title(titleString, 'FontSize', fs_title);
end

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
fs = 1/(tt(2)-tt(1));
% color to plot COI
ccoi = 'w';
if Nw/2/fs > t(1)
    disp_coi(gca,tt,Nw,fs,fs/2,ccoi,min(S(:)));
end

% set default font size back
set(0,'DefaultAxesFontSize',fs_default);
end


function disp_coi(hAx,tReal,Nw,fs,fmax,cc,minLvl)
axes(hAx);
hold on;
% add ones due to imagesc properties
hPatch = patch(tReal(1) + [-1 Nw/2/fs Nw/2/fs -1],...
    [-5 -5 fmax+1 fmax+1], minLvl);
hatchfill(hPatch, 'cross', 45, 10, cc);
plot([Nw/2/fs Nw/2/fs],[-1 fmax+1],cc,'linewidth',1.5);
hPatch = patch(tReal(end) - [-1 Nw/2/fs Nw/2/fs -1],...
    [-5 -5 fmax+1 fmax+1], minLvl);
hatchfill(hPatch, 'cross', 45, 10, cc);
plot(tReal(end)-[Nw/2/fs Nw/2/fs],[-1 fmax+1],cc,'linewidth',1.5);
end


