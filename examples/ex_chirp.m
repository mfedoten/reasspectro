clear all; close all; clc;
rng('default');
addpath('..');
set(0,'DefaultAxesFontSize',12)

% construct a test signal
N  = 1024;          % number of points
T  = 1;             % total time
fs = N/T;           % sampling rate
dt = 1/fs;          % sampling step
t  = 0:dt:T-dt;     % time vector
x  = chirp(t,0.1*fs,t(end),0.4*fs,'quadratic')';
x  = x + 0.2*randn(size(x));

% parameters of spectrogram
df    = 10;                    % frequency resolution
Nw    = floor(fs/df);           % window length
win   = hamming(Nw);            % window function
ovlap = floor(.9*Nw);           % number of overlapping points
nfft  = 2^nextpow2(Nw);         % number of FT points

%% defaults
[RS,f_reas,t_reas,S,f_sp,t_sp] = reasspecgram(x,win,ovlap,nfft,fs);

% plot
figure('Units','Centimeters');
fpos = get(gcf,'Position');
fpos = [0.6*fpos(1) fpos(2) 2*fpos(3) fpos(3)];
set(gcf,'Position',fpos);

h1 = subplot(1,2,1,'Units','centimeters','Position',[2 3 fpos(4)-5 fpos(4)-5]);
h2 = subplot(1,2,2,'Units','centimeters','Position',[fpos(4)+2.5 3 fpos(4)-5 fpos(4)-5]);

plot_spectro(t_sp,f_sp,S,'tReal',t,'Nw',Nw,'hax',h1);
title('Spectrogram','FontSize',16);
plot_spectro(t_reas,f_reas,RS,'tReal',t,'Nw',Nw,'hax',h2);
title('Reassigned spectrogram','FontSize',16);

export_fig(gcf, 'example_defaults', '-png','-transparent','-painters');
