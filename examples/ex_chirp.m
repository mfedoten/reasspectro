clear all; close all; clc;
rng('default');
addpath('..');

% construct a test signal
N  = 1024;          % number of points
T  = 3;             % total time
fs = N/T;           % sampling rate
dt = 1/fs;          % sampling step
t  = 0:dt:T-dt;     % time vector
x  = chirp(t,0.1*fs,t(end),0.4*fs,'quadratic')';
x  = x + 0.5*randn(size(x));

% spectral parameters
df    = 5;                    % desired frequency resolution
Nw    = floor(fs/df);         % window/slepian sequences length
ovlap  = floor(0.9*Nw);       % overlap between adjacent windows
nfft  = 2^nextpow2(Nw);       % number of FT points
win    = gausswin(Nw,2.5);    % window for spectrogram
tapers = [3.5, 6];            % tapers for multitaper NW = 3.5; K = 6

% set additional parameters:
% - return PSD estimate;
% - padd with zeros on both sides of the signal;
% - mean averaging of tapers;
% - interpolate reassigned points for spectrogram (to cover the "holes").
optsMulti = struct('psd',1,'pad','zeros','mean','mean','interp',1);
optsSpec  = struct('psd',1,'pad','zeros','interp',1);

% time-freuqency representations
[RSmulti, fRMulti, tRMulti, Smulti, fMulti, tMulti] = reasmultitapers(x,...
    Nw,tapers,ovlap,nfft,fs,optsMulti);
[RSpec, fRSpec, tRSpec, Spec, fSpec, tSpec] = reasspecgram(x,...
    win,ovlap,nfft,fs,optsSpec);

% plot
figure; 
hf = gcf;
hf.Units = 'centimeters';
hf.Position(3:4) = [2.5*hf.Position(3) 2*hf.Position(4)];
ha1 = subplot(221); ha2 = subplot(222); ha3 = subplot(223); ha4 = subplot(224); 
set([ha1,ha2,ha3,ha4],'Units','centimeters');

% plot non-linear power and linear frequencies using 'image' function, plot COI,
% increase default font size and plot in each subplot
optsPlot = struct('type','image','dbFreq',0,'dbPow',1,'Nw',Nw,'tReal',t,...
    'font',12,'hax',ha1);
plot_spectro(tSpec,fSpec,Spec,optsPlot);
title('Conventional spectrogram','FontSize',15);

optsPlot.hax = ha2;
plot_spectro(tMulti,fMulti,Smulti,optsPlot);
title('Multitapers spectrogram','FontSize',15);

optsPlot.hax = ha3;
plot_spectro(tSpec,fSpec,RSpec,optsPlot);
title('Reassigned spectrogram','FontSize',15);

optsPlot.hax = ha4;
plot_spectro(tRMulti,fRMulti,RSmulti,optsPlot);
title('Multitapers reassigned spectrogram','FontSize',15);