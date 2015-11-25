clear all; close all; %clc;
scrsz = get(0,'ScreenSize');
set(0,'DefaultTextInterpreter','tex');
rng('default');

% construct a test signal
N  = 2560;
T  = 2560;
fs = N/T;
dt = 1/fs;
t  = 0:dt:T-dt;
x  = chirp(t,0,t(end),fs/2,'quadratic')';
x  = x + 0.5*randn(size(x));
% x = sin(2*pi*.10*t)';

% choose window
wLen  = ceil(N/10);
win   = blackman(wLen,'periodic');
ovlap = wLen-1;
% ovlap = ceil(.5*wLen);
nfft  = 2^nextpow2(2*wLen);      % twice zero-padding
% x2 = [zeros(floor(wLen/2),1);x;zeros(floor(wLen/2),1)]; 
x2 = x;
% get spectrogram
tic
[S,F,T,P] = spectrogram(x2,win,ovlap,nfft,fs);
toc

s = struct('psd',1,'crop',99.8,'pad','symmetric');%,'step',[0.008 1],'interp',false);

% get reassigned spectrogram
tic
[RS,fnew,tnew,S2,F2,T2] = reasspecgram(x2,win,ovlap,nfft,fs,s);
toc


% plot
figure('Position',[.05*scrsz(3) .2*scrsz(4) .9*scrsz(3) .6*scrsz(4)]);

subplot(131);
imagesc(T,F,abs(P));
set(gca,'YDir','normal');
title('Spectrogram MATLAB','FontSize',15);
axis square;
colorbar;

subplot(132);
imagesc(tnew,fnew,RS);
set(gca,'YDir','normal');
title('Reassigned pectrogram','FontSize',15);
axis square;
colorbar;

subplot(133);
imagesc(T2,F2,S2);
set(gca,'YDir','normal');
title('Spectrogram','FontSize',15);
axis square;
colorbar;

% % reassigned spectrogram from TFTB
% tic
% [S3,RS3,hat] = tfrrsp(x,1:length(x),nfft,blackman(wLen-1));
% toc
% 
% % plot
% figure('Position',[.05*scrsz(3) .2*scrsz(4) .9*scrsz(3) .6*scrsz(4)]);
% 
% subplot(121);
% imagesc(T,F,S3(1:end/2+1,:));
% set(gca,'YDir','normal');
% title('Spectrogram TFTB','FontSize',15);
% axis tight;
% colorbar;
% 
% subplot(122);
% imagesc(T,F,RS3(1:end/2+1,:));
% set(gca,'YDir','normal');
% title('Reassigned pectrogram TFTB','FontSize',15);
% axis tight;
% colorbar;






% % to obtain spectrograms with sliding window and return the same amount of
% % time-points
% sig = [zeros(floor(wLen/2),1);x;zeros(floor(wLen/2),1)];
% 
% % get spectrogram
% tic
% [S,F,T,P] = spectrogram(sig,win/norm(win),ovlap,nfft,fs);
% toc
% 
% % get reassigned spectrogram
% tic
% [RS,Fhat,That] = reasspecgram(sig,win,ovlap,nfft,fs);
% toc
% 
% % reassigned spectrogram from TFTB
% tic
% [S2,RS2,hat] = tfrrsp(x,1:length(x),nfft,win);
% toc




% % plot
% figure('Position',[.05*scrsz(3) .2*scrsz(4) .9*scrsz(3) .6*scrsz(4)]);
% 
% subplot(121);
% imagesc(T,F,S2(1:end/2+1,:));
% set(gca,'YDir','normal');
% title('Spectrogram TFTB','FontSize',15);
% axis tight;
% colorbar;
% 
% subplot(122);
% imagesc(T,F,RS2(1:end/2+1,:));
% set(gca,'YDir','normal');
% title('Reassigned pectrogram TFTB','FontSize',15);
% axis tight;
% colorbar;
