function [twin,dwin] = reassignment_get_windows(window,fs)
% This function constructs additional windows for spectrogram reassignment.
% 
% INPUT:  window - original window used in spectrogram;
%         fs - samplind rate.
% OUTPUT: twin = t*window(t)     - window multiplied by time vector
%         dwin = d(window(t))/dt - window time derivative
%
%
% (C) Mariia Fedotenkova 2015. 


% ---------------------------- arguments check ----------------------------
% display help if no input provided
if nargin < 1
    help reassignment_get_windows
    return
else
    % force window to be a column vector
    window = window(:);
end
% if sampling frequency is not defined set it to Fs = 1Hz
if nargin < 2 || isempty(fs)
    fs = 1;
end
% ------------------------- end of argument check -------------------------


% construct time vector
Nw = length(window);
if rem(Nw,2)
    tw = (-(Nw-1)/2:(Nw-1)/2)'/fs;  % window is odd
else
    tw = (-Nw/2:Nw/2-1)'/fs;        % window is even
end

twin = window.*tw;                  % window multipied by time
dwin = gradient(window,2*pi/fs);    % approximation of a window derivative
% gradient compures derivative using central difference and single-sided 
% differences for end points. Assuming periodicity of the window we can
% calculate end points also with central differnces
dwin(1) = (window(2) - window(end))/2*fs/2/pi;
dwin(end) = (window(1) - window(end-1))/2*fs/pi/2;

