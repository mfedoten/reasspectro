function x_padded = reassignment_pad_signal(x,Nw,type)
% A function to pad signal from left and right sides.
%
% INPUT:
% x    : original singal;
% Nw   : analysis window length;
% type : type of padding, possible values are:
%        'zeros','const','periodic','symmetric'.
%
% OUTPUT:
% x_paddede : padded signal. The length of the output's signal is now:
%             length(x)+Nw-1
%
% (C) Mariia Fedotenkoava 2016.

% determine if signal is a row vector, if so turn into column
if isrow(x)
    x_padded = x(:);
else
    x_padded = x;
end
% find amount of padding from the left and right sides, it should sum up to
% the length of the window - 1
Nleft  = floor(Nw/2);
Nright = floor((Nw-1)/2);
% chose among different padding types
switch type
    case 'zeros'
        xleft  = zeros(Nleft,1);
        xright = zeros(Nright,1);
    case 'const'
        xleft  = x_padded(1)*ones(Nleft,1);
        xright = x_padded(end)*ones(Nright,1);
    case 'periodic'
        xleft  = x_padded(end-(Nleft-1):end);
        xright = x_padded(1:Nright);
    case 'symmetric'
        xleft  = flipud(x_padded(2:Nleft+1));
        xright = flipud(x_padded(end-Nright:end-1));
    otherwise
        error('%s is not a valid padding time',type);
end
% add padding from left and right
x_padded = [xleft; x_padded; xright];
% if original vector was a row vector, turn it back to a row
if isrow(x)
    x_padded = x_padded.';
end

end

