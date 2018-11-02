function [colormap4imag] = colormap4imagSH(N)
if nargin < 1; N = 256; end

if mod(N,2)~=0; N = 256; disp('error(colormap4complexSH): N is an even number.'); end

sl = N/2; % section length

%   -             0　   　 　 　　　+
%  Green         Glay            Yelow
% [13,255,0] 〜 [128,128,128] 〜 [255,241,0]
scaleR = [linspace( 13,128,sl),linspace(128,255,sl)]';
scaleG = [linspace(255,128,sl),linspace(128,241,sl)]';
scaleB = [linspace(  0,128,sl),linspace(128,  0,sl)]';
colormap4imag = [scaleR, scaleG, scaleB] / 255;

end
