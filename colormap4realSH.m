function [colormap4real] = colormap4realSH(N)
if nargin < 1; N = 256; end

if mod(N,2)~=0; N = 256; disp('error(colormap4complexSH): N is an even number.'); end

sl = N/2; % section length

%   -             0　   　 　 　　　+
%  Blue          Glay             Red
% [0,13,255] ～ [128,128,128] ～ [255,0,13]
scaleR = [linspace(  0,128,sl),linspace(128,255,sl)]';
scaleG = [linspace( 13,128,sl),linspace(128,  0,sl)]';
scaleB = [linspace(255,128,sl),linspace(128, 13,sl)]';
colormap4real = [scaleR, scaleG, scaleB] / 255;

end