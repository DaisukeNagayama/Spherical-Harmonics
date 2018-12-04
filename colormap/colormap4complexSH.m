function [colormap4complex] = colormap4complexSH(N)
if nargin < 1; N = 256; end

if mod(N,4)~=0; N = 256; disp('error(colormap4complexSH): N is multiple of 4.'); end

sl = N/4; % section length

%   0             pi/2              pi          3pi/2          2pi      
%  Blue           Green            Red         Yelow          Blue                  
% [0,13,255] ` [13,255,0] ` [255,0,13] ` [255,241,0] ` [0,13,255]
%   0              64             128           192           256
scaleR = [linspace(  0, 13,sl+1)',linspace( 13,255,sl+1)',linspace(255,255,sl+1)',linspace(255,  0,sl+1)'];
scaleG = [linspace( 13,255,sl+1)',linspace(255,  0,sl+1)',linspace(  0,241,sl+1)',linspace(241, 13,sl+1)'];
scaleB = [linspace(255,  0,sl+1)',linspace(  0, 13,sl+1)',linspace( 13,  0,sl+1)',linspace(  0,255,sl+1)'];
% scaleR = reshape(scaleR(1:sl,:),:,1);
% scaleG = reshape(scaleG(1:sl,:),:,1);
% scaleB = reshape(scaleB(1:sl,:),:,1);
% first aid: reshape error
    scaleR = scaleR(1:sl,:);
    scaleR = [scaleR(:,1);scaleR(:,2);scaleR(:,3);scaleR(:,4)];
    scaleG = scaleG(1:sl,:);
    scaleG = [scaleG(:,1);scaleG(:,2);scaleG(:,3);scaleG(:,4)];
    scaleB = scaleB(1:sl,:);
    scaleB = [scaleB(:,1);scaleB(:,2);scaleB(:,3);scaleB(:,4)];
    
colormap4complex = [scaleR, scaleG, scaleB] / (N-1);

end
