% plot spherical field
% 
% [f,s] = plotSphere(field,plotmode,plotpart,option)
% 
% [f,s] = figure handle and surf handle
% 
% field(theta,phi) = [2d array]
% 
% PLOTPART_DEFAULT = 'abs';
% PLOTMODE_DEFAULT = 'orbit';
% OPTION_DEFAULT   = 'none';
% 
% PLOTPART_LIST = {'abs', 'real', 'imag', 'all'};
% PLOTMODE_LIST = {'orbit', 'sphere'};
% OPTION_LIST   = {'none', 'gridOff'};
% 

function [f,s] = plotSphere(field,plotpart,plotmode,option)
PLOTPART_DEFAULT = 'abs';
PLOTMODE_DEFAULT = 'orbit';
OPTION_DEFAULT   = 'none';

PLOTPART_LIST = {'abs', 'real', 'imag', 'all'};
PLOTMODE_LIST = {'orbit', 'sphere'};
OPTION_LIST   = {'none', 'gridOff'};

if (nargin<2 || isempty(plotpart)); plotpart = PLOTPART_DEFAULT; end
if (nargin<3 || isempty(plotmode)); plotmode = PLOTMODE_DEFAULT; end
if (nargin<4 || isempty(option));   option   = OPTION_DEFAULT; end

if not(any(strcmp(plotpart,PLOTPART_LIST))); plotpart = PLOTPART_DEFAULT; disp('error(plotSphere): plotpart = {''abs'', ''real'', ''imag'', ''all''} or []'); end
if not(any(strcmp(plotmode,PLOTMODE_LIST))); plotmode = PLOTMODE_DEFAULT; disp('error(plotSphere): plotmode = {''sphere'', ''orbit''} or []'); end
if not(any(strcmp(option,  OPTION_LIST)));   option   = OPTION_DEFAULT;   disp('error(plotSphere): option = {''none'', ''gridOff''} or []'); end

if strcmp(plotpart,'real'); field = real(field); end
if strcmp(plotpart,'imag'); field = imag(field); end
if strcmp(plotpart,'all');  plotSphereAll(field,plotmode,option); return;  end

FLAG_COMPLEX = not(isreal(field));

% phi, theta
theta = linspace(0,  pi,size(field,1));
phi   = linspace(0,2*pi,size(field,2));
[phi,theta] = meshgrid(phi,theta);

% amplitude
amplitude = abs(field);

% sph2cart
if strcmp(plotmode,'orbit')
    [x,y,z] = sph2cart(phi, theta-pi/2, amplitude);
elseif strcmp(plotmode,'sphere')
    [x,y,z] = sph2cart(phi, theta-pi/2, ones(size(field)));
end

% color
if FLAG_COMPLEX
    C = angle(field);
else
    C = field;
end

% colormap
if FLAG_COMPLEX
    myColormap = colormap4complexSH(256);
elseif strcmp(plotpart,'imag')
    myColormap = colormap4imagSH(256);
else
    myColormap = colormap4realSH(256);
end

% alpha data
alpha = amplitude/max(max(amplitude));


% plot
scale = max(max(amplitude));
f=figure;
    s = surf(x,y,z,C); 
    if strcmp(plotmode,'orbit'); s.AlphaData = alpha; s.FaceAlpha = 'flat'; end
    if strcmp(option,'gridOff'); s.EdgeColor = 'none'; end
    colormap(myColormap)
    light
    lighting gouraud
    axis equal off
    box off
    % axis([-1 1 -1 1 -1 1]/scale)
    camzoom(1/(scale * 2.0))
    setfig
end
