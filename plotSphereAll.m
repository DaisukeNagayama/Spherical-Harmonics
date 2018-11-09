% plot spherical field with real part and imaginary part
% 
% [f,s] = plotSphereAll(field,plotmode,option)
% 
% f = figure handle
% s = structure of surf handles
%       s.abs: abs part surf
%       s.real: real part surf
%       s.imag: imag part surf
% 
% field(theta,phi) = [2d array]
% 
% PLOTMODE_LIST = {'orbit'(default), 'sphere'};
% OPTION_LIST   = {'none'(default), 'gridoff'};
% 
% % if Input Argument = [], set the default. 
% 
% Nagayama Daisuke, 2018

function [f,s] = plotSphereAll(field,plotmode,option)
PLOTMODE_DEFAULT = 'orbit';
OPTION_DEFAULT   = 'none';

PLOTMODE_LIST = {'orbit', 'sphere'};
OPTION_LIST   = {'none', 'gridoff'};

if (nargin<2 || isempty(plotmode)); plotmode = PLOTMODE_DEFAULT; end
if (nargin<3 || isempty(option));   option   = OPTION_DEFAULT; end

if not(any(strcmp(plotmode,PLOTMODE_LIST))); plotmode = PLOTMODE_DEFAULT; disp('error(plotSphereAll): plotmode = {''sphere'', ''orbit''} or []'); end
if not(any(strcmp(option,  OPTION_LIST)));   option   = OPTION_DEFAULT;   disp('error(plotSphereAll): option = {''none'', ''gridOff''} or []'); end

% phi, theta
theta = linspace(0,  pi,size(field,1));
phi   = linspace(0,2*pi,size(field,2));
[phi,theta] = meshgrid(phi,theta);

% amplitude
amplitude_m = abs(field);
amplitude_r = abs(real(field));
amplitude_i = abs(imag(field));

% sph2cart
[Xm,Ym,Zm] = sph2cart(phi, theta-pi/2, amplitude_m);
[Xr,Yr,Zr] = sph2cart(phi, theta-pi/2, amplitude_r);
[Xi,Yi,Zi] = sph2cart(phi, theta-pi/2, amplitude_i);

% color
Cm = angle(field);
Cr = real(field);
Ci = imag(field);

% colormap
myColormap4Complex = colormap4complexSH(256);
myColormap4real = colormap4realSH(256);
myColormap4imag = colormap4imagSH(256);

% alpha
alpha_m = amplitude_m/max(max(amplitude_m));
alpha_r = amplitude_r/max(max(amplitude_r));
alpha_i = amplitude_i/max(max(amplitude_i));

% plot
f=figure; axis off; hold on;
    axes('position',[0.0500 0 0.2666 1]); 
    s1 = surf(Xm,Ym,Zm,Cm);
    if strcmp(plotmode,'orbit'); s1.AlphaData = alpha_m; s1.FaceAlpha = 'flat'; end
    if strcmp(option,'gridoff'); s1.EdgeColor = 'none'; end
    axis equal off; %rot3d;
    colormap(gca,myColormap4Complex);
    light; lighting phong; camzoom(1.3);
    
    axes('position',[0.3666 0 0.2666 1]); 
    s2 = surf(Xr,Yr,Zr,Cr);
    if strcmp(plotmode,'orbit'); s2.AlphaData = alpha_r; s2.FaceAlpha = 'flat'; end
    if strcmp(option,'gridoff'); s2.EdgeColor = 'none'; end
    axis equal off; %rot3d;
    colormap(gca,myColormap4real);
    light; lighting phong; camzoom(1.3);
    
    axes('position',[0.6833 0 0.2666 1]); 
    s3 = surf(Xi,Yi,Zi,Ci);
    if strcmp(plotmode,'orbit'); s3.AlphaData = alpha_i; s3.FaceAlpha = 'flat'; end
    if strcmp(option,'gridoff'); s3.EdgeColor = 'none'; end
    axis equal off; %rot3d;
    colormap(gca,myColormap4imag)
    light; lighting phong; camzoom(1.3);

    setfig

s(1:3) = struct('surf_handle',NaN);
s(1).surf_handle = s1;
s(2).surf_handle = s2;
s(3).surf_handle = s3;

end



