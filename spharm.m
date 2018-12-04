% This function generates the Spherical Harmonics basis functions of degree
% L and order M.
%
% SYNTAX: [Ymn,THETA,PHI,X,Y,Z]=spharm4(L,M,RES,PLOT_FLAG);
%
% INPUTS:
%
% L         - Spherical harmonic degree, [1x1]
% M         - Spherical harmonic order,  [1x1]
% RES       - Vector of # of points to use [#Theta x #Phi points],[1x2] or [2x1] 
% PLOT_FLAG - Binary flag to generates a figure of the spherical harmonic surfaces (DEFAULT=1)
%
%
% OUTPUTS:
% 
% Ymn   - Spherical harmonics coordinates, [RES(1) x RES(2)]
% THETA - Circumferential coordinates, [RES(1) x RES(2)]
% PHI   - Latitudinal coordinates, [RES(1) x RES(2)]
% X,Y,Z - Cartesian coordinates of magnitude, squared, spherical harmonic surface points, [RES(1) x RES(2)]
%
%
% References, 'spharm4.m' by Daniel Ennis
% 
% Nagayama Daisuke 2018/09/12

function [Ymn,THETA,PHI,Xm,Ym,Zm]=spharm(L,M,RES,PLOT_FLAG)

% Define constants (REQUIRED THAT L(DEGREE)>=M(ORDER))
if nargin==0
  L=3;   % DEGREE
  M=2;   % ORDER
end

if nargin<3
  RES=[55 55];  
  PLOT_FLAG=1;
end

if nargin<4
  PLOT_FLAG=1;
end

if L<M, error('The ORDER (M) must be less than or eqaul to the DEGREE(L).'); end

% theta, phi
THETA = linspace(0,  pi,RES(1));  % Azimuthal/Longitude/Circumferential
PHI   = linspace(0,2*pi,RES(2));  % Altitude /Latitude /Elevation
[PHI,THETA] = meshgrid(PHI,THETA);

% calculation of spherical harmonics
Lmn=legendre(L,cos(THETA));
if L~=0
    Lmn=squeeze(Lmn(abs(M)+1,:,:));
end

a1=((2*L+1)/(4*pi));
a2=factorial(L-abs(M))/factorial(L+abs(M));
C=sqrt(a1*a2);
sign = (-1)^((M+abs(M))/2);

Ymn=sign*C*Lmn.*exp(1i*M*PHI);

[Xm,Ym,Zm]=sph2cart(PHI,THETA-pi/2,abs(Ymn).^2);
[Xr,Yr,Zr]=sph2cart(PHI,THETA-pi/2,real(Ymn).^2);
[Xi,Yi,Zi]=sph2cart(PHI,THETA-pi/2,imag(Ymn).^2);
% [Xp,Yp,Zp]=sph2cart(PHI,THETA-pi/2,angle(Ymn).^2);

% color
Cm = angle(Ymn);
Cr = real(Ymn);
Ci = imag(Ymn);

% colormap
myColormap4Complex = colormap4complexSH(256);
myColormap4real = colormap4realSH(256);
myColormap4imag = colormap4imagSH(256);

% plot
if PLOT_FLAG
f=figure; axis off; hold on;
    axes('position',[0.0500 0 0.2666 1]); 
    surf(Xm,Ym,Zm,Cm);
    axis equal off; %rot3d;
    colormap(gca,myColormap4Complex);
    light; lighting phong; camzoom(1.3);
   
    axes('position',[0.3666 0 0.2666 1]); 
    surf(Xr,Yr,Zr,Cr); 
    axis equal off; %rot3d;
    colormap(gca,myColormap4real);
    light; lighting phong; camzoom(1.3);
    
    axes('position',[0.6833 0 0.2666 1]); 
    surf(Xi,Yi,Zi,Ci); 
    axis equal off; %rot3d;
    colormap(gca,myColormap4imag)
    light; lighting phong; camzoom(1.3);
  
    axes('position',[0 0.9 1 0.1]); axis off;
    t(1)=text(0.50,0.25,'Spherical Harmonics','HorizontalAlignment','Center');
    
    axes('position',[0 0 1 0.2]); axis off;
    t(2)=text(0.20,0.25,['|Y^{',num2str(M),'}_',num2str(L),'|^2'],'HorizontalAlignment','Center');
    t(3)=text(0.50,0.25,['Real(Y^{',num2str(M),'}_',num2str(L),')^2'],'HorizontalAlignment','Center');
    t(4)=text(0.80,0.25,['Imag(Y^{',num2str(M),'}_',num2str(L),')^2'],'HorizontalAlignment','Center');
    
    setfig;
end

return