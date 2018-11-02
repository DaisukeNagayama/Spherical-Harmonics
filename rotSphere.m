% rotate sphere
% 
%    field      : 2dArray(É∆,É”) 
%    elevation  : angle from x-y plane [-pi/2 ~ pi/2]
%    azimuth    : angle in x-y plane   [-pi ~ pi]
%    option     : {'radian', 'angle'}
%    
%    ( theta    : angle from z axis [0 ~ pi] )
%    ( phi      : angle in x-y plane [0 ~ 2*pi] )
%    
%    field_r    : rotated field, 2dArray(É∆,É”) 
%    THETA,PHI  : 
%    [Xq,Yq,Zq] : 
%    
%    [ Example ]
%      plotSH(rotSphere(SH(1,-1), 20, 50, 'angle'));
%      img = SH(1,-1); figure, plotSH(img); figure, plotSH(rotSphere(img, 20, 50));
%    
%    [ Example ]
%      [Ymn,THETA,PHI,Xm,Ym,Zm] = spharm(3,1,[55 55],0);
%      plotSphere(rotSphere(Ymn, pi/4, -pi/4),'all');
%    

function [field_r,THETA,PHI,Xq,Yq,Zq] = rotSphere(field, elevation, azimuth, option)
OPTION_DEFAULT = 'radian';
OPTION_LIST = {'radian', 'angle'};

if (nargin<4 || isempty(option)); option = OPTION_DEFAULT; end
if not(any(strcmp(option,  OPTION_LIST))); option = OPTION_DEFAULT; disp('error(plotSphere): option = {''radian'', ''angle''} or []'); end

[NT,NP] = size(field);
THETA = linspace(0,pi,NT);  PHI = linspace(0,2*pi,NP);
[PHI,THETA] = meshgrid(PHI,THETA);

if strcmp(option,'radian')
    ele = elevation;
    azi = azimuth;
elseif strcmp(option,'angle')
    ele = elevation * pi/180;
    azi = azimuth   * pi/180;
end

[X,Y,Z] = sph2cart(PHI, THETA-pi/2, abs(field));    % orbit
% [X,Y,Z] = sph2cart(PHI, THETA-pi/2, ones(size(field)));   % sphere

% rotating [X,Y,Z]
% parallell for each [X,Y,Z]
[Xq,Yq,Zq] = arrayfun(@(x,y,z) rotatingCartCoordinate(x,y,z,ele,azi), X,Y,Z);

% interpolation the value of Vq
% parallel for each [Xq,Yq,Zq]
field_r = zeros(size(field));
field_r = arrayfun(@(xq,yq,zq) interpolationOfRotatedSphere(field,xq,yq,zq), Xq,Yq,Zq);

end

% X  : scalar
% Y  : scalar
% Z  : scalar
% ele: constant scalar
% azi: constant scalar
function [Xq,Yq,Zq] = rotatingCartCoordinate(X,Y,Z,ele,azi)
point = [X, Y, Z, 1];
getAffineF = getAffineFuncHandles;

% rotate at the azimuth
point = point * getAffineF.RotationZ(azi,0,0,0);

% rotate at the elevation
point = point * getAffineF.RotationX(ele,0,0,0);

Xq = point(1);  Yq = point(2);  Zq = point(3);
end

% V : constant 2dArray(É∆,É”)
% Xq: scalar
% Yq: scalar
% Zq: scalar
% Vq: scalar
function Vq = interpolationOfRotatedSphere(V,Xq,Yq,Zq)
[NT,NP] = size(V);

[azi,ele,~] = cart2sph(Xq,Yq,Zq);
theta = pi/2 - ele;
phi = mod(azi,2*pi);

index_theta = theta * (NT-1)/pi + 1;    % range [1 ~ NT]
index_phi   = phi * (NP-1)/(2*pi) + 1;  % range [1 ~ NP]

T1 = floor(index_theta);  T2 = mod(T1,NT-1)+1;
P1 = floor(index_phi);    P2 = mod(P1,NP-1)+1;

dT1 = index_theta - floor(index_theta); dT2 = 1 - dT1;
dP1 = index_phi   - floor(index_phi);   dP2 = 1 - dP1;

% linear interpolation
Vq = dT2 * dP2 * V(T1, P1) ...
    + dT2 * dP1 * V(T1, P2) ...
    + dT1 * dP2 * V(T2, P1) ...
    + dT1 * dP1 * V(T2, P2);
end