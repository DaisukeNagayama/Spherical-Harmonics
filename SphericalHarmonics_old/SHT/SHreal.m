%%%%%%%%%%%%
% 正規化 実数 球面調和関数 Ylm を 計算する
% 
% 実数球面調和関数は 立方調和関数 や Tesseral spherical harmonics とも呼ぶらしい
% 
% 

function Ylm = SHreal(l,m)

% 球面グリッドの交点の数
polarNum = 41;      % 北極から南極 (0〜π)
azimuthNum = 41;    % 一周 (0〜2π)

% 球面グリッドの作成
theta_ = 0 : pi/(polarNum - 1) : pi;     % polar angle
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi; % azimuth angle
[phi,theta] = meshgrid(phi_,theta_);      % define the grid

% ルジャンドル陪関数の用意
Plm = legendre(l, cos(theta(:,1))); % 各行に Pl0 から Pll までが並んだ行列
Plm = Plm(abs(m) + 1,:)';           % Plm の行を取り出して縦にする
Plm = repmat(Plm,1,size(phi,2));    % 球面グリッドのサイズに合わせる

% 符号と正規化
sign = (-1)^((m+abs(m))/2);
normSH = sqrt( 4*pi / (2*l+1) * factorial(l+abs(m)) / factorial(l-abs(m)) );

% 実数 球面調和関数
if m > 0
    Ylm(:,:) = sqrt(2) * sign/normSH * Plm .* cos(abs(m) * phi);
elseif m < 0 
    Ylm(:,:) = sqrt(2) * sign/normSH * Plm .* sin(abs(m) * phi);
else
    Ylm(:,:) = sign/normSH * Plm;
end

end
